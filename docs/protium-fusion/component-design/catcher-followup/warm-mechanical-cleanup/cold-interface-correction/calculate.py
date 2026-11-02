"""Passive finite-speed cold boundary for the frozen exact-dispatch hot buffer.

This supplement imports only definitions from the archived solver, preserving
its files. Cold panels are quasi-steady; a separate manifold-volume transient
requires its own pneumatic time constant. No steady mean is substituted for
the actual time-dependent cold partial pressures.
"""
from pathlib import Path
import json, math, hashlib
import numpy as np

P = Path(__file__).resolve().parent
ROOT = Path('/home/tkwa/.codex/worktrees/protium-fusion/website')
W = ROOT/'docs/protium-fusion/component-design/catcher-followup/warm-mechanical-cleanup'
G = W.parent/'gas'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read = lambda p: json.loads(Path(p).read_text())
src = W/'buffer.py'
ns = {'__file__': str(src)}
exec(compile(src.read_text().split('outrows=[];')[0], str(src), 'exec'), ns)
old = read(W/'buffer.json')
Sc = ns['Sc'].copy()
old_floor = ns['floor'].copy()
Th, Tc, V, K = 1000., 300., .5, 1e6
joint = read(G/'joint-ports.json')
joint_ref = next(r for r in joint['scenarios']
                 if r['downstream_speed_target_m3_s_per_wall_m2_each_species']==.001
                 and r['shared_mouth_area_fraction']==.0001
                 and r['net_molecular_transmission']==.1)
areas = np.array(joint_ref['required_active_panel_area_fraction_Hgetter_Hecryo'])
panel_ref_speed = np.array([10., 13.])
panel_ref_T = np.array([300., 273.15])
Spanel = areas*panel_ref_speed*np.sqrt(Tc/panel_ref_T)
Seff = Sc/(1+Sc*np.sqrt(Tc/Th)/Spanel)
dt, points = 50., 1024
ts, q, metadata = ns['sources'](dt, points)
checks = []
def check(label, condition):
    assert bool(condition), label
    checks.append(label)

# Reproduce the archived reference before changing any boundary.
base, base_total, base_b, base_s = ns['response'](q, dt, V, K)
for field in ['total_buffer_pressure_min_Pa', 'total_buffer_pressure_max_Pa',
              'total_buffer_pressure_exact_mean_Pa', 'cold_outlet_PV_flow_max_Pa_m_s']:
    check('replay_'+field, abs(base[field]-old['reference'][field]) < 1e-11*max(1, abs(base[field])))
check('provided_H_series_speed', abs(Seff[0]/1.4048324790768463e-6-1)<1e-12)
check('provided_He_series_speed', abs(Seff[1]/9.973280306284756e-7-1)<1e-12)

ns['Sc'] = Seff
cases = []
curves = {}
for name, Peq in [('zero_equilibrium', np.zeros(2)),
                  ('He_equilibrium_1e-6_Pa', np.array([0., 1e-6]))]:
    ns['floor'] = np.sqrt(Th/Tc)*Peq
    row, total, b, s = ns['response'](q, dt, V, K)
    hotflow = Seff[:, None]*(b-ns['floor'][:, None])
    pc = Peq[:, None]+(Tc/Th)*hotflow/Spanel[:, None]
    check(name+'_nonnegative_flow', np.min(hotflow)>0)
    check(name+'_cold_panel_flux', np.max(np.abs(Spanel[:,None]*(pc-Peq[:,None])-(Tc/Th)*hotflow))<1e-18)
    check(name+'_original_hot_conduit_flux', np.max(np.abs(Sc[:,None]*(b-np.sqrt(Th/Tc)*pc)-hotflow))<1e-18)
    # A finite cold reservoir with conductance-only response is a separate
    # possible lag. Here both sides close at the stated quasi-steady boundary.
    row['name'] = name
    row['equilibrium_cold_partial_pressure_Pa'] = Peq.tolist()
    row['source_discretization'] = metadata
    row['changes_from_frozen_reference'] = {
        f: row[f]-base[f] for f in ['total_buffer_pressure_min_Pa',
            'total_buffer_pressure_max_Pa', 'total_buffer_pressure_exact_mean_Pa',
            'total_source_pressure_discretized_max_Pa', 'cold_outlet_PV_flow_max_Pa_m_s']}
    row['changes_from_frozen_reference']['maximum_absolute_source_pressure_change_Pa'] = float(np.max(np.abs(s-base_s)))
    row['changes_from_frozen_reference']['maximum_absolute_buffer_pressure_change_Pa'] = float(np.max(np.abs(b-base_b)))
    row['cold_pressure_total_sampled_min_Pa'] = float(np.min(pc.sum(axis=0)))
    row['cold_pressure_total_sampled_max_Pa'] = float(np.max(pc.sum(axis=0)))
    row['cold_pressure_total_exact_mean_Pa'] = 0.
    row['initial_hot_buffer_gas_mass_kg_m2'] = float(sum(z['initial_buffer_pressure_Pa']*V*m/(ns['KB']*Th) for z,m in zip(row['species'],ns['m'])))
    row['initial_hot_buffer_gas_mass_kg'] = row['initial_hot_buffer_gas_mass_kg_m2']*ns['A']
    for j, z in enumerate(row['species']):
        meanflow = z['cold_out_PV_Pa_m']/metadata['duration_s']
        coldmean = float(Peq[j]+(Tc/Th)*meanflow/Spanel[j])
        reverse = np.sqrt(Th/Tc)*pc[j]/b[j]
        imin = int(np.argmin(b[j]))
        z.update({'panel_speed_at_300K_m3_s_per_wall_m2':float(Spanel[j]),
                  'hot_conduit_speed_without_panel_m3_s_per_wall_m2':float(Sc[j]),
                  'combined_hot_speed_m3_s_per_wall_m2':float(Seff[j]),
                  'finite_record_mean_hot_PV_flow_Pa_m_s':float(meanflow),
                  'cold_pressure_sampled_min_Pa':float(min(pc[j])),
                  'cold_pressure_sampled_max_Pa':float(max(pc[j])),
                  'cold_pressure_exact_mean_Pa':coldmean,
                  'cold_pressure_at_hot_buffer_min_Pa':float(pc[j,imin]),
                  'reverse_flux_fraction_min':float(min(reverse)),
                  'reverse_flux_fraction_max':float(max(reverse)),
                  'reverse_flux_fraction_at_hot_buffer_min':float(reverse[imin]),
                  'initial_buffer_pressure_change_Pa':z['initial_buffer_pressure_Pa']-base['species'][j]['initial_buffer_pressure_Pa'],
                  'initial_source_pressure_change_Pa':z['initial_source_pressure_Pa']-base['species'][j]['initial_source_pressure_Pa']})
        row['cold_pressure_total_exact_mean_Pa'] += coldmean
        if not Peq[j]:
            check(name+f'_constant_passive_reverse_fraction_{j}', np.ptp(reverse)<1e-15)
            check(name+f'_reverse_series_identity_{j}', abs(reverse[0]-(1-Seff[j]/Sc[j]))<1e-15)
    keep = np.unique(np.r_[np.arange(0,len(ts),100),np.argmax(total),np.argmin(total),np.argmax(pc.sum(axis=0))])
    curves[name] = {'time_s':((keep+1)*dt).tolist(), 'buffer_partial_pressure_Pa':b[:,keep].tolist(), 'cold_partial_pressure_Pa':pc[:,keep].tolist(), 'cold_outlet_hot_PV_flow_Pa_m_s':hotflow[:,keep].tolist()}
    cases.append(row)
    if name == 'zero_equilibrium':
        zero_b, zero_s, zero_flow = b.copy(), s.copy(), hotflow.copy()
    else:
        expected = np.sqrt(Th/Tc)*Peq[:,None]
        check('Peq_pressure_shift_exact', np.max(np.abs(b-zero_b-expected))<1e-12)
        check('Peq_source_shift_exact', np.max(np.abs(s-zero_s-expected/K))<1e-12)
        check('Peq_net_flow_invariant',np.max(np.abs(hotflow-zero_flow))<1e-18)

# Regulated pressure has an exact constant shift, but differs physically from
# the passive panel used above. Display it as a comparison only.
regulated_pc = (Tc/Th)*ns['qmean']/Spanel
shift = np.sqrt(Th/Tc)*regulated_pc-old_floor
regulated = {'cold_partial_pressure_Pa':regulated_pc.tolist(),
    'buffer_partial_pressure_shift_Pa':shift.tolist(),
    'source_partial_pressure_shift_Pa':(shift/K).tolist(),
    'total_buffer_pressure_min_Pa':base['total_buffer_pressure_min_Pa']+float(sum(shift)),
    'total_buffer_pressure_max_Pa':base['total_buffer_pressure_max_Pa']+float(sum(shift)),
    'scope':'A separately regulated constant cold pressure, included only to distinguish the exact pressure-shift limit from the passive finite-speed response.'}

out = {'title':'Passive cold-panel closure of the exact-dispatch hot buffer',
    'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [src,W/'buffer.json',G/'inputs.json',G/'pumping-interfaces.json',G/'joint-ports.json']},
    'code_sha256':sha(__file__),
    'parameters':{'hot_T_K':Th,'cold_reservoir_T_K':Tc,'buffer_volume_per_wall_area_m':V,
        'primary_booster_K0':K,'active_panel_area_per_wall_area':areas.tolist(),
        'panel_reference_speed_per_active_area_m_s':panel_ref_speed.tolist(),
        'panel_reference_T_K':panel_ref_T.tolist(), 'panel_speed_at_300K_m_s':Spanel.tolist(),
        'original_hot_conduit_speed_m_s':Sc.tolist(),'combined_passive_hot_speed_m_s':Seff.tolist(),
        'speed_ratio':(Seff/Sc).tolist(),
        'cold_pressure_relation':'Pc=Peq+(Tc/Th)*Qhot/Spanel',
        'combined_speed_relation':'Seff=Sc/[1+Sc*sqrt(Tc/Th)/Spanel]',
        'hot_floor_relation':'sqrt(Th/Tc)*Peq',
        'cold_manifold_assumption':'Negligible pneumatic storage and extra resistance for the declared reference. A separate manifold conductance >=10Spanel adds <=10% to panel-side pressure rise; a finite-volume lag needs its own transient.',
        'species_order':['H2','He']},
    'cases':cases,'regulated_comparison':regulated,
    'unchanged_hardware_heat':{'cold_conduit_hot_radiation_W_m2':old['cold_conduit']['all_admitted_hot_radiation_W_m2_upper'],
        'cold_conduit_solid_conduction_W_m2':old['cold_conduit']['solid_conduction_W_m2'],
        'scope':'The physical shared conduit and panel areas are retained. Sorber cycling, refrigerator heat and chemical separation use the gas-agent handoff; no new credit or duplicated debit is introduced here.'},
    'checks_passed':len(checks)+len(ns['checks']),'checks':checks,'inherited_solver_checks':ns['checks'],
    'limits':['Panel speeds and sqrt(T) scaling retain their existing conditional capture/loading assumptions.',
        'Pc is a reservoir-equivalent pressure at the 300K interface, not the physical temperature of the 5K helium adsorbent.',
        'Peq is an explicit equilibrium/backstream assumption, not a measured irradiated sorber property.',
        'Pressure extrema are sampled exact bin endpoints; the frozen25/50/100s controls bound source-bin sensitivity but are not rerun here.',
        'Source plenum and hot buffer storage are included; extra cold-manifold storage, valve transients and breakthrough are not included.',
        'The direct source here releases 100% of incident H/He into the local buffer and is a capacity response. The archived global-cleanup increment is an independent capacity envelope, not extra recovered mass. A physical coupled history partitions the same incident atoms among local capture, solid hold-up and cavity return.',
        'Separated N/O is an additional species stream; it does not create additional H/He inventory.']}
(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
(P/'curves.json').write_text(json.dumps(curves,separators=(',',':'),allow_nan=False)+'\n')
print(json.dumps({'checks':out['checks_passed'],'cases':[{k:r[k] for k in ['name','total_buffer_pressure_min_Pa','total_buffer_pressure_max_Pa','total_buffer_pressure_exact_mean_Pa','cold_pressure_total_exact_mean_Pa','cold_pressure_total_sampled_max_Pa','cold_outlet_peak_over_continuous_mean','initial_hot_buffer_gas_mass_kg_m2','changes_from_frozen_reference']} for r in cases]},indent=2))
