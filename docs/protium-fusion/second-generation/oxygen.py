#!/usr/bin/env python3
"""Extended-pool absorption, oxygen flows, charge and allocated-work comparisons.

Repeated passes use prescribed host-pulse baths. Pure oxygen pulse hydrodynamics
and isotope sorting hardware are still qualifications, not inferred successes.
"""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
import network as net

P=Path(__file__).resolve().parent
BASE=json.loads((P.parent/'mixture-design/selected-plant/point.json').read_text())
PL=BASE['plant'];CY=BASE['cycle'];MU=1.66053906660e-27;QE=1.602176634e-19;MEV=1.602176634e-13

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--return-probes',action='store_true');a=ap.parse_args()
    data=json.loads((P/'network-32.json').read_text()); rows={(r['stage'],r['isotope']):r for r in data['tracer']}
    states=['o16','o17','o18','f19'];policy={'o16':'n14','o17':'c12','o18':'c12','f19':'c12'}
    selected=[rows[(policy[n],n)] for n in states]
    M=np.array([[r['cold'].get(n,0.) for r in selected] for n in states])
    initial=np.array([1.,0.,0.,0.]);visits=np.linalg.solve(np.eye(4)-M,initial)
    destinations={}
    for v,r in zip(visits,selected):
        for n,c in r['cold'].items():
            if n not in states:destinations[n]=destinations.get(n,0.)+v*c
    success=sum(destinations.get(n,0.) for n in ['c12','c13','n14','n15'])
    beyond=sum(v for n,v in destinations.items() if net.Z[net.NAMES.index(n)]>=10)
    hot={k:sum(v*r['integrated_hot_ledger'][k] for v,r in zip(visits,selected)) for k in net.TRACK_NAMES}
    G,C=net.tracer_generator(net.rate_values(1e6),0.,np.zeros(net.NN+1));unst=np.flatnonzero(-np.diag(G)>1e-18)
    cold_time=np.linalg.inv(-G[np.ix_(unst,unst)])
    cold=np.zeros(len(net.TRACK_NAMES));cold_hold=0.
    for v,r in zip(visits,selected):
        y=np.array([r['hot'].get(net.NAMES[n],0.) for n in net.HEAVY]);residence=cold_time@y[unst]
        cold+=v*(C[len(net.WI):,unst]@residence)
        cold_hold+=v*residence.sum()
    cold=dict(zip(net.TRACK_NAMES,cold.tolist()))
    total={k:hot[k]+cold[k] for k in net.TRACK_NAMES}
    # Stage input energy per catalyst from the final-family account.
    family=json.loads((P.parent/'component-design/final-family-dispatch/sigma350/input.json').read_text())
    U=sum(v*family['loaded_stages'][policy[n]]['input_MeV'] for n,v in zip(states,visits))
    eta=.5457153023359447
    # All positive weak Q is an upper debit on neutrinos; zero is the other energy bound.
    Q=total['total_Q_MeV'];weak=total['weak_Q_MeV']
    E=[eta*(Q-nu+U)-U/.6 for nu in [max(0.,weak),0.]]
    Ndot=PL['helium_per_s'];source=CY['ash_atoms_per_He']['o16']*Ndot
    # Both mass and nuclear-number interpretations are explicit.
    catdot=sum(CY['stage_attempts_per_He'].values())*Ndot
    residence=PL['inventory_residence_assumptions_s']['retained_CNO']
    catinventory=PL['recovered_CNO_kg_s']*residence
    numberinventory=catdot*residence
    invO_mass=.01*catinventory; invO_number=.01*numberinventory
    cases={}
    for label,feed in [('actual_oxygen_production',source),
                       ('one_percent_catalyst_number_throughput',.01*catdot),
                       ('one_percent_catalyst_mass_throughput',.01*PL['recovered_CNO_kg_s']/(16*MU))]:
        stage_attempts={s:0. for s in ['n14','n15','c12','c13']}
        for n,v in zip(states,visits):stage_attempts[policy[n]]+=v*feed
        attempts=sum(stage_attempts.values());large_slots=0.;small_slots=0.
        for stage,flux in stage_attempts.items():
            duration=1.1*family['timing'][stage]['mass999']['occupied_s']
            # Illustrative pure-oxygen equivalent keeps the host H/catalyst number ratio.
            h=family['loaded_stages'][stage]['protons_per_catalyst'];target_A=16+h
            small_N=1.1309733552923254e11/(target_A*MU)
            small_slots+=flux/small_N*duration
            large_slots+=flux/(100*small_N)*duration
        cases[label]=dict(new_O16_nuclei_s=feed,new_O16_kg_s=feed*16*MU,
            repeated_heavy_nucleus_passes_s=attempts,O16_mass_equivalent_pass_kg_s=attempts*16*MU,
            host_attempts_s=stage_attempts,small_shot_serial_calendar_fraction=small_slots,
            large_shot_serial_calendar_fraction=large_slots,
            allocated_prepared_W=U*feed*MEV,allocated_driver_electric_W=U/.6*feed*MEV,
            allocated_incremental_net_W_bounds=[e*feed*MEV for e in E],
            selective_sort_100eV_at50pct_W=attempts*200*QE,
            processing_mass_kg_at1e8_kg_per_kg_s=attempts*16*MU*1e8)
    mix={};mix_sum=0.
    for s in ['n14','n15','c12','c13']:
        r=rows[(s,'o16')];q=1-r['cold'].get('o16',0.);capacity=.01*CY['stage_attempts_per_He'][s]*q
        mix[s]=dict(O16_leaves_O16_probability=q,new_O16_conversion_per_He_at1pct_admixture=capacity)
        mix_sum+=capacity
    out=dict(scope=__doc__,states=states,policy=policy,transient_matrix=M.tolist(),
        attempts_per_introduced_O16=dict(zip(states,visits.tolist())),
        final_nuclei_per_introduced_O16=destinations,quartet_return_probability=success,
        beyond_F_quarantine_fraction=beyond,hot_integrated_ledger=hot,cold_integrated_ledger=cold,
        total_reaction_ledger=total,radioactive_residence_nucleus_s_per_O16=cold_hold,
        allocated_prepared_MeV_per_O16=U,allocated_net_MeV_per_O16_bounds=E,
        replacement_carbon_shadow_energy_J_kg=[-e*MEV/(12*MU) for e in E],
        inventory=dict(catalyst_kg=catinventory,catalyst_nuclei=numberinventory,
            O_at1pct_mass_share_kg=invO_mass,O_at1pct_number_share_kg=invO_number*16*MU,
            actual_O_production_kg_s=source*16*MU,
            O1pct_mass_inventory_net_residence_years=invO_mass/(source*16*MU)/(365.25*86400),
            O1pct_number_inventory_net_residence_years=invO_number/source/(365.25*86400),
            gross_attempt_flow_if1pct_mass_throughput_kg_s=.01*PL['recovered_CNO_kg_s'],
            corresponding_net_O_conversion_kg_s=.01*PL['recovered_CNO_kg_s']/visits.sum()),
        flows=cases,mixed_host_capacity=mix,
        all_hosts_1pct_conversion_per_He=mix_sum,
        uniform_admixture_to_balance_original_O_source=.01*CY['ash_atoms_per_He']['o16']/mix_sum,
        charge=dict(small_equivalent_O_H2_shot_mass_kg=1.1309733552923254e11,
            O_H2_nuclei_per_small_shot=1.1309733552923254e11/(18*MU),
            fully_stripped_positive_charge_C=10*QE*1.1309733552923254e11/(18*MU),
            launch_grain_charge_C_at_0p002656_Ckg=.002656*1.1309733552923254e11,
            circuit='The ion charge and local electron return are equal and opposite. Grain acceleration charge is a distinct earlier delivery circuit.'),
        terminal_loss=dict(old_baseline_heavy_ash_beyond_F_per_He=sum(v for n,v in CY['ash_atoms_per_He'].items() if net.Z[net.NAMES.index(n)]>=10),
            added_recycle_quarantine_per_He=CY['ash_atoms_per_He']['o16']*beyond,
            eventual_material_loss_fraction_for_0p007_budget=.007*PL['fresh_H_kg_s']/PL['recovered_CNO_kg_s']),
        qualifications=['A tracer follows actual saved temperature, density and n/p/alpha abundance histories; bulk oxygen-only shots need their own radiation hydrodynamics.',
            'Equivalent host preparation is an allocated-work comparison. Piggybacking trades fuel displacement and altered hydro for the separate-shot debit.',
            'Cold map includes branching and long lived states. The radioactive-residence inventory exposes their storage cost.',
            'Beyond-F inventory is an operational quarantine boundary; physical return is tested separately. No assumption that every alpha capture or every temporary Ne isotope is irreversible.',
            'Frozen REACLIB oxygen/fluorine rates remain separate from newly read 2025/2026 measurements. Updated oxygen/F19 tables and enlarged-network coupled burns are acceptance gates.'])
    out['watched_reactions_per_restored_input_O16']={key:sum(v*r['watched_events'][key] for v,r in zip(visits,selected)) for key in selected[0]['watched_events']}
    bulk=json.loads((P/'bulk-16.json').read_text())
    base_quarantine=0.
    for row in bulk['cases']:
        if row['case']=='control':
            base_quarantine+=CY['stage_attempts_per_He'][row['stage']]*sum(v for name,v in row['cold_heavy_per_initial_heavy'].items() if name not in ['c12','c13','n14','n15','o16','o17','o18','f19','f16'])
    out['terminal_loss']['enlarged_network_baseline_quarantine_per_He']=base_quarantine
    out['terminal_loss']['baseline_plus_O_return_quarantine_per_He']=base_quarantine+CY['ash_atoms_per_He']['o16']*(1-success)
    mass=dict(zip(net.NAMES,net.MASS))
    out['return_threshold_MeV']=dict(Ne20_gamma_alpha_O16=mass['o16']+mass['he4']-mass['ne20'],
        Ne20_p_alpha_F17=mass['f17']+mass['he4']-mass['ne20']-mass['p'],
        Ne22_p_alpha_F19=mass['f19']+mass['he4']-mass['ne22']-mass['p'])
    probe=json.loads((P/'probe-cycles.json').read_text())
    new=next(r['cycle'] for r in probe['cases'] if r['record']=='pulse-n15-m100-v1-z128-t1.json')
    new_source=new['ash_atoms_per_He']['o16']*Ndot*16
    out['large_N15_16x_source_comparison']=dict(O16_nuclei_s=new_source,O16_kg_s=new_source*16*MU,
        equivalent_repeated_O_pass_kg_s=new_source*visits.sum()*16*MU,
        allocated_incremental_net_W_bounds=[e*new_source*MEV for e in E],
        qualification='Uses the unchanged N14/C12 return baths and the new large-N15 quartet oxygen source. Mixed-pulse hydro feedback remains open.')
    if not a.return_probes and (P/'oxygen-results.json').exists():
        previous=json.loads((P/'oxygen-results.json').read_text())
        if 'quarantined_return_probes' in previous:out['quarantined_return_probes']=previous['quarantined_return_probes']
    if a.return_probes:
        with gzip.open(P/'trajectories.json.gz','rt') as f:st=json.load(f)
        B,_=net.cold_map();probe=[]
        for stage in ['n14','c13']:
            for iso in ['ne20','ne21','ne22']:
                avg=np.zeros(net.NH)
                for i in net.sample(st[stage],8):avg+=net.tracer(st[stage],i,iso)[0]/8
                dec=B@avg
                probe.append(dict(stage=stage,isotope=iso,
                    cold_quartet_probability=float(sum(dec[net.HI[net.NAMES.index(n)]] for n in ['c12','c13','n14','n15'])),
                    cold_extended_pool_probability=float(sum(dec[k] for k,n in enumerate(net.HEAVY) if net.Z[n] in [6,7,8,9]))))
        out['quarantined_return_probes']=probe
    (P/'oxygen-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['attempts_per_introduced_O16','quartet_return_probability','beyond_F_quarantine_fraction','total_reaction_ledger','allocated_prepared_MeV_per_O16','allocated_net_MeV_per_O16_bounds','uniform_admixture_to_balance_original_O_source','inventory']},indent=2))

if __name__=='__main__':main()
