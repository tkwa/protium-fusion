"""Cost envelope for an inner membrane, using the frozen open-loop station.

These are sensitivity calculations. No row certifies a surviving membrane or a
shorter clearing time. All removed usable heat remains in the station's thermal
budget; no cooling credit is granted for radiation into the closed cavity.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path

P = Path(__file__).resolve().parent
OLD = P.parent / 'open-loop-redesign'
spec = importlib.util.spec_from_file_location('frozen_station', OLD/'station.py')
st = importlib.util.module_from_spec(spec)
spec.loader.exec_module(st)
INPUTS = json.loads((OLD/'station-inputs.json').read_text())
CYCLE, RECEIVER = INPUTS['cycle'], INPUTS['receiver']
R = 1.1e9
A = 4*math.pi*R**2
F = .8/(RECEIVER['occupied_s_per_He_at_R1e9']['mass999']*R/1e9)
YEAR = st.YEAR


def evaluate(*, added_kg_m2=0., auxiliary_W_m2=0., lost_gross_fraction=0.,
             renewed_kg_m2_year=0., renewal_J_kg=1e8, rate_fraction=1.,
             renewed_carbon_fraction=0., carbon_return_fraction=1.,
             other_material_return_fraction=1.):
    """Retain original conversion hardware and add declared penalties.

    lost_gross_fraction: fraction of baseline gross electric conversion lost
    through lower-quality heat or an energy-routing change. Retained waste
    heat is still rejected. Added hardware receives the original annual
    1e-4 refurbishment allowance; renewed_kg_m2_year is EXTRA throughput above
    that allowance. Fabrication energy applies to all renewed mass, including
    internally recycled stock; feed applies only to irreversibly lost material.
    """
    assert added_kg_m2 >= 0 and auxiliary_W_m2 >= 0 and renewed_kg_m2_year >= 0
    assert 0 <= lost_gross_fraction <= 1 and 0 < rate_fraction <= 1
    assert 0 <= renewed_carbon_fraction <= 1
    assert 0 <= carbon_return_fraction <= 1 and 0 <= other_material_return_fraction <= 1
    b = st.evaluate(F*rate_fraction,R,CYCLE,RECEIVER,absorber_sigma=5000.,
                    driver_peak_W_kg=100.,service_sigma=1000.)
    rep = 1e-4*1e8/YEAR
    flux = .9*st.SB*450**4
    support_coeff = .001*st.G*1e4/(8*math.pi*R**2)
    stream_coeff = st.G/(2*R*1e8)
    tide = b['tidal_equivalent_mass_kg']
    linear = .2+.001+rep+20*(rep+support_coeff*tide)/flux+stream_coeff*tide
    quadratic = stream_coeff+20*support_coeff/flux
    mb = b['mass_breakdown_kg']
    feedback_keys = {'radiator_panels','replacement_manufacturing','ten_year_spares','support_streams','installed_reserve'}
    fixed = sum(v for k,v in mb.items() if k not in feedback_keys)
    before_replacement = b['net_W']+b['replacement_W']+b['support_W']
    energy_minus_pre_export = b['waste_W']+b['net_W']-before_replacement
    constant0 = fixed+20*(energy_minus_pre_export+b['absorbed_solar_W'])/flux
    renewal_kg_s = renewed_kg_m2_year*A/YEAR
    renewal_W = renewal_kg_s*renewal_J_kg
    other_penalty = auxiliary_W_m2*A + lost_gross_fraction*b['gross_W']
    # 1 kg/W additional renewal manufacturing + 20 kg/m2 radiator panels.
    constant = constant0+A*added_kg_m2+renewal_W+20*(renewal_W+other_penalty)/flux
    discr = (1-linear)**2-4*quadratic*constant
    if discr < 0:
        return {'mass_solution':False}
    mass = 2*constant/(1-linear+math.sqrt(discr))
    support = support_coeff*mass*(mass+tide)
    power = before_replacement-rep*mass-support-renewal_W-other_penalty
    waste = b['waste_W']+b['net_W']-power
    radiator = (waste+b['absorbed_solar_W'])/flux
    reconstructed_mass = fixed+A*added_kg_m2+renewal_W+20*radiator+rep*mass+.001*mass+stream_coeff*mass*(mass+tide)+.2*mass
    carbon_extra = renewal_kg_s*renewed_carbon_fraction*(1-carbon_return_fraction)
    other_extra = renewal_kg_s*(1-renewed_carbon_fraction)*(1-other_material_return_fraction)
    carbon_ceiling = b['net_material_C_allowance_kg_s']
    max_carbon_escape = carbon_ceiling/(renewal_kg_s*renewed_carbon_fraction) if renewal_kg_s*renewed_carbon_fraction else None
    return {
        'mass_solution':True, 'positive_net_power':power>0,
        'rate_fraction':rate_fraction, 'added_kg_m2':added_kg_m2,
        'auxiliary_W_m2':auxiliary_W_m2,'lost_gross_fraction':lost_gross_fraction,
        'renewed_kg_m2_year':renewed_kg_m2_year,'renewal_J_kg':renewal_J_kg,
        'mass_kg':mass,'net_W':power,'net_W_kg':power/mass,'net_W_m2':power/A,
        'net_W_kg_relative_to_baseline':power/mass/BASELINE_W_KG,
        'extra_renewal_kg_s':renewal_kg_s,'extra_renewal_W':renewal_W,
        'extra_renewal_manufacturing_kg':renewal_W,'other_penalty_W':other_penalty,
        'waste_W':waste,'radiator_area_fraction':radiator/A,
        'net_extra_C_feed_kg_s':carbon_extra,'net_extra_noncarbon_feed_kg_s':other_extra,
        'total_material_C_allowance_kg_s':carbon_ceiling,
        'carbon_allowance_used_by_extra_renewal_fraction':carbon_extra/carbon_ceiling,
        'max_unreturned_fraction_if_all_C_reserve_assigned_here':min(1,max_carbon_escape) if max_carbon_escape is not None else None,
        'mass_relative_residual':reconstructed_mass/mass-1,
        'energy_relative_residual':(power+waste-b['waste_W']-b['net_W'])/b['gross_W'],
    }

BASELINE = json.loads((OLD/'operating-cases.json').read_text())['cases']['high_source_rating_surface']
BASELINE_W_KG = BASELINE['net_W_kg']

if __name__ == '__main__':
    zero = evaluate()
    assert abs(zero['net_W_kg']/BASELINE_W_KG-1)<1e-12
    rows = {'baseline_reproduction':zero}
    for sigma in [.0022,.1,1.,10.,100.,1000.]:
        rows[f'mass_only_{sigma:g}_kg_m2']=evaluate(added_kg_m2=sigma)
    for q in [.1,1.,5.,10.,20.]:
        rows[f'auxiliary_only_{q:g}_W_m2']=evaluate(auxiliary_W_m2=q)
    for loss in [.01,.1,.25,.5]:
        rows[f'conversion_loss_only_{loss:g}']=evaluate(lost_gross_fraction=loss)
    for fraction in [.5,.1]:
        rows[f'cadence_fraction_{fraction:g}']=evaluate(rate_fraction=fraction)
    intervals=BASELINE['pulse_intervals_s']
    frequencies={
        'C13_only':YEAR/intervals['c13'],
        'large_pulses':sum(YEAR/intervals[k] for k in ['c12','c13','n14']),
        'all_pulses':sum(YEAR/t for t in intervals.values()),
    }
    for material,density,carbon in [('silica',2200.,0.),('carbon',2000.,1.)]:
      for thickness in [1e-8,1e-6]:
       for cadence,n in frequencies.items():
        sigma=density*thickness
        rows[f'renew_{material}_{thickness:g}_m_{cadence}']=evaluate(
            added_kg_m2=sigma,renewed_kg_m2_year=sigma*n,
            renewed_carbon_fraction=carbon,carbon_return_fraction=0.,other_material_return_fraction=0.)
    for row in rows.values():
        assert abs(row['mass_relative_residual'])<1e-12
        assert abs(row['energy_relative_residual'])<1e-12
    thresholds=[]
    for sigma in [0.,10.,100.,1000.]:
        lo,hi=0.,30.
        for _ in range(80):
            q=(lo+hi)/2
            if evaluate(added_kg_m2=sigma,auxiliary_W_m2=q)['net_W_kg']>=.001:lo=q
            else:hi=q
        thresholds.append({'added_kg_m2':sigma,'max_extra_average_electric_W_m2_for_1mWkg':lo})
    # The loss budget is shared with the existing wall and other carbon parts.
    out={
        'scope':'Conditional fixed-radius station penalties; no row qualifies thermal survival, stopping, membrane life, or cadence. Carbon return requirements allocate the ENTIRE material reserve to the new membrane and therefore are necessary, not sufficient.',
        'radius_m':R,'area_m2':A,'helium_per_s':F,
        'fresh_C_H_mass_budget':.007,'fuel_C_H_mass_ratio':.0056,
        'baseline_net_W_m2':BASELINE['net_W']/A,
        'pulse_counts_per_year':frequencies,
        'one_mW_kg_penalty_thresholds':thresholds,
        'cases':rows,
        'sources_sha256':{str(f.relative_to(P.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [OLD/'station.py',OLD/'station-inputs.json',OLD/'operating-cases.json']},
        'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks':{'zero_extension_matches_frozen_station':True,'mass_energy_residuals_below_1e_12':True},
    }
    (P/'station-extension.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'baseline_W_m2':out['baseline_net_W_m2'],'counts':frequencies,'cases':{k:{x:v[x] for x in ['net_W_kg','positive_net_power','extra_renewal_W','max_unreturned_fraction_if_all_C_reserve_assigned_here']} for k,v in rows.items()}},indent=2))
