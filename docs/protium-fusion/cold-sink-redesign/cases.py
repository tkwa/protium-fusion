"""Small auditable cold-sink family; no radius or pulse-rate optimization."""
import hashlib
import importlib.util
import json
from pathlib import Path
from station import evaluate, finite_engine, direct_marginal_bound, MEV

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    bp=ROOT/'open-loop-redesign/operating-cases.json'
    ip=ROOT/'open-loop-redesign/station-inputs.json'
    hp=HERE/'heat-engine.json'
    base=json.loads(bp.read_text())['cases']['high_source_rating_surface']
    inp=json.loads(ip.read_text());cycle=inp['cycle']
    heat=json.loads(hp.read_text());heatflux=heat['baseline']['captured_heat_W_m2']
    cycle60=heat['conditional_regenerative_Brayton_cycles']['six_stage_93pct_radiator200K']
    UA60=(cycle60['UA_recuperator_W_K_m2']+cycle60['UA_generator_cooling_W_K_m2'])/heatflux
    cycle65=heat['conditional_regenerative_Brayton_cycles']['six_stage_97pct_radiator200K']
    UA65=(cycle65['UA_recuperator_W_K_m2']+cycle65['UA_generator_cooling_W_K_m2'])/heatflux
    cases={}
    def add(name,group,**kwargs):
        r=evaluate(base,cycle,**kwargs)
        r['group']=group
        cases[name]=r
        return r
    historical=dict(heat_efficiency=.4,radiator_K=450.,routine_MJ_kg=100.,
        environment_absorbed_W_m2=0.,bypass_W_m2=0.,pump_fraction_hot=0.,
        include_incremental_hx=False,hx_internal_UA_per_hot_W_K=0.)
    add('historical_450K_40pct_sun','historical_control',**(historical|{'direct_solar_absorbed_W_m2':13.61}))
    add('historical_197K_60pct_sun','historical_control',**(historical|{'direct_solar_absorbed_W_m2':13.61,'radiator_K':197.,'heat_efficiency':.6}))
    add('historical_197K_60pct_shaded','historical_control',**(historical|{'radiator_K':197.,'heat_efficiency':.6}))
    add('shading_only_450K_40pct_oldcost','isolated_effect_control',**historical)
    add('shading_plus_cold_197K_40pct_oldcost','isolated_effect_control',**(historical|{'radiator_K':197.}))
    # Working 60% is package electrical conversion including the separately
    # studied cycle's own pump and bearing allowances. Extra pumps are zero here.
    common=dict(heat_efficiency=.6,carnot_quality=.8,radiator_K=200.,receiver_K=1000.,
        routine_MJ_kg=2.,environment_absorbed_W_m2=11.0254,environment_basis='whole_shell',
        radiator_available_fraction=.9,bypass_W_m2=1.,hx_internal_UA_per_hot_W_K=UA60,
        hx_W_K_kg=10.,thermal_converter_We_kg=10.,pump_fraction_hot=0.)
    add('working_200K_60pct_IR11_bypass1','working_reference',**common)
    for t in [197.,200.,220.]:
      for bypass in [1.,10.]:
        name=f'working_{int(t)}K_60pct_IR11_bypass{int(bypass)}'
        if name not in cases:add(name,'working_geometry_family',**(common|{'radiator_K':t,'bypass_W_m2':bypass}))
    for q in [0.,1.,10.,30.,100.]:
        add(f'IR_{int(q)}_200K_60pct','environment_sensitivity',**(common|{'environment_absorbed_W_m2':q}))
    for cost in [3.1,10.,100.]:
        add(f'cost_{cost:g}MJ_200K_60pct','manufacturing_sensitivity',**(common|{'routine_MJ_kg':cost}))
    add('40pct_450K_IR11_bypass1','lower_efficiency_comparison',**(common|{'heat_efficiency':.4,'radiator_K':450.}))
    add('40pct_200K_IR11_bypass1','lower_efficiency_comparison',**(common|{'heat_efficiency':.4}))
    add('65pct_200K_phi85','higher_efficiency_requirement',**(common|{'heat_efficiency':.65,'carnot_quality':.85,'hx_internal_UA_per_hot_W_K':UA65}))
    add('70pct_200K_phi90','higher_efficiency_requirement',**(common|{'heat_efficiency':.7,'carnot_quality':.9,'hx_internal_UA_per_hot_W_K':UA65}))
    add('60pct_250K_phi80_boundary','numerical_boundary_control',**(common|{'radiator_K':250.}))
    add('80pct_150K_phi95','high_efficiency_radiator_control',**(common|{'heat_efficiency':.8,'carnot_quality':.95,'radiator_K':150.,'hx_internal_UA_per_hot_W_K':UA65}))
    for rating in [1.,100.]:
        add(f'HX_rating_{rating:g}','exchanger_rating_sensitivity',**(common|{'hx_W_K_kg':rating}))
    for rating in [1.,100.]:
        add(f'converter_rating_{rating:g}','converter_rating_sensitivity',**(common|{'thermal_converter_We_kg':rating}))
    for pump in [.005,.02]:
        add(f'extra_pumps_{pump:g}','incremental_external_pump_sensitivity',**(common|{'pump_fraction_hot':pump}))
    add('full_outer_area_available','radiator_geometry_sensitivity',**(common|{'radiator_available_fraction':1.}))
    add('escape_view90','radiator_geometry_sensitivity',**(common|{'escape_view_factor':.9}))
    add('IR_charged_only_on_radiators','environment_boundary_sensitivity',**(common|{'environment_basis':'radiator'}))
    delayed=5.678502856908327
    source=cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He']
    add('all_delayed_heat_to_cold_200K','heat_routing_sensitivity',**(common|{'hot_eligible_fraction':1-delayed/source}))
    add('all_delayed_heat_to_cold_220K','heat_routing_sensitivity',**(common|{'hot_eligible_fraction':1-delayed/source,'radiator_K':220.}))
    # Reproduce two original evaluated points without editing original code.
    oldpath=ROOT/'open-loop-redesign/station.py'
    spec=importlib.util.spec_from_file_location('frozen_station',oldpath)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    checks=[]
    for name,t,eta in [('historical_450K_40pct_sun',450.,.4),('historical_197K_60pct_sun',197.,.6)]:
        o=old.evaluate(base['helium_per_s'],base['radius_m'],cycle,inp['receiver'],
            radiator_T=t,hot_T=1000,absorber_sigma=5000,driver_peak_W_kg=100,service_sigma=1000,
            heat_efficiency=eta,replacement_J_kg=1e8)
        for key in ['mass_kg','net_W','net_W_kg','radiator_area_fraction']:
            err=cases[name][key]/o[key]-1
            assert abs(err)<1e-12
            checks.append({'check':name+':'+key,'relative_error':err,'passed':True})
    assert not cases['60pct_250K_phi80_boundary']['feasible']
    checks.append({'check':'60pct at 250K, phi=.8 rejects zero finite-HX margin','passed':True})
    assert abs(cases['historical_450K_40pct_sun']['mass_kg']/base['mass_kg']-1)<1e-12
    for name,c in cases.items():
        if 'closure_relative_residuals' in c:
            assert max(abs(x) for x in c['closure_relative_residuals'].values())<1e-10
    # Direct converter never appears in the ranked/integrated thermal family.
    K_MeV=19.9104911015317
    captured_K=.99*K_MeV*MEV*base['helium_per_s']
    ref=cases['working_200K_60pct_IR11_bypass1']
    bounds=[direct_marginal_bound(ref,captured_K,d) for d in [.5,.7,.85]]
    data={'scope':'Bounded conditional thermal station family at fixed pulse rate and radius. 60% package conversion is the working premise; higher efficiencies and apparatus ratings are explicit scenarios. No global optimum or demonstrated plant is claimed.',
        'environment_reference':'Zero direct sunlight. Explicit11.0254 W/m² absorbed IR per full sphere is a gray compact Solar-IR benchmark, not a measured local sky. Zero-IR rows are optimistic controls.',
        'hardware_scope':'Original generic service, source, storage and insulation allocations retained. Thermal converter replaces the old gross-electric converter allocation. HX and optional extra pump hardware are incremental. Factory mass scales with repaired throughput independently of its electricity cost.',
        'additional_internal_UA_source':{'60pct_W_K_per_hot_W':UA60,'65plus_W_K_per_hot_W':UA65,
             'meaning':'Recuperator plus generator-cooling conductance from independent conditional six-stage gas-cycle cases is added to the optimized external-exchanger minimum. The 220 K and 70% cases are temperature/conductance requirements, not complete matching solved gas cycles.'},
        'cases':cases,'direct_marginal_energy_bounds':bounds,
        'direct_scope':'No direct-converter point is ranked. Peak handled power, buffering/export and added equipment must be resolved before turning these bounds into a plant.',
        'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [bp,ip,hp,oldpath]},
        'code_hashes':{p.name:sha(p) for p in [Path(__file__),HERE/'station.py']},
        'checks':checks,'all_closure_checks_passed':True}
    (HERE/'cases.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    for name,c in cases.items():
        print(name,{'feasible':c['feasible'],'mW_kg':round(c.get('net_W_kg',0)*1e3,6),
              'rad_fraction':round(c.get('radiator_area_fraction',0),6),'reasons':c['reasons']})

if __name__=='__main__':main()
