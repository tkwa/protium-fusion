#!/usr/bin/env python3
"""Conditional architecture accounts; pulse-size extrapolations are requirements.

All source energy, including recovered driver energy, enters once. Mass ratings
are explicit acceptance allocations. No Monte Carlo whole-design interval.
"""
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path

P=Path(__file__).resolve().parent; ROOT=P.parents[2]
BASE=ROOT/'docs/protium-fusion/mixture-design/selected-plant/point.json'
DISPATCH=ROOT/'docs/protium-fusion/component-design/final-family-dispatch/sigma350/input.json'
B=json.loads(BASE.read_text()); D=json.loads(DISPATCH.read_text()); REF=B['plant']; C=B['cycle']
AREA=REF['area_m2']; YEAR=365.25*86400; SB=5.670374419e-8
spec=importlib.util.spec_from_file_location('engine',P.parent/'cold-sink-redesign/heat_engine.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
with gzip.open(P/'trajectories.json.gz','rt') as f:T=json.load(f)

def cycle(quality='central'):
    strong=quality=='strong'
    return engine.best_cycle(ts=1600.,tr=400.,th=1500.,tc=450.,nc=6,nt=6,
        ec=.97 if strong else .93,et=.97 if strong else .93,
        er=.98 if strong else .95,eg=.99 if strong else .98,dp=.02 if strong else .03,q=1.)

def cold_requirements(area_factor=16.,guarded=True,helium_Carnot_fraction=.3):
    """Areal auxiliary requirements from the newer warm-buffer component account.
    A dedicated 200 K radiation intercept can be cooled against the 400 K sink.
    Guide and 5 K refrigerators reject directly to 400 K; their work is not
    routed through a second inefficient refrigerator.
    """
    f=16.;cop65old=.15*65/(200-65);cop65=.15*65/(400-65)
    ir65old=(.31945831170032074+.36922118719844016)*cop65old
    photon65old=.3956152524360582*cop65old
    shield=200. if guarded else 400.
    guide=(ir65old*(shield/200.)**4+f*photon65old)/cop65
    cop5=helium_Carnot_fraction*5/(400-5)
    helium=.004433445869545912*area_factor*(shield/200.18862258959285)**4/cop5
    loads=dict(guide65K=guide,helium5K=helium,
        helium_bank_cycle=.031050635753950923*area_factor*(2.5 if guarded else 32.),
        warm_motors=1.13*f,defect_repair=.18940849643963228*f,
        face_renewal=.02326076076076076,monthly_rework=.22815423226100842,
        gas_selective_work=.004938396589206173*f,downstream_reserve=.005*f,
        getter_compression=5.7994072258631276e-5*f)
    if guarded:
        cop200=.3*200/(400-200)
        guide_area=ir65old/(.01*SB*(200**4-65**4))
        loads['guide200K_intercept']=guide_area*.01*SB*(400**4-200**4)/cop200
        loads['helium200K_intercept']=4.8667e-5*area_factor*.9*SB*(400**4-200**4)/cop200
        loads['hot_guide_shield200K_intercept']=.26051007830783396*(1600**4-200**4)/(1000**4-200**4)/cop200
    return dict(guarded_200K_intercepts=guarded,helium_active_area_multiplier=area_factor,
        helium_Carnot_fraction=helium_Carnot_fraction,components_W_m2=loads,
        named_electric_W_m2=sum(loads.values()),
        qualification='Linear photon/particle throughput, unchanged guide interception and chosen gas collection/repair laws. Area multiplication is a capacity scenario. New temperature-dependent leakage, bank cycling, capture and refrigerator apparatus need design. Coolant header pumping is inside the engine package pump allocation.')

def schedule(rate=16.,n15_mass=100.,overlap=False,centres=1,days=370.,name=None,new_cycle=None):
    horizon=days*86400.; phases={'c12':.37,'c13':.05,'n14':.12,'n15':.5};raw=[];windows={}
    for n in phases:
        scale=n15_mass if n=='n15' else 1.
        correction=C['stage_attempts_per_He'][n]/new_cycle['stage_attempts_per_He'][n] if new_cycle else 1.
        period=REF['pulse_interval_s'][n]*scale/rate*correction
        w={k:D['timing'][n]['mass999'][k]*1.1 for k in ['incoming_s','outgoing_s']}
        windows[n]=w
        for j in range(int(horizon/period)+2):
            t=(j+phases[n])*period
            if t<=horizon:raw.append((t,n))
    raw.sort();free=[0.]*centres;events=[]
    for nominal,n in raw:
        w=windows[n];lane=min(range(centres),key=lambda i:free[i])
        # 'overlap' is a purely kinematic lower bound: outgoing/incoming interactions unsolved.
        start=max(nominal-w['incoming_s'],free[lane]);focus=start+w['incoming_s'];end=focus+w['outgoing_s']
        free[lane]=focus if overlap else end
        events.append(dict(stage=n,nominal_focus_s=nominal,launch_s=start,focus_s=focus,
                           capture_end_s=end,delay_s=focus-nominal,centre=lane))
    occupancy=sum((w['incoming_s']+(0 if overlap else w['outgoing_s']))/
        (REF['pulse_interval_s'][n]*(n15_mass if n=='n15' else 1)/rate*
         (C['stage_attempts_per_He'][n]/new_cycle['stage_attempts_per_He'][n] if new_cycle else 1.)) for n,w in windows.items())
    if name:
        with (P/(name+'-events.csv')).open('w') as f:
            wr=csv.DictWriter(f,fieldnames=events[0].keys());wr.writeheader();wr.writerows(events)
    return dict(rate_multiplier=rate,n15_mass_multiplier=n15_mass,centres=centres,
        overlapping_legs=overlap,event_count=len(events),load_per_centre=occupancy/centres,
        max_delay_s=max(e['delay_s'] for e in events),finish_overrun_s=max(0,max(e['capture_end_s'] for e in events)-horizon),
        stable_average_load=occupancy/centres<1,
        minimum_focus_gap_s=min(b-a for a,b in zip(sorted(e['focus_s'] for e in events),sorted(e['focus_s'] for e in events)[1:])),
        counts={n:sum(e['stage']==n for e in events) for n in windows})

def account(name,f=16.,n15_mass=100.,eta=None,direct=0.,peak_rating=1000.,storage_J_kg=1e5,
            service=1000.,bypass=6.55,thermal_rating=10.,hx_rating=10.,hot=1600.,cold=400.,
            direct_storage='sufficient',availability=1.,source_peak_factor=1.,quality='central',new_cycle=None,aux_floor=0.):
    ec=cycle(quality);eta=eta if eta is not None else ec['net_converter_efficiency_after_declared_losses']
    # Per receiver area. Prepared electricity is debited once, even if paid by direct output.
    cy=new_cycle or C
    H=REF['helium_per_s']/AREA*f*cy['H_nonnu_MeV_per_He']*1.602176634e-13
    U=REF['helium_per_s']/AREA*f*cy['U_prepared_MeV_per_He']*1.602176634e-13
    NU=REF['helium_per_s']/AREA*f*cy['neutrino_MeV_per_He']*1.602176634e-13;driver=U/.6
    supplied=H+U;eligible=REF['parameters']['hot_eligible_fraction']*supplied
    kinetic=cy['bulk_kinetic_MeV_per_He']/(cy['H_nonnu_MeV_per_He']+cy['U_prepared_MeV_per_He'])*supplied
    direct_gross=direct*kinetic;storage_loss=.1*direct_gross
    thermal=eligible-direct_gross-bypass
    assert thermal>0
    gross_h=eta*thermal; gross_d=direct_gross-storage_loss
    margin=gross_h+gross_d-driver
    # Covers pump/rework/processing requirements once; new direct modules get extra 1% input.
    ordinary=max(.1*margin,aux_floor);direct_aux=.01*kinetic if direct else 0.
    before=margin-ordinary-direct_aux
    stage_K={n:T[n]['ballistic']['kinetic_J_exact_nodal']*(n15_mass if n=='n15' else 1) for n in T}
    peak_K=max(T[n]['ballistic']['peak_kinetic_flux_W_m2']/(1.1**2)*(n15_mass if n=='n15' else 1) for n in T)
    # Exact final-family precursor B1 source rating supersedes its hot-start rating.
    source_peak=max(REF['peak_source_electric_W'],2.59427e21/1.1*n15_mass)*source_peak_factor
    max_prepared=max(T[n]['prepared_J']*(n15_mass if n=='n15' else 1) for n in T)
    store_driver=max(REF['stored_electric_J'],2*max_prepared/.6)
    out_store=(sum(stage_K.values()) if direct_storage=='sufficient' else max(stage_K.values()))*direct/AREA
    if direct_storage=='sufficient':out_store+=direct_gross*8951.35
    refmass={k:v/AREA for k,v in REF['mass_breakdown_kg'].items()}
    masses={k:refmass[k] for k in ['receiver','insulation']}
    masses['processing_structure_thermal_service']=service
    masses['peak_driver_hardware']=source_peak/100./AREA
    masses['two_pulse_driver_storage']=store_driver/1e5/AREA
    masses['average_driver_conditioning']=driver/10.
    masses['working_and_ready_inventory']=(REF['working_inventory_kg'] if 'working_inventory_kg' in REF else REF['mass_breakdown_kg']['working_inventory'])*f/AREA
    masses['fresh_and_flight_inventory']=(refmass['inflight_fuel_inventory']+refmass['fresh_feed_buffer'])*f
    masses['thermal_converter']=gross_h/thermal_rating
    masses['heat_exchangers']=ec['UA_total_W_K_m2']*thermal/hx_rating
    masses['direct_converter_input_peak']=peak_K/peak_rating if direct else 0.
    masses['direct_output_storage']=out_store/storage_J_kg
    # Charge/discharge hardware independently rated. Conservative sum can overlap converter conditioner.
    masses['direct_store_charger']=direct*peak_K/1e4 if direct else 0.
    masses['direct_store_discharge']=direct_gross/1000. if direct else 0.
    m0=REF['mass_kg']/AREA;mass=m0;support0=REF['support_W']/AREA
    for _ in range(100):
        support=support0*(mass/m0)**2;repair=mass*1e-4/YEAR*2e6
        power=before-repair-support
        waste=H-power # closed plant: nuclear non-neutrino supply = export + waste
        radarea=(waste+11.0254)/(.9*SB*cold**4)
        masses['radiator_panels']=20*radarea
        masses['replacement_factory']=mass*1e-4/YEAR*1e8
        masses['ten_year_spares']=.001*mass;masses['installed_reserve']=.2*mass
        masses['support_streams']=refmass['support_streams']*(mass/m0)**2
        updated=sum(masses.values())
        if abs(updated-mass)<1e-9:mass=updated;break
        mass=updated
    cold_components=dict(engine=(1-eta)*thermal,missed_hot=supplied-eligible,bypass=bypass,
        driver_loss=driver-U,ordinary=ordinary,direct_auxiliary=direct_aux,
        output_storage_loss=storage_loss,repair=repair,support=support)
    assert abs(sum(cold_components.values())-waste)<1e-7
    assert abs(sum(masses.values())-mass)<1e-6
    assert eta<1-cold/hot
    rates=json.loads((P/'rates.json').read_text())
    amass={n:m/931.49410242*1.66053906660e-27 for n,m in zip(rates['species'],rates['mass_MeV'])}
    count_rate=REF['helium_per_s']*f
    flows=dict(gross_fuel_kg_s=cy['gross_fuel_kg_per_He']*count_rate,
        recovered_CNO_kg_s=sum(cy['returned_atoms_per_He'].get(n,0.)*amass[n] for n in cy['states'])*count_rate,
        recovered_H_kg_s=cy['returned_atoms_per_He']['p']*amass['p']*count_rate,
        fresh_H_kg_s=cy['fresh_hydrogen_kg_per_He']*count_rate,
        fresh_C_quartet_boundary_kg_s=cy['fresh_carbon_kg_per_He']*count_rate,
        retained_O16_kg_s=cy['ash_atoms_per_He'].get('o16',0.)*amass['o16']*count_rate)
    return dict(name=name,rate_multiplier=f,n15_mass_multiplier=n15_mass,eta_heat=eta,
        direct_fraction_of_all_K=direct,availability=availability,net_W_m2=power,
        installed_kg_m2=mass,net_W=power*AREA,mass_kg=mass*AREA,net_W_kg=power/mass,
        calendar_W_kg=availability*power/mass,gain_over_reference=(power/mass)/REF['net_W_kg'],
        radiator_fraction=radarea,radiator_fits=radarea<=.9,
        mass_breakdown_kg_m2=masses,cold_heat_W_m2=cold_components,
        energy_W_m2=dict(nuclear_nonnu=H,neutrino=NU,prepared=U,driver_electric=driver,
            eligible_before_direct=eligible,kinetic=kinetic,direct_gross=direct_gross,
            direct_after_storage=gross_d,engine_hot=thermal,thermal_electric=gross_h,
            net=power,waste=waste,environment=11.0254),
        equipment=dict(source_peak_electric_W=source_peak,source_average_electric_W=driver*AREA,
            driver_storage_J=store_driver,receiver_peak_K_W_m2=peak_K,
            direct_input_peak_W=peak_K*AREA if direct else 0.,
            direct_output_store_J=out_store*AREA,largest_K_pulse_J=max(stage_K.values()),
            engine_heat_W=thermal*AREA,engine_electric_W=gross_h*AREA,
            exchanger_UA_W_K=ec['UA_total_W_K_m2']*thermal*AREA),flows_kg_s=flows)

def main():
    cases=[account('thermal16'),account('thermal16_services700',service=700.),
        account('hybrid16',direct=.7),account('hybrid16_10kWkg',direct=.7,peak_rating=1e4),
        account('hybrid16_store1MJkg',direct=.7,storage_J_kg=1e6),
        account('thermal20',f=20.),account('thermal16_strong_engine',quality='strong'),
        account('thermal16_availability90',availability=.9),
        account('thermal16_source_x2',source_peak_factor=2.)]
    if (P/'probe-cycles.json').exists():
        probe=json.loads((P/'probe-cycles.json').read_text())
        new=next(x['cycle'] for x in probe['cases'] if x['record']=='pulse-n15-m100-v1-z128-t1.json')
        cases +=[account('spherical_large_N15_16',new_cycle=new),
                 account('spherical_large_N15_16_services700',new_cycle=new,service=700.),
                 account('hybrid16_both_strong_ratings',direct=.7,peak_rating=1e4,storage_J_kg=1e6)]
        finer=next((x['cycle'] for x in probe['cases'] if x['record']=='pulse-n15-m100-v1-z256-t1.json'),None)
        if finer:cases +=[account('spherical_large_N15_256_16',new_cycle=finer)]
        cases +=[account('spherical16_unguarded_cold_area16',new_cycle=new,aux_floor=cold_requirements(guarded=False)['named_electric_W_m2']),
                 account('spherical16_guarded_He_fridge15pctCarnot',new_cycle=new,aux_floor=cold_requirements(helium_Carnot_fraction=.15)['named_electric_W_m2'])]
    sched=[schedule(f,g) for f in [1.,4.,8.,16.,20.] for g in [1.,10.,30.,100.]]
    sched +=[schedule(16.,1.,overlap=True),schedule(16.,1.,centres=11),schedule(16.,100.,name='large-n15-16x')]
    if (P/'probe-cycles.json').exists():sched +=[schedule(16.,100.,new_cycle=new,name='spherical-large-n15-16x')]
    out=dict(scope=__doc__,central_cycle=cycle(),strong_cycle=cycle('strong'),cases=cases,schedules=sched,
        cold_equipment_requirements=[cold_requirements(area_factor=a,guarded=g) for a in [1.,16.] for g in [False,True]],
        target_W_kg=16*.0208,reference=B['plant']['net_W_kg'],
        input_sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [BASE,DISPATCH,P/'trajectories.json.gz']},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (P/'plant-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for c in cases:print(c['name'],round(c['net_W_kg'],4),round(c['gain_over_reference'],2),round(c['radiator_fraction'],3))
    print('cycle',cycle()['net_converter_efficiency_after_declared_losses'],cycle('strong')['net_converter_efficiency_after_declared_losses'])

if __name__=='__main__':main()
