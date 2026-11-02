"""Conditional 60%/50 K ion commissioning ledger, including an 80 AU sky.

Uses the saved source family; no new plasma, material or nuclear-data fit.
The mass/power reserve and particle-processing costs are selected requirements.
Run with docs/protium-fusion and /tmp/protium-python on PYTHONPATH.
"""
from pathlib import Path
import gzip, hashlib, json, math, sys
from copy import deepcopy
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'phase-two/plant-audit'))
from network_plant import plant, YEAR, SB
from nonneutral_screen import MEV, C
from catalyst_catalog import atomicmass
from heat_transport import thermal_power, exchanger, optimize_loop, additional_cost_case

CONFIG = dict(heat=.6, radiator_K=50.)
RESERVE_MASS_KG = 1e14
RESERVE_POWER_W = 1e11
SOLAR_1AU_W_M2 = 1361.
DISTANCE_AU = 80.
checks = []

def near(name, x, y, rtol=2e-12, atol=0.):
    assert abs(x-y) <= atol + rtol*max(abs(x),abs(y),1.), (name,x,y)
    checks.append(name)

def close(r, zero_extra, absorption, fixed_other):
    a=r['assumptions']; old=r['mass_breakdown_kg']
    masses={k:v for k,v in old.items() if k not in
            ('radiators','replacement_manufacturing_equipment','replacement_inventory')}
    masses['additional_integration_reserve']=RESERVE_MASS_KG
    extra_aux=max(0.,zero_extra['other_auxiliary_W']-r['other_auxiliary_W']) if fixed_other else 0.
    bus=r['plant_bus_net_W']-RESERVE_POWER_W-extra_aux
    emit=a['radiator_emissivity']*SB*(a['radiator_K']**4-a['ambient_K']**4)
    incident=SOLAR_1AU_W_M2/DISTANCE_AU**2
    absorbed_flux=absorption*incident/4  # outward sphere: projected area = emitting area / 4
    usable_flux=emit-absorbed_flux
    assert usable_flux>0
    k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR
    spare=a['annual_replacement_fraction']*a['replacement_inventory_years']
    radkg=a['radiator_areal_mass']/usable_flux
    denominator=1-k*(radkg+a['manufacturing_kg_W'])-spare
    mass=(sum(masses.values())+(r['nuclear_after_capture_W']-bus)*radkg)/denominator
    maintenance=k*mass; net=bus-maintenance
    internal=r['nuclear_after_capture_W']-net
    area=internal/usable_flux
    masses.update(radiators=area*a['radiator_areal_mass'],
                  replacement_manufacturing_equipment=a['manufacturing_kg_W']*maintenance,
                  replacement_inventory=spare*mass)
    near('mass closure',sum(masses.values()),mass)
    near('global heat closure',area*emit,internal+area*absorbed_flux)
    near('electrical closure',net+maintenance,bus)
    if not absorption and not fixed_other:
        oldclose=additional_cost_case(r,{'pipe_and_fluid_mass_kg':RESERVE_MASS_KG,
                                        'pump_electric_W':RESERVE_POWER_W})
        near('existing additional_cost_case mass',mass,oldclose['total_mass_kg'])
        near('existing additional_cost_case power',net,oldclose['net_W'])
    return dict(extra_processing_MeV_per_primary=a['extra_processing_MeV_per_primary'],
        auxiliary_policy='hold zero-extra auxiliary draw fixed' if fixed_other else 'original 10% scaling',
        bolometric_absorptivity=absorption,net_W=net,total_mass_kg=mass,net_mW_kg=1000*net/mass,
        bus_before_replacement_W=bus,replacement_manufacturing_W=maintenance,
        restored_generic_auxiliary_W=extra_aux,extra_processing_W=r['extra_processing_W'],
        generic_other_auxiliary_W=r['other_auxiliary_W']+extra_aux,
        nuclear_after_capture_W=r['nuclear_after_capture_W'],internal_waste_W=internal,
        radiator_area_m2=area,radiator_sphere_radius_m=math.sqrt(area/(4*math.pi)),
        radiator_emission_flux_W_m2=emit,solar_directional_flux_W_m2=incident,
        external_absorbed_flux_mean_W_m2=absorbed_flux,external_absorbed_W=area*absorbed_flux,
        emitted_W=area*emit,usable_radiator_flux_mean_W_m2=usable_flux,
        allowed_internal_flux_subsolar_W_m2=emit-absorption*incident,
        allowed_internal_flux_antisolar_W_m2=emit,
        construction_energy_J=mass*a['embodied_J_kg'],
        construction_energy_payback_years=mass*a['embodied_J_kg']/net/YEAR,
        annual_replacement_kg=a['annual_replacement_fraction']*mass,
        mass_breakdown_kg=masses)


def main():
    energy=json.loads((ROOT/'network-ion-budget.json').read_text())
    family=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
    baseline=plant(energy,family,**CONFIG)
    source={x:plant(energy,family,**CONFIG,extra_processing_MeV_per_primary=x) for x in (1.,1.2,1.5,1.75,2.)}
    cases=[close(r,baseline,alpha,fix) for r in source.values() for alpha in (0.,1.) for fix in (False,True)]
    selected=next(x for x in cases if x['extra_processing_MeV_per_primary']==1.75 and
                  x['bolometric_absorptivity']==1 and x['auxiliary_policy'].startswith('hold'))
    r=source[1.75];rate=r['primary_pn_rate_s'];n=energy['network']
    for other in source.values():
        near('same source rate',other['primary_pn_rate_s'],rate)
        near('same nuclear power',other['nuclear_after_capture_W'],r['nuclear_after_capture_W'])
    # A two-leg ideal He pipe screen, assessed at the complete selected radiator radius.
    # The source engine rejects only part of the plant's heat. A second comparison
    # sends ALL internally generated waste through the cold loop and therefore
    # overcounts its own friction already included in the 100 GW reserve.
    h=thermal_power(r)
    x=exchanger(h['hot_input_W'],eta=.6,relative_Carnot=.85,Th=300.,Tr=50.)
    loop_r=deepcopy(r);loop_r['radiator_area_m2']=selected['radiator_area_m2']
    loops=[]
    for f in (.005,.01):
        for all_waste in (False,True):
            loop_h=deepcopy(h)
            if all_waste:loop_h['engine_cold_rejection_W']=selected['internal_waste_W']
            z=optimize_loop(loop_r,loop_h,x,span_fraction=.5,darcy_f=f)
            z['load_scope']='all internal waste; self-friction deliberately counted twice' if all_waste else 'engine rejection only'
            loops.append(z)
    inventories=[]
    for row in n['decay_inventories']:
        tau=row['mean_life_s']; m=rate*row['kg_per_primary_pn_s']
        if not m:continue
        za=divmod(int(row['nuclide']),1000)
        near('inventory activity '+row['nuclide'],m/(atomicmass(*za)*MEV/C**2)/tau,
             rate*row['decays_per_primary_pn'])
        inventories.append(dict(row,inventory_kg=m,equilibrium_activity_s=rate*row['decays_per_primary_pn'],
            half_life_days=tau*math.log(2)/86400,mean_life_years=tau/YEAR,
            time_to_90_percent_years=tau*math.log(10)/YEAR,
            time_to_99_percent_years=tau*math.log(100)/YEAR,
            buildup=[dict(years=t,fraction=1-math.exp(-t*YEAR/tau),
                     inventory_kg=m*(1-math.exp(-t*YEAR/tau))) for t in (0,1,2,3,5,6)]))
    near('decay inventory sum',sum(z['inventory_kg'] for z in inventories),r['mass_breakdown_kg']['cold_decay_inventory'])
    v=next(z for z in inventories if z['nuclide']=='23049')
    ti49_mass=atomicmass(22,49)*MEV/C**2
    startup=dict(cold_decay_inventory_total_kg=sum(z['inventory_kg'] for z in inventories),
        V49_inventory_kg=v['inventory_kg'],V49_flow_mass_kg_s=v['inventory_kg']/v['mean_life_s'],
        Ti49_maximum_bridge_kg=v['equilibrium_activity_s']*v['mean_life_s']*ti49_mass,
        Ti49_initial_bridge_consumption_kg_s=v['equilibrium_activity_s']*ti49_mass,
        steady_Ti49_boundary_makeup_kg_year=r['capture_bank']['Ti49_makeup_kg_year'],
        V49_decay_buildup=v['buildup'],all_decay_inventories=inventories,
        capture_bank=r['capture_bank'],
        neutral_store_note='Neutral/cold means outside the MeV ion plasma. Store heat at the warm process temperature, not on the 4.2 K coil circuit.',
        startup_scope='Independent inventory fill law under fixed source flows; not a coupled transient isotope-network solution. Stable bridge replaces the future decay inventory, so its mass is not added twice.',
        weak_energy_note='The underlying ledger assigns all weak ground-state Q and a known-isomer margin to escaping neutrinos, so no delayed weak-decay heat is credited. Actual electromagnetic decay energy requires warm-store heat handling and a later partition refinement.')
    provenance={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
        ROOT/'network_plant.py',ROOT/'network-ion-budget.json',ROOT/'network-channel-family.json.gz',
        ROOT/'catalyst_network.py',ROOT/'capture_bank_screen.py',ROOT/'nonneutral_screen.py',
        ROOT/'bundle_solenoid.py',ROOT/'phase-two/plant-audit/heat_transport.py',Path(__file__)]}
    out=dict(scope=__doc__,selected_engineering_requirements=dict(CONFIG,
        additional_integration_mass_kg=RESERVE_MASS_KG,additional_integration_power_W=RESERVE_POWER_W,
        extra_processing_cases_MeV_per_primary=[1.,1.2,1.5,1.75,2.],relative_internal_Carnot=.85,
        solar_1AU_W_m2=SOLAR_1AU_W_M2,distance_AU=DISTANCE_AU,
        environmental_geometry='outward spherical emitting area; full solar luminosity envelope, no credit for swarm interception; uniform 50 K requires sector-specific coolant loading'),
        retained_plant_assumptions=r['assumptions'],source_rate_s=rate,source_field_T=r['source_field_T'],
        source_loss_MeV_per_primary=r['source_loss_MeV_per_primary'],
        stored_field_energy_J=r['stored_field_energy_J'],refrigeration_W=r['refrigeration_W'],
        cases=cases,selected_conservative_case=selected,engine=h,exchanger=x,pipe_screens=loops,
        startup=startup,validation=dict(check_count=len(checks),all_passed=True),input_sha256=provenance,
        primary_sources=[dict(url='https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html',use='1361 W/m2 solar irradiance at Earth'),
          dict(url='https://www.isotopes.gov/sites/default/files/2021-02/V-49.pdf',use='DOE V49 product information: 330 day half-life, Ti49 daughter, production route'),
          dict(url='https://pubmed.ncbi.nlm.nih.gov/10724407/',use='Primary half-life measurement 328.6 +/- 3.3 d; saved network uses NUBASE-derived 330 d')])
    path=HERE/'selected-ion-point.json'
    path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(selected=selected,engine=h,exchanger=x,
        pipe_screens=[{k:z[k] for k in ('load_scope','darcy_f','pipe_and_fluid_mass_kg','pump_electric_W','diameter_per_pipe_m','velocity_m_s')} for z in loops],
        startup={k:startup[k] for k in ('cold_decay_inventory_total_kg','V49_inventory_kg','Ti49_maximum_bridge_kg','Ti49_initial_bridge_consumption_kg_s','steady_Ti49_boundary_makeup_kg_year')},
        capture_stock_kg=r['capture_bank']['total_listed_inventory_lower_bound_kg'],checks=len(checks)),indent=2))

if __name__=='__main__':main()
