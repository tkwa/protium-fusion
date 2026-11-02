"""Conditional whole-system budget for a very large thermal pulse collector.

This is a requirements allocation, not an engineered or validated plant. In
particular, the fixed driver and processing masses are explicit allowances,
and the 30% driver efficiency is a hypothesis. The calculation asks whether
those assumptions can close below Jupiter mass at 1 mW/kg, including support
losses, recovery losses, pulse storage and gross-rated balance of plant.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from cycle_ledger import ledger,nuc,MU,MEV

ROOT=Path(__file__).parent
G=6.67430e-11;SIGMA_SB=5.670374419e-8;MJ=1.89813e27

def budget(driver_mass=1e25,eta_driver=.3,collection_loss=1e-4,
           stream_loss_per_orbit=.001,prepared_energy_multiplier=1.):
    radius=2e10;area=4*np.pi*radius**2;specific_power=.001
    eta_heat=.4;energy_capture=.99;auxiliary_fraction=.1
    stream_speed=1e4;storage_specific_energy=1e5
    d=ledger(collection_loss=collection_loss,factory_file='compression-He-64-r30-burning.json')
    u=d['fuel_preparation_MeV_per_helium']*prepared_energy_multiplier
    # Escaping prepared energy was supplied by the recirculating electrical
    # system and must also be deducted from the nuclear energy available to it.
    q=d['available_heat_MeV_per_helium']*energy_capture-(1-energy_capture)*u
    generated=eta_heat*(q+u);driver=u/eta_driver
    net=generated-driver
    if net<=0: return {'closes':False,'reason':'Negative cycle energy before other auxiliaries',
      'driver_efficiency':eta_driver,'prepared_energy_multiplier':prepared_energy_multiplier,'net_MeV':net}
    input_records=[json.loads((ROOT/s['source']).read_text()) for s in d['stage_data']]
    input_records.append(json.loads((ROOT/d['factory_data']['source']).read_text()))
    max_pulse_electric=max(r['initial_energy_J'] for r in input_records)/eta_driver*prepared_energy_multiplier
    fixed={'absorber':area*1e4,'insulation':area*100,
           'processing_and_inventory_allowance':area*1000,
           'support_guides_controls_and_redundancy_allowance':area*100,
           'driver_facility_allowance':driver_mass,
           'pulse_electrical_energy_storage':max_pulse_electric/storage_specific_energy}
    gen_kg_per_W=.1;power_conditioning_kg_per_W=.1
    radiator_T=300.;radiator_emissivity=.9;radiator_areal_mass=10.
    rad_flux=radiator_emissivity*SIGMA_SB*radiator_T**4
    def state(mass):
        stream_energy=G*mass**2/(4*radius)
        stream_mass=2*stream_energy/stream_speed**2
        orbit=2*np.pi*radius/stream_speed
        support_loss=stream_loss_per_orbit*stream_energy/orbit
        delivered=specific_power*mass
        throughput=(delivered+support_loss)/((1-auxiliary_fraction)*net*MEV)
        gross=generated*MEV*throughput;drive=driver*MEV*throughput
        heat=q*MEV*throughput-delivered
        variable={'circulating_support_streams':stream_mass,
          'thermal_power_conversion':gen_kg_per_W*gross,
          'driver_average_power_conditioning':power_conditioning_kg_per_W*drive,
          'radiators':radiator_areal_mass*heat/rad_flux}
        return sum(fixed.values())+sum(variable.values()),variable,{
          'net_useful_power_W':delivered,'helium_atoms_per_s':throughput,
          'gross_electric_generation_W':gross,'driver_electric_power_W':drive,
          'waste_heat_W':heat,'support_loss_W':support_loss,
          'other_auxiliary_reserve_W':auxiliary_fraction*net*MEV*throughput,
          'support_stream_kinetic_energy_J':stream_energy,'support_orbit_s':orbit,
          'hydrogen_consumption_kg_s':throughput*d['hydrogen_atoms_per_net_helium']*nuc('p').A*MU}
    lo=sum(fixed.values());grid=np.geomspace(lo,MJ,300)
    brackets=[(a,b) for a,b in zip(grid,grid[1:]) if (state(a)[0]-a)*(state(b)[0]-b)<0]
    if not brackets: return {'closes':False,'reason':'No whole-mass fixed point below Jupiter',
      'driver_mass_allowance_kg':driver_mass,'driver_efficiency':eta_driver}
    mass=brentq(lambda m:state(m)[0]-m,*brackets[0],xtol=1e8)
    total,variable,operations=state(mass)
    periods={}
    for st,r in zip(d['stage_data'],input_records):
        a=sum(nuc(k).A*v for k,v in st['fuel_input'].items())
        atoms=r['mass_kg']/(a*MU)
        periods[st['cat']]=atoms/(operations['helium_atoms_per_s']*d['stage_attempts_per_net_helium'][st['cat']])
    r=input_records[-1]
    carbon_atoms=r['mass_kg']*r['output_mass_fractions']['c12']/(12*MU)
    periods['carbon_factory']=carbon_atoms/(operations['helium_atoms_per_s']*d['new_carbon_per_net_helium'])
    operations['radiator_area_m2']=operations['waste_heat_W']/rad_flux
    operations['radiator_area_over_collector_area']=operations['radiator_area_m2']/area
    operations['mean_absorber_heat_flux_W_m2']=(q+u)*MEV*operations['helium_atoms_per_s']/area
    operations['mean_interval_between_all_pulses_s']=1/sum(1/t for t in periods.values())
    warnings=[]
    if net<1.: warnings.append('Small energy margin; model/driver uncertainty can reverse the sign')
    if operations['radiator_area_over_collector_area']>.5:
        warnings.append('Radiator area exceeds half the collector sphere; separate view-factor and support geometry is required')
    return {'closes':True,'scope':'Conditional requirements allocation, not a validated plant',
      'parameters':{'collector_radius_m':radius,'absorber_areal_mass_kg_m2':1e4,
       'thermal_conversion_efficiency':eta_heat,'driver_efficiency':eta_driver,
       'energy_capture_fraction':energy_capture,'catalyst_collection_loss_per_stage':collection_loss,
       'other_auxiliary_fraction_of_preaux_net':auxiliary_fraction,
       'prepared_energy_multiplier':prepared_energy_multiplier,
       'support_stream_speed_m_s':stream_speed,'support_fractional_energy_loss_per_orbit':stream_loss_per_orbit,
       'installed_storage_specific_energy_J_kg':storage_specific_energy,
       'gross_generator_kg_W':gen_kg_per_W,'driver_conditioning_kg_W':power_conditioning_kg_per_W,
       'radiator_T_K':radiator_T,'radiator_emissivity':radiator_emissivity,'radiator_kg_m2':radiator_areal_mass},
      'mass_budget_kg':{**fixed,**variable},'total_mass_kg':mass,'total_mass_Jupiter':mass/MJ,
      'mass_fixed_point_relative_residual':total/mass-1,
      'net_specific_power_W_kg':operations['net_useful_power_W']/mass,
      'cycle_heat_MeV_per_helium':q,'prepared_energy_MeV_per_helium':u,
      'preaux_net_MeV_per_helium':net,'max_electric_pulse_energy_J':max_pulse_electric,
      'operations':operations,'pulse_intervals_s':periods,
      'allocation_warnings':warnings,
      'unproven_requirements':['Driver creates the specified low-entropy convergent flow at the stated coupling and within its mass allocation',
       'Pulse remains sufficiently symmetric during compression and expansion',
       'Collector captures 99% of deposited energy and recovers 99.99% of circulating catalyst per processing stage',
       'Distributed absorber, insulation and processing function at their allocated areal masses',
       'Active orbital support realizes the ideal momentum balance with its stated guide mass and losses',
       'Pulse storage achieves the installed specific energy and power requirements',
       'Gross-rated generators, conditioning and radiators achieve the stated mass coefficients']}

if __name__=='__main__':
    cases=[budget(),budget(driver_mass=1e26),budget(eta_driver=.2),
           budget(collection_loss=1e-3),budget(stream_loss_per_orbit=.01),
           budget(prepared_energy_multiplier=1.5),budget(prepared_energy_multiplier=2.)]
    (ROOT/'plant-budget.json').write_text(json.dumps({'cases':cases},indent=2,allow_nan=False)+'\n')
    for c in cases:
        print({k:c[k] for k in ['closes','total_mass_kg','total_mass_Jupiter','preaux_net_MeV_per_helium'] if k in c})
        if c['closes']: print('A interval days',c['pulse_intervals_s']['n14']/86400,
          'support/output',c['operations']['support_loss_W']/c['operations']['net_useful_power_W'])
