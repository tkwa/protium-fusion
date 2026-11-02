"""Additional thermal cavity-source robustness at the fixed final350 schedule.

Read-only reuse of the frozen gas packet. No total release probability is fitted.
"""
from pathlib import Path
import json,math,hashlib,sys
import numpy as np
P=Path(__file__).resolve().parent;G=P.parent
sys.path.insert(0,str(G));import chamber
k=1.380649e-23

def steady_column(radius_m, temperature_K, source_atomic_kg_s, gas_particle_mass_kg, effective_capture):
 """Uniform Maxwellian radial mass column; arrays may contain any species."""
 q=np.asarray(source_atomic_kg_s,dtype=float);m=np.asarray(gas_particle_mass_kg,dtype=float);g=np.asarray(effective_capture,dtype=float)
 if np.any(q<0) or np.any(m<=0) or np.any(g<=0) or np.any(g>1):raise ValueError('nonphysical gas source, mass, or effective capture')
 speed=np.sqrt(8*k*temperature_K/(math.pi*m));area=4*math.pi*radius_m**2
 col=4*radius_m/area*q/(g*speed)
 return {'species_column_kg_m2':col.tolist(),'total_column_kg_m2':float(sum(col)),'species_pressure_Pa':(col*k*temperature_K/(m*radius_m)).tolist(),'mean_speed_m_s':speed.tolist(),'steady_inventory_kg':(col*area/3).tolist(),'effective_capture_residence_s':(4*radius_m/(3*g*speed)).tolist()}

def main():
 d=json.loads((G/'inputs.json').read_text());caps=json.loads((G/'pellet-interfaces.json').read_text());base=json.loads((G/'chamber-dt150.0-q128.json').read_text());route=json.loads((G/'chamber-dt150.0-q128-cap0.1_eta100000.0.json').read_text());R=d['R_m'];area=4*math.pi*R*R
 Q=np.array([d['plant']['recovered_H_kg_s'],d['plant']['helium_per_s']*d['atomic_mass_kg']['he4']]);mass=np.array([2*d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']]);unit=steady_column(R,1000,Q,mass,[1,1]);ideal=steady_column(R,1000,route['cavity_source_mass_rate_atomic_kg_s'],mass,[.01,.01]);checks=[]
 def check(name,ok,value=None):
  checks.append({'name':name,'passed':bool(ok),'value':value});assert ok,(name,value)
 check('full source steady column matches frozen',math.isclose(unit['total_column_kg_m2'],sum(base['stationary'][0]['unit_capture_columns_kg_m2']),rel_tol=1e-12))
 # Reuse complete frozen scenarios atg.01/.1. A single linear response atg.03
 # uses the already-archived actual arrival bins, without new arrival extraction.
 profiles={g:next(x for x in base['cases'] if x['model']=='well_mixed_volume_balance' and x['T_K']==1000 and x['H_atoms_per_gas_particle']==2 and x['material_release_time_s']==0 and x['effective_wall_capture_g']==g) for g in [.01,.1]}
 trace=np.load(G/'chamber-dt150.0-q128.npz');t=trace['time_s'];source=trace['source_H_He_kg'];tau=np.array(unit['effective_capture_residence_s']);M=np.array([chamber.uniform_response(t,source[j],Q[j],.03/tau[j],0) for j in range(2)]);cols=3*M/area;press=cols*k*1000/(mass[:,None]*R)
 profiles[.03]=chamber.fields_stats(t,cols,press,d)
 ci=[x for x in caps['C13_cases'] if x['radius_m']==.02 and x['fill_fraction']==.8];limits=[]
 for x in ci:
  v=x['v_m_s'];rho=x['core_mass_per_capsule_volume_kg_m3'];a=x['radius_m'];beta=3/(4*rho*a);before=x['full_flight_IR_J_kg_core']+x['field_energy_J_kg_core'];budget=[h-before for h in x['enthalpy20to300_J_kg_core']];limit=[max(0,2*b/(beta*v*v)) for b in budget];added=[max(0,z-ideal['total_column_kg_m2']) for z in limit]
  limits.append({'mass_speed_quantile':x['mass_velocity_quantile'],'speed_m_s':v,'flight_s':x['flight_s'],'IR_plus_field_J_kg_core':before,'enthalpy_material_cases_J_kg_core':x['enthalpy20to300_J_kg_core'],'remaining_bulk_heat_budget_J_kg_core':budget,'zero_gas_satisfies_material_cases':[b>=0 for b in budget],'maximum_total_gas_column_kg_m2':limit,'maximum_added_column_above_ideal100nm_g0p01_kg_m2':added,'maximum_total_common_HHe_cavity_fraction_by_g':{str(g):[g*z/unit['total_column_kg_m2'] for z in limit] for g in [.01,.03,.1]},'maximum_added_common_HHe_cavity_fraction_atg0p01':[.01*z/unit['total_column_kg_m2'] for z in added],'maximum_added_H_only_cavity_fraction_atg0p01':[.01*z/unit['species_column_kg_m2'][0] for z in added],'maximum_added_He_only_cavity_fraction_atg0p01':[.01*z/unit['species_column_kg_m2'][1] for z in added],'uniform_drag_delay_at_total_limit_s_Cm1':[beta*z*x['flight_s']/2 for z in limit]})
 for L in limits:
  x=next(x for x in ci if x['mass_velocity_quantile']==L['mass_speed_quantile']);beta=3/(4*x['core_mass_per_capsule_volume_kg_m3']*x['radius_m']);frac={}
  for g in [.01,.03,.1]:
   e=next(e for e in profiles[g]['source_profile_pellet_encounters']['c13'] if e['source_mass_quantile']==L['mass_speed_quantile']);heat_per_unit_f=.5*beta*e['maximum_column_kg_m2']*x['v_m_s']**2
   frac[str(g)]=[max(0,b)/heat_per_unit_f for b in L['remaining_bulk_heat_budget_J_kg_core']]
  L['maximum_total_common_HHe_fraction_from_actual_chronology_by_g']=frac
 rows=[]
 for g in [.01,.03,.1]:
  encounters={x['source_mass_quantile']:x for x in profiles[g]['source_profile_pellet_encounters']['c13']}
  for f in [.001,.003,.01,.03,.1,1.]:
   ss=steady_column(R,1000,Q*f,mass,[g,g]);pellets=[]
   for x in ci:
    q=x['mass_velocity_quantile'];e=encounters[q];col=f*e['maximum_column_kg_m2'];beta=3/(4*x['core_mass_per_capsule_volume_kg_m3']*x['radius_m']);heat=.5*beta*col*x['v_m_s']**2;total=x['full_flight_IR_J_kg_core']+x['field_energy_J_kg_core']+heat
    pellets.append({'mass_speed_quantile':q,'maximum_scheduled_encounter_column_kg_m2':col,'gas_bulk_heat_J_kg_core_full_retention':heat,'total_IR_field_gas_J_kg_core':total,'within300K_material_cases':[total<=h for h in x['enthalpy20to300_J_kg_core']],'uniform_drag_delay_s_Cm1':beta*col*x['flight_s']/2,'near_launch_layer_delay_s_Cm1':beta*col*x['flight_s']})
   check(f'fraction scaling f{f} g{g}',math.isclose(ss['total_column_kg_m2'],unit['total_column_kg_m2']*f/g,rel_tol=1e-12))
   rows.append({'total_thermal_cavity_fraction_of_full_H_and_He_return_each':f,'effective_capture_H2_He':g,'steady':ss,'actual_dispatch_diffuse_peak_column_kg_m2':f*profiles[g]['peak_radial_column_kg_m2'],'actual_dispatch_diffuse_mean_column_kg_m2':f*profiles[g]['mean_radial_column_kg_m2'],'C13_2cm_capsules':pellets})
 # Atomic H sensitivity holds escaped H MASS fixed. No chemical adsorption or
 # recombination energy is silently used to remove atoms from the chamber.
 atomic=[]
 for gH in [.00797,.01,.03,.1]:
  a=steady_column(R,1000,route['cavity_source_mass_rate_atomic_kg_s'],[d['atomic_mass_kg']['p'],mass[1]],[gH,.01]);atomic.append({'chemical_form':'H+He','H_effective_capture':gH,'He_effective_capture':.01,**a})
 for L in limits:
  for h,col,budget in zip(L['enthalpy_material_cases_J_kg_core'],L['maximum_total_gas_column_kg_m2'],L['remaining_bulk_heat_budget_J_kg_core']):
   if budget>0:
    x=next(x for x in ci if x['mass_velocity_quantile']==L['mass_speed_quantile']);b=3/(4*x['core_mass_per_capsule_volume_kg_m3']*x['radius_m']);check('enthalpy at limit '+str(L['mass_speed_quantile']),abs((L['IR_plus_field_J_kg_core']+.5*b*col*L['speed_m_s']**2)/h-1)<1e-12)
 out={'schema':1,'family':d['family'],'scope':'Additional thermal-source robustness, not an inferred total retention probability. All full-source masses and chronology come from final350. No new hydro, receiver spectrum or plant optimization.','criterion':{'capsule_radius_m':.02,'fill_fraction':.8,'core_density_kg_m3':ci[0]['core_density_kg_m3'],'temperature_comparison_K':300,'enthalpy_interval_type':'Two measured-material/isotope-scaled cases; not a subjective credible interval.','heat_partition':'Conservative full incident gas kinetic energy percoremass. Local coating deposition, reflection, conduction, ablation and directional errors are separate.'},'full_return_atomic_kg_s':Q.tolist(),'unit_fullsource_g1':unit,'ideal_mobile100nm_routing_g0p01':ideal,'enthalpy_limits':limits,'total_source_fraction_scenarios':rows,'atomic_H_sensitivity':atomic,'population_bookkeeping':['The finite-cap PDE fraction describes an already-mobile ideal channel population. It does not include prompt energetic backscatter, damaged pores, seams, source exhaust, or all thermal reemission.','The tabulated f is the TOTAL fraction of each full returnedH/He stream that becomes thermal cavitygas. An additional source should be added to the existing source, not substituted for it.','Additional thermal sources add linearly throughSigma=4R/A sum_i[q_i/(g_i vbar_i)]. Species-specific sources andcapturecan be supplied to steady_column().','A reflected energetic ion can cross the cavity and reimplant in another retained wall. It is not automatically a1000K gas particle; only its later thermal liberation contributes toq_i, with the appropriate delay.','Gas collected or held in a closed process loop remains station inventory. Heating and ionization alone do not receive removal credit.','Theg=.01/.03/.1 rows are effective physical cleanup requirements, not demonstrated pump operation. The warm cleanup field/beam-tail interface is still qualified separately.','Nitrogen andhydrocarbon chemistry use the separate frozen chemistry packet and are not silently included inH2/He.'],'model_limits':['Steady uniform Maxwellian wallbalance is used even where gas collisions matter. Directed plumes or transport-limited localized pumps need spatial kinetics.','Chronological numbers assume instantaneous thermal release at the archived wall-arrival bins and a steady reservoir before the schedule starts. This is one explicit thermal-source history, not an upperbound on every possible plume.','Theg=.01/.1 chronology is rescaled from frozen cases. Oneg=.03 linear inventory response reuses the same archived arrivals; no fulltrace extraction or hydro is repeated.','The1st-percentile and slower2cm grains exceed the300K fullflightIR comparison without gas. Contained melting or protected delivery is an existing separate design branch.','Increasing grainradius affects source accelerationlength andhardware; these numbers do not preserve the previous stationrating automatically.'],'checks':checks,'source_sha256':{f:hashlib.sha256((G/f).read_bytes()).hexdigest() for f in ['inputs.json','pellet-interfaces.json','chamber-dt150.0-q128.json','chamber-dt150.0-q128.npz','chamber-dt150.0-q128-cap0.1_eta100000.0.json','chamber.py']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'limits':[x for x in limits if x['mass_speed_quantile'] in [.5,.999]],'checks_passed':len(checks)},indent=2))
if __name__=='__main__':main()
