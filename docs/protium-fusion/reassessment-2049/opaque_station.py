"""Conditional complete small-pulse station allocation; no hardware validation.

Inputs are prescribed hot inflows. A B2 failed endpoint is permitted only when
its incompleteness is preserved in the output. All installed factory mass is
counted in this single station; no future-customer allocation is made.
"""
import argparse, hashlib, json, math, sys
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--b2',default='scaled-hot-B2-64-R30-tight.json');args=ap.parse_args();args.root=args.root.resolve();args.inputs=args.inputs.resolve()
sys.path.insert(0,str(args.root));from cycle_ledger import ledger,nuc,MU,MEV
G=6.67430e-11;SB=5.670374419e-8;YEAR=365.25*86400;LSUN=.9*3.828e26;VENUS=4.8675e24;AU=1.495978707e11
names={'n14':('scaled-hot-A-64-R30.json',1),'n15':('scaled-hot-B1-64-R30.json',1),'c12':(args.b2,2),'c13':('scaled-hot-C13-64-R30.json',1)}
files={cat:(str(args.inputs/name),protons) for cat,(name,protons) in names.items()};factory=args.inputs/'scaled-hot-He-64-R3.json'
records={cat:json.loads(Path(file).read_text()) for cat,(file,_) in files.items()};records['factory']=json.loads(factory.read_text())
source_status={k:{'source':(names[k][0] if k!='factory' else factory.name),'success':v['success'],'termination':v.get('effective_termination_reason',v.get('termination_reason')),'mass_kg':v['mass_kg'],'initial_energy_J':v['initial_energy_J'],'full_duration_elapsed_sound_times':v.get('elapsed_sound_times'),'monotone_final_velocity':all(b>a for a,b in zip(v['profile']['velocity_m_s'],v['profile']['velocity_m_s'][1:]))} for k,v in records.items()}
# Some original failed endpoints have a misleading final_time field: never
# promote an unsuccessful integration on the strength of that field.
for s in source_status.values():
 if not s['success']:s['effective_termination']='solver_failure'

def case(loss=1e-5,eta_heat=.4,mass=1e24,eta_driver=.6,prepared_multiplier=1.):
 d=ledger(collection_loss=loss,files=files,factory_file=str(factory))
 H=d['available_heat_MeV_per_helium']+d['escaped_radiation_MeV_per_helium'];U=d['fuel_preparation_MeV_per_helium']*prepared_multiplier
 capture=.99;aux=.1;R=1e9;A=4*math.pi*R**2;speed=1e4;guide_loss=.001;storage=1e5;rad_T=330.;rad_eps=.9
 # Closed station also captures radiation that left the modeled plasma;
 # apply collection loss once to all non-neutrino nuclear and prepared heat.
 captured=capture*(H+U);net=eta_heat*captured-U/eta_driver;q=captured-U
 replacement=mass*1e-4*1e8/YEAR
 stream_E=G*mass**2/(4*R);stream_mass=2*stream_E/speed**2;orbit=2*math.pi*R/speed;support=guide_loss*stream_E/orbit
 P=mass*.001
 if net<=0:return {'positive_cycle':False,'collection_loss_per_stage':loss,'eta_heat':eta_heat,'prepared_multiplier':prepared_multiplier,'H_MeV_He':H,'U_MeV_He':U,'preaux_net_MeV_He':net}
 rate=(P+replacement+support)/((1-aux)*net*MEV)
 gross=eta_heat*captured*MEV*rate;drive=U/eta_driver*MEV*rate;auxW=aux*net*MEV*rate
 heat=q*MEV*rate-P
 rad_flux=rad_eps*SB*rad_T**4;rad_area=heat/rad_flux
 maxE=max(x['initial_energy_J'] for x in records.values())/eta_driver*prepared_multiplier
 periods={}
 circulating_kg_s=0
 for st in d['stage_data']:
  r=records[st['cat']];a=sum(nuc(k).A*v for k,v in st['fuel_input'].items())
  per_pulse=r['mass_kg']/(a*MU);attempts=d['stage_attempts_per_net_helium'][st['cat']]*rate
  periods[st['cat']]=per_pulse/attempts;circulating_kg_s+=a*MU*attempts
 r=records['factory'];C=r['mass_kg']*r['output_mass_fractions']['c12']/(12*MU)
 periods['factory']=C/(rate*d['new_carbon_per_net_helium'])
 circulating_kg_s+=r['mass_kg']/periods['factory']
 # Dense explicit allowances; the reserve is installed project capacity,
 # rather than a claim that these engineering masses have been derived.
 budget={'absorber_10000kg_m2':A*1e4,'insulation_100kg_m2':A*100,
 'driver_guidance_and_peak_power_hardware_allowance':2.5e23,
 'isotope_separation_decay_recovery_inventory_allowance':1e23,
 'structural_guides_controls_redundancy_allowance':1e23,
 'thermal_manifolds_insulating_bridges_and_local_exchangers_allowance':1e23,
 'pulse_storage_two_max_pulses_at100kJ_kg':2*maxE/storage,
 'support_streams':stream_mass,'gross_power_conversion_0.1kg_W':.1*gross,
 'average_driver_conditioning_0.1kg_W':.1*drive,'radiator_panels20kg_m2':20*rad_area,
 'replacement_manufacture_1kg_W':replacement,'ten_year_replacement_stock':mass*1e-4*10}
 base=sum(budget.values());budget['unassigned_installed_design_reserve']=mass-base
 # Source-schedule similarity only; source storage supplies peak power.
 old=json.loads((args.root/'phase-two/driver/delivery-screen.json').read_text())['cases']
 lookup={('A',128):'n14',('B1',128):'n15',('B2',128):'c12',('C13',64):'c13',('He',64):'factory'}
 peak={}
 for o in old:
  k=lookup.get((o['stage'],o['zones']))
  if k and (o['stage']!='He' or o['tag']=='-r30'):
   peak[k]=o['launch_schedule']['electrical_peak_at_assumed_30pct_conversion_W']*(records[k]['mass_kg']/o['mass_kg'])*(2e10/R)*(.3/eta_driver)*prepared_multiplier
 # Local radial momentum delivery and practical refurbishment scale.
 rates={k:records[k]['mass_kg']/periods[k] for k in periods}
 return {'positive_cycle':True,'conditional_on_prescribed_inflow_and_B2_completion':True,
 'parameters':{'mass_kg':mass,'collector_radius_m':R,'specific_power_W_kg':.001,'driver_efficiency':eta_driver,'thermal_efficiency':eta_heat,'capture':capture,'loss_per_catalyst_stage':loss,'aux_fraction':aux,'replacement_fraction_per_year':1e-4,'replacement_energy_J_kg':1e8,'radiator_temperature_K':rad_T,'radiator_emissivity':rad_eps,'prepared_multiplier':prepared_multiplier},
 'H_MeV_He':H,'plasma_escaped_prompt_included_in_H_MeV_He':d['escaped_radiation_MeV_per_helium'],'U_MeV_He':U,'factory_U_MeV_He':d['factory_input_MeV_per_helium'],'preaux_net_MeV_He':net,'nuclear_to_prepared_H_over_U':H/U,
 'mass_budget_kg':budget,'mass_sum_kg':sum(budget.values()),'allocated_before_reserve_kg':base,
 'power':{'useful_W':P,'gross_generation_W':gross,'driver_W':drive,'other_auxiliary_W':auxW,'support_W':support,'replacement_W':replacement,'waste_W':heat,'captured_heat_W':captured*MEV*rate,'energy_closure_relative':(gross-drive-auxW-support-replacement-P)/P},
 'operations':{'helium_atoms_per_s':rate,'hydrogen_kg_s':rate*d['hydrogen_atoms_per_net_helium']*MU,'new_carbon_kg_s':rate*d['new_carbon_per_net_helium']*12*MU,'all_stage_material_throughput_kg_s':circulating_kg_s,'pulse_intervals_s':periods,'mean_interval_all_pulses_s':1/sum(1/t for t in periods.values()),'stage_material_flow_kg_s':rates,'maximum_electric_pulse_J':maxE,'source_peak_similarity_W':peak,'max_source_peak_per_driver_allowance_W_kg':max(peak.values())/2.5e23,'storage_peak_specific_power_W_kg':max(peak.values())/budget['pulse_storage_two_max_pulses_at100kJ_kg'],'support_stream_energy_J':stream_E,'radiator_area_m2':rad_area,'radiator_area_fraction_of_outer_sphere':rad_area/A,'mean_captured_heat_flux_W_m2':captured*MEV*rate/A,'mean_waste_flux_W_m2':heat/A,'radiator_emission_W_m2':rad_flux,'carbon_factory_all_hardware_counted_here':True},
 'ledger_residuals':{k:d[k] for k in ['baryon_residual_per_hydrogen','mass_defect_energy_relative_residual']}}

cases=[case(),case(eta_heat=.45),case(loss=1e-4),case(loss=1e-4,eta_heat=.45),case(prepared_multiplier=1.25),case(eta_driver=.5),case(mass=1.5e24)]
# Logistics includes transport separately from 100MJ/kg fabrication.
M=1e24;V=3e4;D=(5-.723)*AU;K=.5*M*V**2
out={'scope':'Conditional whole-installed-station requirements and selected-history arithmetic, not demonstrated hardware or a future forecast. B2 status remains explicit. The complete small helium factory is counted, with no shared-future-customer mass discount.','input_status':source_status,'input_sha256':{Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [x[0] for x in files.values()]+[str(factory)]},'cases':cases,'logistics':{'mass_kg':M,'fraction_of_Venus':M/VENUS,'destination_AU':5,'cruise_m_s':V,'Venus_source_AU':.723,'outer_carbon_source_AU':20,'outer_carbon_straight_cruise_years':(20-5)*AU/V/YEAR,'carbon_mass_reference_kg':4*math.pi*(1e9)**2*1e4,'route_accounting':'Total transport kinetic/electric energy includes all installed mass at the stated speed once, irrespective of source distance. Venus and outer-carbon routes have distinct cruise durations. Carbon lift/processing energy is additional.','straight_cruise_years':D/V/YEAR,'cruise_kinetic_J':K,'launch_electric_J_at50pct':K/.5,'launch_full_swarm_years_at50pct':K/.5/LSUN/YEAR,'construction_energy_J_at100MJ_kg':M*1e8,'construction_full_swarm_days':M*1e8/LSUN/86400,'additional_braking_energy': 'Mechanical braking must remove the cruise kinetic energy; reuse is possible but not credited to construction electricity. Receiver and transit geometry remain allocated engineering.'},'validation':{}}
checks={'mass_closure':all(abs(c['mass_sum_kg']/c['parameters']['mass_kg']-1)<1e-12 for c in cases if c['positive_cycle']),'electric_closure':all(abs(c['power']['energy_closure_relative'])<1e-12 for c in cases if c['positive_cycle']),'baseline_positive_reserve':cases[0]['mass_budget_kg']['unassigned_installed_design_reserve']>0,'baseline_radiator_fits':cases[0]['operations']['radiator_area_fraction_of_outer_sphere']<1,'baseline_all_pulse_periods_positive':all(x>0 for x in cases[0]['operations']['pulse_intervals_s'].values()),'failed_inputs_not_silenced':all(s['success'] or s.get('effective_termination')=='solver_failure' for s in source_status.values())}
assert all(checks.values()),checks
out['validation']={k:bool(v) for k,v in checks.items()};checks=out['validation']
args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'baseline':cases[0],'logistics':out['logistics'],'validation':checks},indent=2))
