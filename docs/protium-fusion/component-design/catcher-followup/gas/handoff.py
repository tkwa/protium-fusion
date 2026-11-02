"""Named, non-overlapping conditional gas loads for parent plant integration."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
files=['joint-ports.json','regeneration.json','pumping-interfaces.json','routed-pellets.json','collision-screen.json','source-size.json']
j,r,p,g,c,s=[json.loads((P/f).read_text()) for f in files]
interface=j['representative_interface'];reg=r['representative']
getter=next(x for x in p['getter_regeneration_scenarios'] if x['working_H_mass_per_getter_mass']==.001 and x['cycle_s']==3600)
aux={'He_200K_view_refrigeration':interface['He_refrigeration_W_m2_etaCarnot0p5_at200K_sink'],
 'He_hot_photon_leakage_1e-4_refrigeration':next(x['additional_refrigeration_W_m2_etaCarnot0p5'] for x in interface['hot_photon_leakage_into5K_sensitivity'] if x['fraction_of_admitted1000K_aperture_energy_reaching5K']==1e-4),
 'He_5to85K_cycle_CA1000_3600s':reg['electric_cooldown_work_W_per_wall_m2_at_etaCarnot0p5'],
 'H_getter_sensible_regeneration_if_electrically_heated':getter['sensible_regeneration_W_m2'],
 'gas_300K_isothermal_compression_to1bar_50percent':p['gas_thermodynamic_scales']['isothermal300K_compression_to1bar'][0]['electrical_W_m2_at_free_energy_efficiency0p5'],
 'He_gas200to5K_refrigeration_50percent':2*p['gas_thermodynamic_scales']['He_refrigeration200to5K_reversible_W_m2']}
selected=next(x for x in g['cases'] if x['routing_case']=='cap0.1_eta100000.0' and x['global_capture_g']==.01)
fridge_names=['He_200K_view_refrigeration','He_hot_photon_leakage_1e-4_refrigeration','He_5to85K_cycle_CA1000_3600s','He_gas200to5K_refrigeration_50percent']
efficiency_cases=[]
for eta in [.15,.3,.5]:
 scaled={key:value*(.5/eta if key in fridge_names else 1.) for key,value in aux.items()}
 efficiency_cases.append({'fraction_of_Carnot':eta,'named_local_electrical_W_m2':scaled,'sum_named_local_electrical_W_m2':sum(scaled.values())})
out={'status':'conditional_component_interface','chronology':'Existing final350 dispatch and He rate unchanged. No new lower-rate or warmer-cold-side schedule is silently substituted.',
 'local_pump':{'joint_port':interface,'heat_bypass_hot_to_cold_W_m2':interface['hot_to_cold_bypass_upper_W_m2'],'named_electrical_W_m2':aux,'sum_named_local_electrical_W_m2':sum(aux.values()),'regeneration':reg,'getter_selected_capacity_scenario':getter,
 'hardware_mass_status':'Active-area and two-bank mass scales are provided, but total manifolds/refrigerator/control/repair mass is not engineered. A chosen service allocation is an acceptance target, not a result.',
 'additional_unassigned_requirements':['loaded-sorber equilibrium and mixture poisoning','actual baffle molecular transmission and hot photon leakage','solid supports/conduction and valve heat leaks','refrigerator installed mass and efficiency','gas chemical separation/adsorption work','front damage, permeation and repairs']},
 'cryogenic_efficiency_uncertainty':{'central_fraction_of_Carnot':.3,'subjective_80percent_interval':[.15,.5],'criterion':'Electrical work to provide the specified5K loads with200K ultimate rejection, conditional on large industrial helium refrigeration and the separately stated optical/cycle loads.','basis':'CERN-ATS-2013-004 official indexed primary report gives230W/W at4.5K,28.5percentofCarnot. The central value and interval are research judgments for this future architecture, not a measured5K/200K apparatus or a consequence of manufacturing efficiency.','source':'https://cds.cern.ch/record/1507611/files/CERN-ATS-2013-004.pdf','access':'Official indexed fulltext passage checked; direct repository fetch reported blocked by parent.','cases':efficiency_cases},
 'residual_cavity':{'selected_ideal_mobile_routing':selected,'global_capture_g':.01,'cleanup_credit':'Must come from a separate physical global cleanup device. The retained gas barrier does not provide this capture. The warm-ion packet currently owns that device and its loads.','warm_cleanup_electrical_and_mass_included_here':False},
 'grain_size_changes':s['source_mass_scaling_sensitivities'],
 'energy_rules':['Debit the single hot-to-cold bypass from engine input and credit cold rejection once.','Do not count200K-to5K radiation as a second independent heat source; add refrigerator electrical work once.','He bank regeneration heat is drawn from the200K waste reservoir and later returned; only cooling work is added. H getter sensible heating is conservatively charged electrically in this named example.','The local named electrical sum excludes global chamber cleanup, repair and unquantified component losses. It is not a complete installed apparatus power guarantee.','Changes to grain size affect source hardware and kinematics separately; orbital support guides are not scaled.'],
 'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'plant-handoff.json').write_text(json.dumps(out,indent=2)+'\n');print('local named W/m2',sum(aux.values()))
