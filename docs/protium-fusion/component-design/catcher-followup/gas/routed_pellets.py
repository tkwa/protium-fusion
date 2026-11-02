"""Apply packet-specific mobile front routing to capsule and timing margins."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
caps=json.loads((P/'pellet-interfaces.json').read_text());rows=[];sources=['pellet-interfaces.json']
for tag in ['cap0.1_eta100000.0','cap0.01_eta100000.0']:
 fn=f'chamber-dt150.0-q128-{tag}.json';sources.append(fn);d=json.loads((P/fn).read_text())
 for g in [.01,.1,1.]:
  r=next(x for x in d['cases'] if x['model']=='well_mixed_volume_balance' and x['T_K']==1000 and x['H_atoms_per_gas_particle']==2 and x['effective_wall_capture_g']==g and x['material_release_time_s']==0)
  arrivals={x['source_mass_quantile']:x for x in r['source_profile_pellet_encounters']['c13']};cr=[]
  for c in caps['C13_cases']:
   if c['fill_fraction']!=.8:continue
   a=c['radius_m'];rho=c['core_mass_per_capsule_volume_kg_m3'];q=c['mass_velocity_quantile'];e=arrivals[q];col=e['maximum_column_kg_m2'];v=c['v_m_s'];beta=3/(4*rho*a);heat=.5*beta*col*v*v;total=c['full_flight_IR_J_kg_core']+c['field_energy_J_kg_core']+heat
   cr.append({'radius_m':a,'fill_fraction':.8,'mass_velocity_quantile':q,'velocity_m_s':v,'worst_encounter_column_kg_m2':col,'gas_heat_J_kg_core_full_retention':heat,'IR_plus_field_J_kg_core':total-heat,'total_bulk_energy_J_kg_core':total,'enthalpy20to300_material_cases_J_kg':c['enthalpy20to300_J_kg_core'],'within300K_material_cases':[total<=h for h in c['enthalpy20to300_J_kg_core']],'uniform_arrival_delay_s_Cm1':beta*col*c['flight_s']/2,'near_launch_layer_arrival_delay_s_Cm1':beta*col*c['flight_s'],'coherent_crossflow3000m_s_displacement_m_Cm1':beta*col*3000*c['flight_s']})
  rows.append({'routing_case':tag,'global_capture_g':g,'front_source_fractions_H_He':d['cavity_source_fraction_H_He'],'mean_column_kg_m2':r['mean_radial_column_kg_m2'],'peak_column_kg_m2':r['peak_radial_column_kg_m2'],'mean_pressure_Pa':r['mean_pressure_Pa'],'peak_pressure_Pa':r['peak_pressure_Pa'],'C13_capsule_rows':cr})
out={'scope':'Exact packet-specific ideal already-mobile finite-cap front exits, actual final350 chronology, diffuse steady-initialized cavity, and actual C13CH2 capsule IR/field model. Added reflection, defects, permeability, backpressure and unknown population fractions remain separate. Globalg is conditional physical cleanup; the front cap does not supply it. Bulk energy and geometric momentum deposition do not qualify local skin ablation or directional plumes.','cases':rows,'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in sources},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'routed-pellets.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([{k:r[k] for k in ['routing_case','global_capture_g','mean_column_kg_m2','mean_pressure_Pa']} for r in rows],indent=2))
