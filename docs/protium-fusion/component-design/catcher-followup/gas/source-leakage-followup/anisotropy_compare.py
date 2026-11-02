"""Conditional eta1e4 front routing plus additional plenum/thermal leakage."""
from pathlib import Path
import json,hashlib
import numpy as np
from calculate import steady_column
P=Path(__file__).resolve().parent;G=P.parent
d=json.loads((G/'inputs.json').read_text());cap=json.loads((G/'pellet-interfaces.json').read_text());a=json.loads((P/'anisotropy-inputs.json').read_text());chosen=max(a['row']['resolutions'],key=lambda r:r['depth_resolution_parameter']);f0=np.array([chosen['flow_weighted_front_exit_fraction'][s] for s in ['p','he4']]);Q=np.array([d['plant']['recovered_H_kg_s'],d['plant']['helium_per_s']*d['atomic_mass_kg']['he4']]);mass=[2*d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']];plenum=.001
f1=f0+(1-f0)*plenum;base=steady_column(d['R_m'],1000,Q*f1,mass,[.01,.01]);slope=steady_column(d['R_m'],1000,Q*(1-f1),mass,[.01,.01])['total_column_kg_m2'];rows=[];checks=[]
for delta in [0,.001,.003,.004]:
 f=f1+(1-f1)*delta;s=steady_column(d['R_m'],1000,Q*f,mass,[.01,.01]);pellets=[]
 for x in cap['C13_cases']:
  if x['fill_fraction']!=.8 or x['radius_m'] not in [.01,.02] or x['mass_velocity_quantile'] not in [.5,.999]:continue
  beta=3/(4*x['core_mass_per_capsule_volume_kg_m3']*x['radius_m']);before=x['full_flight_IR_J_kg_core']+x['field_energy_J_kg_core'];heat=.5*beta*s['total_column_kg_m2']*x['v_m_s']**2
  maxcol=[max(0,2*(h-before)/(beta*x['v_m_s']**2)) for h in x['enthalpy20to300_J_kg_core']]
  pellets.append({'radius_m':x['radius_m'],'mass_speed_quantile':x['mass_velocity_quantile'],'IR_plus_field_J_kg_core':before,'gas_heat_J_kg_core':heat,'total_energy_J_kg_core':before+heat,'within300K_material_cases':[before+heat<=h for h in x['enthalpy20to300_J_kg_core']],'additional_post_plenum_leak_fraction_limit_from_delta0':[max(0,(c-base['total_column_kg_m2'])/slope) for c in maxcol]})
 rows.append({'additional_fraction_of_post_plenum_retained_gas_delta':delta,'total_thermal_cavity_fractions_H_He':f.tolist(),'steady':s,'capsules':pellets})
 checks.append(abs(s['total_column_kg_m2']-(base['total_column_kg_m2']+slope*delta))/s['total_column_kg_m2']<1e-12)
assert all(checks)
out={'criterion':'Same1000K,globalg=.01,300K material comparison as the leakage follow-up. No minimum anisotropy or total retention probability is inferred.','source_geometry':{k:a['row'][k] for k in ['cap_um','pitch_um','active_depth_um','anisotropy','maximum_absolute_refinement_difference']},'mean_ideal_mobile_front_fractions_H_He':f0.tolist(),'plenum_leak_fraction_of_locally_retained_gas':plenum,'definition':'f_total_i=f_front_i+(1−f_front_i)*f_plenum+(1−f_front_i)*(1−f_plenum)*delta. Thus a .001 plenum leak acts on the remaining population, not on gas already in the cavity. delta is an additional conditional thermal-release path.','cases':rows,'scope':['Only flow-weighted steady gas columns are evaluated for this eta1e4 comparison; the packet-specific eta1e4 transient is not integrated.','The .003 additional leakage case passes the median and99.9percent-speed2cm bulk comparisons. It does not qualify coating, plumes, fastreflection, gas chemistry, or the global cleanup apparatus.','The samecase is more restrictive for1cm because the fullflight IR dose is larger.','Eta1e4 suffices for this conditional heat screen; neither eta1e4 nor eta1e5 is established as a minimum material requirement. Performance must be assessed by total species-resolved cavity source, effective cleanup, and capsule/delivery acceptance.'],'checks_passed':len(checks),'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['anisotropy-inputs.json','calculate.py']},'frozen_source_sha256':{f:hashlib.sha256((G/f).read_bytes()).hexdigest() for f in ['inputs.json','pellet-interfaces.json']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'anisotropy-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows,indent=2))
