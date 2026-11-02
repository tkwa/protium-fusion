"""Read saved research outputs and verify source identity and cross-ledger scales."""
from pathlib import Path
import json,hashlib,math
P=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();checks=[]
def check(name,condition):checks.append({'name':name,'pass':bool(condition)})
for stem,code in [('screen','screen'),('capillary','capillary'),('placement','placement'),('cathode','cathode'),('coil-array','coil_array'),('switching','switching')]:
 x=json.loads((P/(stem+'.json')).read_text());check(stem+' code hash',x['code_sha256']==sha(P/(code+'.py')))
 for key in ['input_sha256','sources_sha256','source_sha256']:
  v=x.get(key)
  if isinstance(v,dict):
   for path,h in v.items():
    pp=Path(path)
    if pp.is_absolute():check(stem+' source '+str(pp),pp.is_file() and sha(pp)==h)
co=json.loads((P/'coil-array.json').read_text());sw=json.loads((P/'switching.json').read_text());ca=json.loads((P/'cathode.json').read_text());pl=json.loads((P/'placement.json').read_text());e=1.602176634e-19;area=4*math.pi*(1.1e9)**2
check('named base field matches access model',co['parameters']['Bmin_front_T']==pl['selected_actual_duty']['B_T'])
for r in co['rows']:
 if r['status']!='finite_fixed_point':continue
 check('coil mass '+str((r['pitch_m'],r['Je_A_mm2'],r['shield_radial_column_kg_m2'])),abs((r['cable_with_semicircular_returns_kg_m2']+r['graphite_halo_kg_m2']+r['PE_buffer_kg_m2'])/r['subtotal_cable_shield_buffer_kg_m2']-1)<1e-12)
 check('cold COP '+str(r['pitch_m']),abs(r['photon_refrigeration_electric_W_m2']*co['parameters']['COP65to200']/r['mean_photon_heat65K_W_m2']-1)<1e-12)
 check('finite radiation view '+str(r['pitch_m']),0<r['outward_isotropic_outer_intercept_fraction']<=1)
for r in ca['rows']:
 if 'acceleration_V'in r:
  check('routing partition '+str((r['incident_ion'],r['acceleration_V'],r['Dparallel_over_Dnormal'])),abs(r['mobile_collection_fraction_at_mean_projected_depth']+r['mobile_front_fraction_at_mean_projected_depth']-1)<1e-14)
  check('fragment acceleration energy',abs(r['fragment_energy_MeV']*1e6*r['fragments_per_incident_ion']/r['acceleration_V']-1)<1e-14)
for r in ca['residual_load_scenarios']:
 check(r['species']+' voltage work',abs(r['300kV_one_pass_electric_W_m2']/r['one_pass_current_A_m2']/300000-1)<1e-14)
for r in sw['rows']:
 check('two-way field cycle '+str(r['pitch_m']),abs(r['unsigned_charge_and_discharge_throughput_W_m2']/(2*r['field_energy_J_m2']*sw['event_rate_s'])-1)<1e-14)
 for rr in r['ramp_ratings']:
  check('duty affects target',rr['g_with_original_cloud_unchanged']<.01 and rr['required_Bmin_T_to_restore_g001_in_same_rotation_screen']>co['parameters']['Bmin_front_T'])
check('full benchmark has less live time',sw['live_full_resolved_with10s_recovery']<sw['live_mass999_with10s_recovery'])
check('main-mass duty exact dispatch equality',abs(sw['live_mass999_with10s_recovery']-pl['selected_actual_duty']['duty'])<1e-14)
out={'scope':__doc__,'code_sha256':sha(__file__),'checks':checks,'number':len(checks),'passed':sum(v['pass']for v in checks),'independent_review':'pulse_transport independently verified gas/capillary normalizations and the periodic current-array field and energy coefficient; communicated review, not a separate archived test suite.','limits':'Passing accounting/source checks does not qualify field/cloud stability, cathode reflection/straggling, hot barrier behavior, actual magnet cable or pulse thermal access.'}
(P/'validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='checks'}));assert all(v['pass']for v in checks),[v for v in checks if not v['pass']]
