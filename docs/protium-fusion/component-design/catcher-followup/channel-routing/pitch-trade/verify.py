"""Independent pitch-trade accounting and frozen-input checks, no PDE rerun."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent;d=json.loads((P/'results.json').read_text());checks=[]
def check(name,x):checks.append({'name':name,'pass':bool(x)})
def near(a,b):return math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-18)
for name,h in d['inputs_sha256'].items():check('input '+name,hashlib.sha256(Path(name).read_bytes()).hexdigest()==h)
fine=json.loads((P.parent/'fine-candidates.json').read_text());check('frozen fine solver version',fine['solver_sha256']==hashlib.sha256((P.parent/'calculate.py').read_bytes()).hexdigest())
for s in d['exact_similarity']:
 check('eta/p2'+str(s['pitch_um']),near(s['eta_to_equal_20um_1e4']/s['pitch_um']**2,25)and near(s['eta_to_equal_20um_1e5']/s['pitch_um']**2,250))
src=json.loads((P.parents[1]/'stopping/surface-stopping.json').read_text())
for r in d['new_refinements']:
 for v in r['versions']:
  for species in ['p','he4']:
   flow={name:sum(x['species'][species]['packet_nuclei_m2'])/x['period_s']for name,x in src['stages'].items()};avg=sum(flow[name]*v['stages'][name][species]['front_fraction']for name in flow)/sum(flow.values());check('weighted front '+str(r['pitch_um'])+str(v['depth_resolution_parameter'])+species,near(avg,v['flow_weighted_front_fraction'][species]))
for r in d['geometry_and_pumping']:
 label=str(r['pitch_um'])+'/'+str(r['downstream_speed_m_s']);check('slit '+label,near(r['slot_fraction'],.1/r['pitch_um']))
 for i in [0,1]:
  C=1/(1/r['slit_C_m_s_H2_He'][i]+1/r['via_C_m_s_H2_He'][i]);S=1/(1/C+1/r['downstream_speed_m_s']);check('series '+label+str(i),near(C,r['channel_series_C_m_s_H2_He'][i])and near(S,r['full_series_speed_m_s_H2_He'][i]));check('open fraction '+label+str(i),near(r['equivalent_front_open_fraction_for_added_leak0p001_H2_He'][i]/r['slot_fraction'],r['per_slit_cap_area_equivalent_transmission_for_added_leak0p001_H2_He'][i]));check('relaxation '+label+str(i),near(r['local_plenum_plus_slit_void_m3_m2']/S,r['local_void_relaxation_s_H2_He'][i]))
 vf=math.pi*.0005**2/.01**2;mass=2000*((.175+r['added_physical_thickness_m_preserving350kg_m2'])*(1-vf)-50e-6*.75-150e-6*r['slot_fraction']);check('350kg/m2 preserved '+label,near(mass,350.));check('mean partial pressures '+label,near(sum(r['mean_source_partial_pressure_Pa_H2_He']),r['mean_source_total_pressure_Pa']))
for r in d['direct_cold_panel_comparison']:
 active=[10*math.sqrt(1000/300),13*math.sqrt(1000/273.15)]
 for i in [0,1]:check('cold port '+str(r['Sp_m_s'])+str(i),near(1/(1/r['port_C_H2_He_m_s'][i]+1/(active[i]*r['active_area_fraction_Hgetter_Hecryo'][i])),r['Sp_m_s']))
 check('cold refrigeration '+str(r['Sp_m_s']),near(r['He_refrigeration_W_m2_at0p5Carnot_200K_sink'],78*r['He_200K_to5K_heat_W_m2']))
check('all model checks',all(x['pass']for x in d['checks']));result={'scope':__doc__,'checks':checks,'number':len(checks),'passed':sum(c['pass']for c in checks)};(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result['passed'],'/',result['number']);assert all(c['pass']for c in checks)
files=[p for p in sorted(P.rglob('*'))if p.is_file()and p.name!='manifest.json'and '__pycache__'not in p.parts];(P/'manifest.json').write_text(json.dumps({'scope':'New channel-pitch trade only. No original archive hashes changed.','files':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files}},indent=2)+'\n')
