"""Extract archived sufficient statistics plus independent per-event energy arrays.
Downstream reweight.py needs only basis.json and a source spectrum, never raw runs.
"""
from pathlib import Path
import argparse,json,hashlib,gzip,statistics,math
from reweight import evaluate,positron_response_grid

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def extract(raw,out,spectra):
 out.mkdir(parents=True,exist_ok=True);(out/'events').mkdir(exist_ok=True);basis={'schema':1,'columns':{},'source_sha256':sha(raw/'transport_v3.py'),'scope':'Geant411.4.1 option4 EM with standard2gamma in-flight annihilation; natural carbon2000kg/m3,1000K, closed sphereR1e6m. Isotropic photon source25um below inner face; radial prompt and positron source1um inside cavity.10um production cuts,100eV electron tracking floor.','raw_inputs':{}};names=[]
 for S in [250,500,750,1000,2000]:
  names += [f'final-S{S}-isotropic',f'final-S{S}-prompt20']
 for S in [500,1000]:names += [f'final-S{S}-positron']
 for name in names:
  p=raw/(name+'.json');d=json.loads(p.read_text());assert d['code_sha256']==basis['source_sha256'];S=d['parameters']['sigma'];node=basis['columns'].setdefault(str(int(S)),{'sigma_kg_m2':S,'responses':[]});events={'source_file_sha256':sha(p),'parameters':d['parameters'],'columns':['deposited_MeV','escaped_calorimetric_MeV','vacuum_deposit_MeV','residual_inside_kinetic_MeV'],'records':[]}
  basis['raw_inputs'][name]={'path':str(p),'sha256':sha(p),'log_sha256':sha(raw/(name+'.log')),'source_sha256':d['code_sha256'],'parameters':d['parameters'],'creator_counts':d['secondary_creator_counts'],'energy_closure_policy':d['energy_closure_policy']}
  for r in d['records']:
   vals=[[e[k] for k in events['columns']] for e in r['events']];E=r['incident_calorimetric_energy_MeV'];leak=[v[1]/E for v in vals];res=[sum(v)-E for v in vals];N=len(vals);assert N==r['histories'];mean=statistics.mean(leak);se=statistics.stdev(leak)/math.sqrt(N);assert abs(mean-r['leakage_fraction'])<1e-14 and abs(se-r['sampling_standard_error'])<1e-14;assert all(len(v)==2 for e in r['events'] for v in e['annihilation_photons_by_parent'].values());assert max(abs(v) for v in res)<1e-4*E and statistics.mean(abs(v) for v in res)<1e-7*E
   evindex=len(events['records']);events['records'].append({'kinetic_energy_MeV':r['incident_energy_MeV'],'calorimetric_energy_MeV':E,'particle':r['incident_particle'],'values':vals})
   node['responses'].append({'particle':r['incident_particle'],'mode':d['parameters']['mode'],'energy_MeV':r['incident_energy_MeV'],'calorimetric_energy_MeV':E,'histories':N,'leakage_fraction':mean,'deposited_fraction':r['deposited_fraction'],'sampling_standard_error':se,'sampling_80pct_half_width':1.2815515655446004*se,'max_relative_event_energy_residual':max(abs(v) for v in res)/E,'mean_absolute_relative_energy_residual':statistics.mean(abs(v) for v in res)/E,'mean_signed_relative_energy_residual':statistics.mean(res)/E,'strict1e_7_exceedances':sum(abs(v)>1e-7*E for v in res),'residual_kinetic_MeV':sum(v[3] for v in vals),'escape_by_particle_MeV':r['escape_by_particle_MeV'],'event_archive':f'events/{name}.json.gz','event_archive_record':evindex,'raw_name':name})
  ep=out/'events'/f'{name}.json.gz'
  with ep.open('wb') as f:
   with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0,compresslevel=6) as g:g.write(json.dumps(events,separators=(',',':'),allow_nan=False).encode())
  basis['raw_inputs'][name]['event_archive_sha256']=sha(ep)
 write(out/'basis.json',basis)
 allresults={}
 for label,path in spectra.items():
  s=json.loads(path.read_text());write(out/(label+'-spectrum.json'),s);rows=[]
  for S in [250,500,750,1000,2000]:
   rows.append({'sigma_kg_m2':S,'scenarios':{f'{sc}:{pr}':evaluate(basis,s,S,sc,pr,hot_bypass_ppm=6.4) for sc in ['nominal_lines_hardbeta','broad_D1p10'] for pr in ['EM20','all_to_cold']}})
  allresults[label]={'spectrum_source_path':str(path),'spectrum_sha256':sha(path),'thermal_boundary_MeV_per_He':s['thermal_boundary_MeV_per_He'],'rows':rows,'positron_response_grids':[positron_response_grid(basis,s,S) for S in [500,1000]]}
 write(out/'weighted-results.json',allresults)
 return basis,allresults
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--current-spectrum',type=Path,required=True);ap.add_argument('--source-spectrum',type=Path,required=True);a=ap.parse_args();b,r=extract(a.raw,a.out,{'current128':a.current_spectrum,'source128':a.source_spectrum});print(json.dumps({k:[{'sigma':v['sigma_kg_m2'],'broad_loss':v['scenarios']['broad_D1p10:EM20']['leakage_fraction'],'broad_fullprompt':v['scenarios']['broad_D1p10:all_to_cold']['leakage_fraction'],'half80':v['scenarios']['broad_D1p10:EM20']['sampling_80pct_half_width']} for v in d['rows']] for k,d in r.items()},indent=2))
