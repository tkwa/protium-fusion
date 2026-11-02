"""Frozen v3: real option4 EM, standard two-photon in-flight annihilation.
The replaced 11.4.1 three-photon branch has a separately preserved energy defect.
"""
from pathlib import Path
import subprocess,os,json,concurrent.futures,time,re,hashlib
P=Path(__file__).resolve().parent;script=P/'transport_v3.py';frozen_sha=hashlib.sha256(script.read_bytes()).hexdigest();env=os.environ.copy();env.update(PYTHONFAULTHANDLER='1',PYTHONDONTWRITEBYTECODE='1',PYTHONPYCACHEPREFIX=str(P/'cache'))
jobs=[]
for i,S in enumerate([250,500,750,1000,2000]):
 jobs.append(dict(name=f'final-S{S}-isotropic',sigma=S,mode='isotropic',energies='.511,2.312593,2.5438523,5.14304,5.657344',histories=20000,seed=123751+7919*i,depth_um=25))
 jobs.append(dict(name=f'final-S{S}-prompt20',sigma=S,mode='radial',energies='20',histories=5000,seed=298123+7919*i,depth_um=0))
for i,S in enumerate([500,1000]):jobs.append(dict(name=f'final-S{S}-positron',sigma=S,mode='radial',particle='e+',energies='.5,1,1.8,4.12',histories=10000,seed=971523+11111*i,depth_um=0))
jobs.append(dict(name='final-S500-cut1',sigma=500,mode='isotropic',energies='2.5438523,5.657344',histories=25000,seed=383337,depth_um=25,cut_um=1))
jobs.append(dict(name='final-S500-R1e4',sigma=500,mode='isotropic',energies='2.5438523,5.657344',histories=25000,seed=383337,depth_um=25,radius_m=10000))
jobs.append(dict(name='final-S500-surface',sigma=500,mode='isotropic',energies='2.5438523,5.657344',histories=25000,seed=383337,depth_um=0))
def run(spec):
 spec=spec.copy();name=spec.pop('name');out=P/(name+'.json');args={'radius_m':1000000,'cut_um':10}|spec;cmd=[str(P/'venv/bin/python'),'-u',str(script)]
 for k,v in args.items():cmd+=['--'+k.replace('_','-'),str(v)]
 cmd+=['--out',str(out)];t=time.monotonic()
 with (P/(name+'.log')).open('w') as f:proc=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1200)
 result={'name':name,'returncode':proc.returncode,'seconds':time.monotonic()-t,'output':str(out),'source_sha256':frozen_sha,'command':cmd}
 if proc.returncode==0:
  log=(P/(name+'.log')).read_text();tab=log.split('Material : ReceiverGraphite')[-1];cut=args['cut_um'];assert re.search(rf'Range cuts.*gamma\s+{cut} um\s+e-\s+{cut} um\s+e\+\s+{cut} um',tab);assert 'COMPLETE' in log and not any(x in log for x in ['GeomNav1002','GeomNav0003','FatalException','Fatal Python'])
  d=json.loads(out.read_text());assert d['code_sha256']==frozen_sha;assert all(r['histories']==args['histories'] and r['residual_inside_kinetic_MeV']==0 for r in d['records']);assert all(len(v)==2 for r in d['records'] for e in r['events'] for v in e['annihilation_photons_by_parent'].values());result['actual_cut_table_verified']=True;result['all_annihilation_multiplicities_two']=True;result['sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
 else:result['failure_preserved']=True
 (P/(name+'.status.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,jobs))
(P/'final-campaign-status.json').write_text(json.dumps(results,indent=2)+'\n')
assert all(r['returncode']==0 for r in results)
