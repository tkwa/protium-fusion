from pathlib import Path
import concurrent.futures,subprocess,os,json,hashlib,time,re
P=Path(__file__).resolve().parent;original=P.parent/'coupled-photon-transport';runner=original/'transport_v3.py';expected='f1f95c9c10bb2bea9302639f22e6202bd522ca1206a0ef3e34e6d298510b853e';assert hashlib.sha256(runner.read_bytes()).hexdigest()==expected
env=os.environ.copy();env.update(PYTHONFAULTHANDLER='1',PYTHONDONTWRITEBYTECODE='1',PYTHONPYCACHEPREFIX=str(original/'cache'))
jobs=[]
for i,S in enumerate([200,350]):
 jobs += [(f'S{S}-broad',S,'isotropic','2.5438523,5.657344',20000,25,681023+i*1717),(f'S{S}-prompt20',S,'radial','20',5000,0,681923+i*1717)]
def run(job):
 name,S,mode,E,N,depth,seed=job;out=P/(name+'.json');cmd=[str(original/'venv/bin/python'),'-u',str(runner),'--data',str(original/'data'),'--sigma',str(S),'--mode',mode,'--energies',E,'--histories',str(N),'--depth-um',str(depth),'--cut-um','10','--radius-m','1000000','--seed',str(seed),'--out',str(out)];t=time.monotonic()
 with (P/(name+'.log')).open('w') as f:r=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 assert r.returncode==0,(name,r.returncode);d=json.loads(out.read_text());assert d['code_sha256']==expected;assert all(len(v)==2 for x in d['records'] for e in x['events'] for v in e['annihilation_photons_by_parent'].values());log=(P/(name+'.log')).read_text();assert re.search(r'Range cuts.*gamma\s+10 um\s+e-\s+10 um\s+e\+\s+10 um',log.split('Material : ReceiverGraphite')[-1]);z={'name':name,'command':cmd,'seconds':time.monotonic()-t,'output':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'source_sha256':expected,'actual_cuts_verified':True};print(json.dumps(z),flush=True);return z
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(run,jobs))
(P/'run-provenance.json').write_text(json.dumps({'rows':rows},indent=2)+'\n')
