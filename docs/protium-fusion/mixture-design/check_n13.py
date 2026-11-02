from pathlib import Path
import sys,json
import numpy as np
sys.path[:0]=['/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/b2-completion','/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion',str(Path(__file__).resolve().parent)]
import radial_completion as radial
import updated_network_n13 as n
import pulse_model as pm
args,meta=n.setup(mode='identity');rng=np.random.default_rng(5432);T=np.geomspace(1e7,1e10,20);rho=np.ones(20)*1e5;y=rng.uniform(.001,.1,(20,pm.NN));y/=y@pm.AA[:,None]
a=radial.network(T,rho,y);b=n.kernel(T,rho,y,*args);err=max(float(np.max(abs(x-z)/np.maximum(abs(x),np.max(abs(x))*1e-15))) for x,z in zip(a,b));assert err<2e-12
args,meta=n.setup();co=args[-1];j,h=args[-2];rows=[]
for t in [.3,.5,.7,1,1.25,1.5,1.75,2,3,10]:
 tf=np.array([1,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)]);new=float(np.exp(co@tf).sum());old=pm.RATES[j].eval(t*1e9)
 rows.append({'T9':t,'Li2006_Eq17':new,'JINA_lg06':old,'ratio':new/old})
out={'checks':{'identity_kernel_max_relative':err},'rates':rows,'metadata':meta};Path(__file__).with_name('n13-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']))
