from pathlib import Path
import sys,json
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path[:0]=[str(ROOT),str(ROOT/'b2-completion'),str(ROOT/'mixture-design')]
import radial_completion as r,updated_network_n15 as u,updated_network as c,pulse_model as pm
rng=np.random.default_rng(16);t=np.geomspace(.005,10,20)*1e9;rho=np.geomspace(1,3e5,20);y=rng.uniform(.0001,.1,(20,pm.NN));y/=y@pm.AA[:,None]
a,m=u.setup(mode='identity');ref=r.network(t,rho,y);got=u.kernel(t,rho,y,*a);err=max(float(np.max(abs(x-z)/np.maximum(abs(x),np.max(abs(x))*1e-15))) for x,z in zip(ref,got));assert err==0
node=balance=0.
d=json.loads((u.HERE/'nacreii-n15.json').read_text())
for curve in ['adopted','low','high']:
 a,m=u.setup(n15_curve=curve)
 for row in d['rows']:
  T=row['T9'];rr=c.evaluate_rates(T,a[0],a[1],pm.NR,*a[10:15]);old=rr.copy();out=u.apply_n15(rr,T,*a[15:],True);j,k=a[15]
  node=max(node,abs(out[j]/row[curve]-1))
  if old[k]>1e-290:balance=max(balance,abs((out[j]/old[j])/(out[k]/old[k])-1))
assert node<2e-12 and balance<2e-12
out={'identity_max_relative':err,'table_nodes_max_relative':node,'paired_balance_max_relative':balance,'metadata':m}
Path(__file__).with_name('n15-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='metadata'}))
