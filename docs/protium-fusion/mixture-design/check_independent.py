from pathlib import Path
import json,sys
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path[:0]=[str(ROOT),str(ROOT/'mixture-design')]
import pulse_model as pm,updated_network as c,updated_network_independent as u
rng=np.random.default_rng(351);temp=np.geomspace(7e6,1e10,20);rho=np.geomspace(1,3e5,20);y=rng.uniform(.0001,.1,(20,pm.NN));y/=y@pm.AA[:,None]
exact=True
for q in ['minus1sigma','median','plus1sigma']:
 ac,_=c.setup(quantile=q);au,_=u.setup(quantile={name:q for name,rev in u.PAIRS});rc=c.kernel(temp,rho,y,*ac);ru=u.kernel(temp,rho,y,*au);exact=exact and all(np.array_equal(x,z) for x,z in zip(rc,ru))
assert exact
source=json.loads((u.HERE/'updated-carbon-rates.json').read_text());err=bal=0.
for q12,q13 in [('minus1sigma','plus1sigma'),('plus1sigma','minus1sigma')]:
 qs=dict(zip([n for n,r in u.PAIRS],[q12,q13]));args,meta=u.setup(quantile=qs)
 for row in source['tables'][u.PAIRS[0][0]]['rows']:
  T=row[0];rr=u.evaluate_rates(T,args[0],args[1],pm.NR,*args[10:]);old=u.evaluate_rates(T,args[0],args[1],pm.NR,*args[10:-1],False)
  for k,(name,rev) in enumerate(u.PAIRS):
   j,h=args[10][k];target=next(r[source['columns'].index(qs[name])] for r in source['tables'][name]['rows'] if r[0]==T);err=max(err,abs(rr[j]/target-1))
   if old[h]>1e-290:bal=max(bal,abs((rr[j]/old[j])/(rr[h]/old[h])-1))
assert err<2e-12 and bal<2e-12
out={'co_shifted_kernel_bit_identical_to_current':exact,'independent_table_nodes_max_relative':err,'paired_equilibrium_ratio_max_relative':bal,'metadata':meta}
Path(__file__).with_name('independent-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='metadata'}))
