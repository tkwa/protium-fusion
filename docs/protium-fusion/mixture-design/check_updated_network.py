from pathlib import Path
import sys,json,hashlib
import numpy as np
sys.path[:0]=['/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/b2-completion','/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion',str(Path(__file__).resolve().parent)]
import radial_completion as radial
import updated_network as u
import pulse_model as pm
rng=np.random.default_rng(38234)
t=np.geomspace(7e6,1e10,20);rho=np.geomspace(1,3e5,20);y=rng.uniform(.0001,.1,(20,pm.NN));y/=y@pm.AA[:,None]
args,meta=u.setup(mode='identity');ref=radial.network(t,rho,y);got=u.kernel(t,rho,y,*args)
err=max(float(np.max(abs(a-b)/(np.maximum(abs(a),np.max(abs(a))*1e-15)))) for a,b in zip(ref,got))
assert err<2e-12,err
args,meta=u.setup();d=json.loads((u.HERE/'updated-carbon-rates.json').read_text());maxnode=0.;maxbalance=0.
for row in d['tables'][u.PAIRS[0][0]]['rows']:
 T=row[0];r=u.evaluate_rates(T,args[0],args[1],pm.NR,*args[10:]);old=u.evaluate_rates(T,args[0],args[1],pm.NR,*args[10:-1],False)
 for k,(a,b) in enumerate(u.PAIRS):
  j,h=args[10][k];newtable=next(v[4] for v in d['tables'][a]['rows'] if v[0]==T)
  maxnode=max(maxnode,abs(r[j]/newtable-1))
  if old[h]>1e-290:maxbalance=max(maxbalance,abs((r[j]/old[j])/(r[h]/old[h])-1))
assert maxnode<2e-12 and maxbalance<2e-12
# Dense interpolation sensitivity; this is a numerical sensitivity, not experimental error.
lin,_=u.setup(interpolation='linear');maxinterp=0.
for T in np.geomspace(.3,2.,1000):
 r=u.evaluate_rates(T,args[0],args[1],pm.NR,*args[10:]);l=u.evaluate_rates(T,lin[0],lin[1],pm.NR,*lin[10:]);maxinterp=max(maxinterp,float(np.max(abs(r[args[10][:,0]]/l[args[10][:,0]]-1))))
out={'checks':{'identity_kernel_max_relative':err,'table_nodes_max_relative':maxnode,'paired_balance_max_relative':maxbalance,'pchip_vs_linear_rate_max_fraction_0p3_to_2_GK':maxinterp},'metadata':meta}
Path(__file__).with_name('updated-network-checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']))
