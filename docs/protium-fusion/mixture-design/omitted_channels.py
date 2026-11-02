"""First-order omitted-channel exposure on saved trajectories; no feedback."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
sys.path[:0]=['/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion']
import pulse_model as pm
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
allrates=pm.LIB.get_rates();have=set(pm.NUCS)
miss=[r for r in allrates if all(n in have for n in r.reactants) and any(n not in have for n in r.products)]
print('Omitted entry channels',len(miss),flush=True)

def one(p):
 d=json.loads(p.read_text());z=np.load(d['trace_path']);m=z['cell_mass_kg'];ts=z['time_s'];st=z['state'].reshape(len(m),int(z['variables_per_zone']),len(ts));nv=d['state_variables']
 T=np.maximum(st[:,nv.index('T/T0')]*d['T0_K'],1e6);r=np.vstack([np.zeros(len(ts)),st[:,nv.index('r/R0')]*d['R0_m']]);rho=m[:,None]/(4*np.pi/3*np.diff(r**3,axis=0));Y=np.maximum(st[:,:pm.NN],0)
 dt=np.r_[.5*np.diff(ts)[0],.5*(np.diff(ts)[:-1]+np.diff(ts)[1:]),.5*np.diff(ts)[-1]];ncat=d['mass_kg']/pm.MU*sum(v for n,v in d['initial_Y'].items() if n!='p');weights=m[:,None]/pm.MU*dt/ncat
 t=T/1e9;tf=np.stack([np.ones_like(t),1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
 rows=[]
 for rr in miss:
  rate=sum(np.exp(np.clip(np.einsum('i,izt->zt',s.a,tf),-745,700)) for s in rr.sets)
  f=rate*rr.prefactor*(rho/1000)**rr.dens_exp
  for n in rr.reactants:f*=Y[:,pm.NUCS.index(n)]
  if rr.use_ye_weighting:f*=np.einsum('i,zit->zt',pm.ZZ,Y)
  count=float(np.sum(f*weights))
  if count>1e-16:rows.append({'name':rr.fname,'reaction':str(rr),'source_label':rr.src,'Q_MeV':rr.Q,'events_per_initial_catalyst_unperturbed':count,'omitted_products':[str(n).lower() for n in rr.products if n not in have]})
 # This transfer channel is absent from this REACLIB release. A constant
 # one-barn envelope above the atomic-mass threshold is a sensitivity,
 # not a measured cross section or a rigorous unitarity bound.
 from pynucastro import Nucleus
 eth=-(Nucleus('c13').mass+Nucleus('p').mass-Nucleus('c12').mass-Nucleus('d').mass)
 x=eth*pm.MEV/(1.380649e-23*T);mu=13/14*pm.MU
 sv=1e-28*np.sqrt(8*1.380649e-23*T/(np.pi*mu))*(1+x)*np.exp(-x)
 f=(rho/pm.MU)*sv*Y[:,pm.SPECIES.index('c13')]*Y[:,pm.SPECIES.index('p')]
 transfer=float(np.sum(f*weights))
 rows.sort(key=lambda x:x['events_per_initial_catalyst_unperturbed'],reverse=True)
 return {'case':p.stem,'source_sha256':sha(p),'trace_sha256':sha(d['trace_path']),'channels':rows,'sum_exposure':sum(x['events_per_initial_catalyst_unperturbed'] for x in rows),'C13_pd_C12_one_barn_above_threshold_exposure':transfer,'C13_pd_C12_threshold_cm_MeV':eth}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--record',type=Path,action='append',required=True);a.add_argument('--output',type=Path,required=True);a=a.parse_args()
 rows=[one(p) for p in a.record];out={'scope':'Entry flows from existing species into omitted species, evaluated on unmodified saved trajectories. Not a closed enlarged network or abundance prediction; threshold for escalation is appreciable exposure. Library coverage is not a guarantee of all physical channels.','code_sha256':sha(__file__),'omitted_entry_channel_count':len(miss),'cases':rows};a.output.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps([{'case':x['case'],'sum':x['sum_exposure'],'top':x['channels'][:6]} for x in rows],indent=2))
