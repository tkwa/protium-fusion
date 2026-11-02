from pathlib import Path
import sys,json,numpy as np
from scipy.integrate import simpson
P=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/component-design/assembly-audit');sys.path.insert(0,str(P))
from initial_screen import load,sha
from backward_isentrope import Adiabats
from feos_contact import FEOS
profile,inp,base,r,v,m,T,rho,*_=load('B2');ad=Adiabats('B2',rho,T,m,FEOS(),1201);ad.ids=np.repeat(ad.ids,2)
root=Path('/home/tkwa/src/protium-fusion/output/assembly-audit/B2-backward-feos256');d=np.load(root.with_suffix('.npz'));j=json.loads(root.with_suffix('.json').read_text());events={q['cell']:q['reverse_time_s'] for q in j['contact_events']};mass=d['mass_kg'];w=mass/sum(mass);active=np.array([q in events for q in range(len(mass))]);out={}
for sign in [-1,1]:
 minimum=1e99;negcells=set();allminimum=1e99;allnegcells=set();minwhere=None;negative_gradient_envelope=[];jump_mach=0.
 for n,t in enumerate(d['time_s']):
  r=d['radius_m'][n];v=d['reverse_velocity_m_s'][n];rho=d['density_kg_m3'][n];cs=ad.sound(rho);rc=.5*(r[1:]+r[:-1]);vc=.5*(v[1:]+v[:-1]);g=np.diff(vc+sign*cs)/np.diff(rc)
  allminimum=min(allminimum,float(min(g)));negative_gradient_envelope.append(max(0.,float(-min(g))));jump_mach=max(jump_mach,float(max(np.maximum(-np.diff(vc+sign*cs),0.)/(.5*(cs[1:]+cs[:-1])))))
  for k in np.where(g<0)[0]:allnegcells.update([int(k),int(k+1)])
  allowed=np.array([active[k] and active[k+1] and t<=min(events[k],events[k+1]) for k in range(len(g))]);good=g[allowed]
  if len(good) and min(good)<minimum:
   minimum=float(min(good));idx=np.where(allowed)[0][np.argmin(good)];minwhere={'time_s':float(t),'interface':int(idx),'enclosed_mass_fraction':float(sum(w[:idx+1]))}
  for k in np.where(allowed&(g<0))[0]:negcells.update([int(k),int(k+1)])
 out[str(sign)]={'minimum_gradient_full_per_s':allminimum,'minimum_gradient_postcontact_nonseed_per_s':minimum,'postcontact_min_location':minwhere,'mass_touching_negative_gradient_any_time':float(sum(w[list(allnegcells)])),'mass_touching_negative_postcontact_gradient':float(sum(w[list(negcells)])),'negative_postcontact_cells':sorted(negcells),'sampled_global_compressive_gradient_time_integral':float(simpson(negative_gradient_envelope,x=d['time_s'])),'maximum_compressive_characteristic_jump_over_local_cs':jump_mach}
record={'scope':'Resolved acoustic slopes u±c on nested256 FEOS continuation; restricted version includes only interfaces whose two neighboring nonseed parcels are between their contact-density crossings and hot snapshot. Finite negative slopes are a steepening diagnostic, not a proven shock.','code_sha256':sha(__file__),'trace_sha256':sha(root.with_suffix('.npz')),'results':out};(P/'characteristic-refinement.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
