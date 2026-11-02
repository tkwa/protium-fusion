"""Radiation and nuclear first-order screens on fixed reversible source histories.
Radiation redistributes energy between shells and escapes at the outer boundary.
Entropy shifts are perturbative diagnostics; they are not a coupled new solution.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import simpson
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT),str(ROOT/'b2-completion'),str(ROOT/'mixture-design')]
from initial_screen import load,sha,quant
from feos_contact import FEOS
import radial_completion as radial
import pulse_model as pm
import updated_network
def integral_cutoff(t,values,cutoff):
 idx=t<cutoff;tt=np.r_[t[idx],cutoff];vv=np.r_[values[idx],np.interp(cutoff,t,values)];return float(simpson(vv,x=tt)) if len(tt)>2 else float(np.trapz(vv,x=tt))

rate=updated_network.install(radial,quantile='median');feos=FEOS();cases={}
for stage in ['B1','B2']:
 stem=Path('/home/tkwa/src/protium-fusion/output/assembly-audit')/f'{stage}-backward-feos128';d=np.load(stem.with_suffix('.npz'));meta=json.loads(stem.with_suffix('.json').read_text());t=d['time_s'];r=d['radius_m'];rho=d['density_kg_m3'];T=d['temperature_K'];m=d['mass_kg'];w=m/sum(m);n=len(m)
 p,i,b,*_=load(stage);parent=np.array(i['source_refinement']['parent_cell_for_each_new_cell']);src=json.loads(Path(i['source_refinement']['source_initial']).read_text());wm=np.array(src['mass_fraction_of_total']);active=np.isin(parent,np.where(np.cumsum(wm)-.5*wm>=.001)[0]);cvhot=np.array([(feos.state(stage,rho[0,j],T[0,j]*1.0001)['u_J_kg']-feos.state(stage,rho[0,j],T[0,j]*.9999)['u_J_kg'])/(.0002*T[0,j]) for j in range(n)])
 er=7.5657e-16*T**4;rc=.5*(r[:,1:]+r[:,:-1]);densface=.5*(rho[:,1:]+rho[:,:-1]);eface=.5*(er[:,1:]+er[:,:-1]);gradient=np.diff(er,axis=1)/np.diff(rc,axis=1);radiation=[]
 for kappa in [.001,.01,.03]:
  R=np.abs(gradient)/(kappa*densface*np.maximum(eface,1e-300));lim=(2+R)/(6+3*R+R**2);fi=-299792458.*lim*gradient/(kappa*densface);fo=299792458.*er[:,-1]/(4+3*kappa*rho[:,-1]*(r[:,-1]-rc[:,-1]));lum=np.c_[np.zeros(len(t)),4*np.pi*r[:,1:-1]**2*fi,4*np.pi*r[:,-1]**2*fo];qrad=np.diff(lum,axis=1)/m
  out=simpson(lum[:,-1],x=t);s_signed=simpson(-qrad/T,x=t,axis=0);s_abs=simpson(np.abs(qrad)/T,x=t,axis=0);eqfrac=s_abs/cvhot
  face_abs=simpson(np.sum(np.abs(lum[:,1:-1]),axis=1),x=t)
  events={e['cell']:e for e in meta['contact_events']};crossed=np.array([j in events for j in range(n)]);partial=np.array([integral_cutoff(t,np.abs(qrad[:,j])/T[:,j],events[j]['reverse_time_s'])/cvhot[j] if j in events else 0. for j in range(n)]);signedpartial=np.array([integral_cutoff(t,-qrad[:,j]/T[:,j],events[j]['reverse_time_s'])/cvhot[j] if j in events else 0. for j in range(n)])
  radiation.append({'postcontact_segment':{'scope':'Each parcel is integrated only from its first backward contact-density crossing to the hot snapshot. Fluxes still use the auxiliary whole-fluid neighbors; contact boundary geometry must be resolved for a physical trajectory.','crossed_mass_fraction':float(sum(w[crossed])),'max_abs_equivalent_hot_dT_over_T':float(max(partial[crossed])),'mass_quantiles_abs_equivalent_hot_dT_over_T':quant(partial[crossed],w[crossed]),'mass_quantiles_signed_equivalent_hot_dT_over_T':quant(signedpartial[crossed],w[crossed]),'total_mass_above_1pct':float(sum(w[crossed&(partial>.01)])),'total_mass_above_10pct':float(sum(w[crossed&(partial>.1)]))},'opacity_m2_kg':kappa,'outer_escaped_J':float(out),'outer_escaped_fraction_prepared':float(out/b['initial_energy_J']),'sum_absolute_internal_face_transfer_J':float(face_abs),'internal_face_transfer_fraction_prepared':float(face_abs/b['initial_energy_J']),'internal_transfer_scope':'Spatial redistribution counted by crossed interfaces; not additional external energy loss.','max_equivalent_hot_abs_temperature_fraction_all':float(max(eqfrac)),'max_equivalent_hot_abs_temperature_fraction_nonseed':float(max(eqfrac[active])),'mass_quantiles_equivalent_hot_abs_temperature_fraction_nonseed':quant(eqfrac[active],w[active]),'nonseed_mass_fraction_above_1pct_equivalent_hot_T':float(sum(w[active&(eqfrac>.01)])),'nonseed_mass_fraction_above_10pct_equivalent_hot_T':float(sum(w[active&(eqfrac>.1)])),'first_order_outer_escape_downsample_relative':float(simpson(lum[::2,-1],x=t[::2])/out-1),'cells':[{'cell':j,'mass_fraction':float(w[j]),'seed':bool(not active[j]),'integrated_forward_entropy_change_J_kg_K':float(s_signed[j]),'integrated_absolute_entropy_change_J_kg_K':float(s_abs[j]),'equivalent_hot_abs_dT_over_T':float(eqfrac[j])} for j in range(n)]})
 y=np.array([i['initial_Y'].get(name,0.) for name in pm.SPECIES]);yp=np.broadcast_to(y,(len(t)*n,len(y))).copy();dy,q,qn,qc=radial.network(np.maximum(T.ravel(),7e6),rho.ravel()/1000,yp);dy=dy.reshape(len(t),n,len(y));q=q.reshape(len(t),n)*1e-4;qn=qn.reshape(len(t),n)*1e-4;cat='n15' if stage=='B1' else 'c12';ci=pm.SPECIES.index(cat);catloss=simpson(np.abs(dy[:,:,ci]),x=t,axis=0)/y[ci];heat=simpson(q-qn,x=t,axis=0)
 cases[stage]={'source_trace_sha256':sha(stem.with_suffix('.npz')),'source_record_sha256':sha(stem.with_suffix('.json')),'radiation':radiation,'nuclear_fixed_composition':{'maximum_catalyst_fraction_changed':float(max(catloss)),'massmean_catalyst_fraction_changed':float(w@catloss),'nuclear_deposition_J':float(m@heat),'nuclear_deposition_fraction_prepared':float(m@heat/b['initial_energy_J']),'temperature_floor_K':7e6,'scope':'Median-carbon network on fixed initial composition and reverse-isentropic histories, with a conservative 7 MK temperature floor for cold tails below the updated carbon table. This is a perturbative diagnostic, not a full reactive trajectory or a formal low-temperature nuclear-data bound.'},'radiation_inventory_max_fraction_prepared':float(max(np.sum(er/rho*m,axis=1))/b['initial_energy_J'])}
 print(stage,{k:v for k,v in cases[stage].items() if k!='radiation'},flush=True)
 print([{k:v for k,v in row.items() if k!='cells'} for row in radiation],flush=True)
out={'scope':__doc__,'code_sha256':sha(__file__),'rate_update':rate,'cases':cases}
(P/'omission-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
