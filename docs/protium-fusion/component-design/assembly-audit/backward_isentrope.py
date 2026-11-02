"""Reversible nonnuclear spherical screen from exact source hot profiles.

Fixed entropy per cell; FEOS pressure on tabulated isentropes. No radiation,
viscosity, nuclear energy, grain breakup or contact injection. Reverse time is
represented by reversed velocity and outward evolution. This tests a smooth
fluid continuation, not upstream cold launch or physical reverse shocks.
"""
from pathlib import Path
import argparse,json,hashlib,math,time
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from initial_screen import load,sha,quant,P,ROOT
from feos_contact import FEOS

class Adiabats:
 def __init__(self,stage,rho,T,m,feos,nodes=601):
  self.stage=stage;self.feos=feos;self.m=m;self.w=m/sum(m)
  initial=[feos.state(stage,a,b) for a,b in zip(rho,T)]
  self.initial=initial;entropy=np.array([a['s_J_kg_K'] for a in initial]);self.groups=[];self.ids=np.zeros(len(m),int)
  # Nested cells with the same entropy share a table; the seed keeps its actual
  # entropy. Group tolerance is numerical and is explicitly reported.
  for j,s in enumerate(entropy):
   match=[k for k,g in enumerate(self.groups) if abs(s/g['entropy']-1)<1e-10]
   if match:self.ids[j]=match[0]
   else:self.ids[j]=len(self.groups);self.groups.append({'entropy':float(s),'cells':[]})
   self.groups[self.ids[j]]['cells'].append(j)
  self.x=np.linspace(math.log(10.),math.log(max(rho)*1.3),nodes);self.lo=math.exp(self.x[0]);self.hi=math.exp(self.x[-1]);self.gl,self.gw=leggauss(8)
  for g in self.groups:
   states=[]
   for x in self.x:
    den=math.exp(x)
    target=g['entropy']
    def residual(lt):return feos.state(stage,den,math.exp(lt))['s_J_kg_K']-target
    anchor=g['cells'][0];guess=T[anchor]*(den/rho[anchor])**.6;left=math.log(max(10.,guess/2));right=math.log(min(3e9,guess*2))
    for iteration in range(25):
     fl=residual(left);fr=residual(right)
     if fl<=0<=fr:break
     if fl>0:left=max(math.log(10.),left-math.log(2))
     if fr<0:right=min(math.log(3e9),right+math.log(2))
    else:raise ValueError('Isentrope temperature bracket failed')
    temp=math.exp(brentq(residual,left,right,xtol=1e-11))
    st=feos.state(stage,den,temp)
    if st['p_Pa']<=0:raise ValueError('FEOS isentrope has nonpositive pressure')
    states.append(st)
   pp=np.array([z['p_Pa'] for z in states]);uu=np.array([z['u_J_kg'] for z in states]);tt=np.array([z['T_K'] for z in states])
   g['lp']=PchipInterpolator(self.x,np.log(pp));g['lT']=PchipInterpolator(self.x,np.log(tt));g['uFEOS']=PchipInterpolator(self.x,uu)
   dx=np.diff(self.x);mid=.5*(self.x[:-1]+self.x[1:]);gx=mid[:,None]+.5*dx[:,None]*self.gl
   du=.5*dx*np.sum(self.gw*np.exp(g['lp'](gx)-gx),axis=1)
   g['uwork']=np.r_[0.,np.cumsum(du)];g['T_range_K']=[float(min(tt)),float(max(tt))];g['pressure_monotone']=bool(all(np.diff(pp)>0));g['max_entropy_residual']=float(max(abs(z['s_J_kg_K']/g['entropy']-1) for z in states));g['work_vs_FEOS_delta_u_relative']=float(max(abs((g['uwork']-g['uwork'][0])-(uu-uu[0])))/max(abs(uu[-1]-uu[0]),1.))
  self.offset=np.array([a['u_J_kg'] for a in initial])-self.internal_from_work(rho,offset=False)
 def pressure(self,rho):
  if min(rho)<self.lo or max(rho)>self.hi:raise ValueError('Isentrope density domain exceeded')
  out=np.empty(len(rho));x=np.log(rho)
  for k,g in enumerate(self.groups):idx=self.ids==k;out[idx]=np.exp(g['lp'](x[idx]))
  return out
 def temperature(self,rho):
  out=np.empty(len(rho));x=np.log(rho)
  for k,g in enumerate(self.groups):idx=self.ids==k;out[idx]=np.exp(g['lT'](x[idx]))
  return out
 def sound(self,rho):
  out=np.empty(len(rho));x=np.log(rho)
  for k,g in enumerate(self.groups):idx=self.ids==k;out[idx]=np.sqrt(np.maximum(np.exp(g['lp'](x[idx]))/rho[idx]*g['lp'].derivative()(x[idx]),0.))
  return out
 def internal_from_work(self,rho,offset=True):
  out=np.empty(len(rho));x=np.log(rho)
  for k,g in enumerate(self.groups):
   idx=self.ids==k;xx=x[idx];j=np.clip(np.searchsorted(self.x,xx)-1,0,len(self.x)-2);left=self.x[j];dx=xx-left;gx=.5*(xx+left)[:,None]+.5*dx[:,None]*self.gl
   out[idx]=g['uwork'][j]+.5*dx*np.sum(self.gw*np.exp(g['lp'](gx)-gx),axis=1)
  return out+self.offset if offset else out
 def dump(self):
  return {'density_range_kg_m3':[self.lo,self.hi],'nodes_per_adiabat':len(self.x),'entropy_group_relative_tolerance':1e-10,'entropy_groups':[{'entropy_J_kg_K':g['entropy'],'cells':g['cells'],'T_range_K':g['T_range_K'],'pressure_monotone':g['pressure_monotone'],'max_entropy_residual':g['max_entropy_residual'],'work_vs_FEOS_delta_u_relative':g['work_vs_FEOS_delta_u_relative']} for g in self.groups],'internal_energy_scope':'Each cell is anchored to its exact initial FEOS u. Changes integrate the same interpolated P/rho² used by momentum, making work exchange conservative. FEOS energy-difference discrepancy measures EOS/table consistency.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['B1','B2'],required=True);ap.add_argument('--duration',type=float,required=True);ap.add_argument('--subdivide',type=int,choices=[1,2],default=1);ap.add_argument('--nodes',type=int,default=601);ap.add_argument('--max-step',type=float,default=2e-7);ap.add_argument('--rtol',type=float,default=1e-9);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
 start=time.monotonic();profile,inp,base,r,v,m,T,rho,ye,yi=load(a.stage);refinement={'subdivide':a.subdivide}
 if a.subdivide==2:
  oldm=m.copy();oldr=r.copy();oldv=v.copy();oldrho=rho.copy();M=sum(m);q=np.r_[0.,np.cumsum(m)/M];q[-1]=1.;mid=(.5*(q[:-1]**(1/3)+q[1:]**(1/3)))**3;newq=np.sort(np.r_[q,mid]);r=np.interp(newq,q,oldr**3)**(1/3);m=np.diff(newq)*M;T=np.repeat(T,2);v=np.interp(r,oldr,oldv);rho=m/(4*np.pi/3*np.diff(r**3));oldnm=np.r_[.5*(oldm[:-1]+oldm[1:]),.5*oldm[-1]];newnm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];scale=math.sqrt((oldnm@(oldv[1:]**2))/(newnm@(v[1:]**2)));v*=scale;inp=json.loads(json.dumps(inp));inp['source_refinement']['parent_cell_for_each_new_cell']=np.repeat(inp['source_refinement']['parent_cell_for_each_new_cell'],2).tolist();refinement.update(velocity_scale=scale,max_relative_density_error=float(max(abs(rho/np.repeat(oldrho,2)-1))))
 n=len(m);M=sum(m);w=m/M;nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];feos=FEOS();ad=Adiabats(a.stage,rho,T,m,feos,a.nodes)
 Rscale=r[-1];Vscale=max(abs(v));tscale=Rscale/Vscale
 def physical(z):
  rr=np.r_[0.,z[:n]*Rscale];vv=np.r_[0.,z[n:]*Vscale]
  if min(np.diff(rr))<=0:raise ValueError('Shell crossing in trial state')
  den=m/(4*np.pi/3*np.diff(rr**3));return rr,vv,den
 def rhs(t,z):
  rr,vv,den=physical(z);pp=ad.pressure(den);force=4*np.pi*rr[1:]**2*(pp-np.r_[pp[1:],0.]);accel=force/nm
  return np.r_[vv[1:]*tscale/Rscale,accel*tscale/Vscale]
 z0=np.r_[r[1:]/Rscale,-v[1:]/Vscale]
 sol=solve_ivp(rhs,(0,a.duration/tscale),z0,method='DOP853',rtol=a.rtol,atol=a.rtol*1e-2,max_step=a.max_step/tscale,dense_output=True)
 ts=np.linspace(0,sol.t[-1],801);zz=sol.sol(ts);rs=[];vs=[];ds=[];ps=[];temps=[];us=[];kin=[];intern=[];divs=[];charmin=[]
 for z in zz.T:
  rr,vv,den=physical(z);pp=ad.pressure(den);uu=ad.internal_from_work(den);cs=ad.sound(den);rc=.5*(rr[1:]+rr[:-1]);vc=.5*(vv[1:]+vv[:-1]);dv=4*np.pi*np.diff(rr*rr*vv)*den/m
  rs.append(rr);vs.append(vv);ds.append(den);ps.append(pp);temps.append(ad.temperature(den));us.append(uu);kin.append(.5*nm@(vv[1:]**2));intern.append(m@uu);divs.append(dv);charmin.append([min(np.diff(vc+sgn*cs)/np.diff(rc)) for sgn in [-1,1]])
 rs=np.array(rs);vs=np.array(vs);ds=np.array(ds);ps=np.array(ps);temps=np.array(temps);us=np.array(us);ener=np.array(kin)+np.array(intern);initialE=ener[0];divs=np.array(divs)
 contact=inp['contact'];target=contact['rho_kg_m3'];parent=np.array(inp['source_refinement']['parent_cell_for_each_new_cell']);srcprofile=json.loads(Path(inp['source_refinement']['source_initial']).read_text());srcmass=np.array(srcprofile['mass_fraction_of_total']);active_parent=np.where(np.cumsum(srcmass)-.5*srcmass>=.001)[0];active=np.isin(parent,active_parent)
 events=[]
 for j in np.where(active)[0]:
  residual=ds[:,j]-target;changes=np.where(residual[:-1]*residual[1:]<=0)[0]
  for k in changes:
   root=brentq(lambda tt:physical(sol.sol(tt))[2][j]-target,ts[k],ts[k+1],xtol=1e-13);rr,vv,den=physical(sol.sol(root));tmp=ad.temperature(den)
   events.append({'cell':int(j),'mass_fraction':float(w[j]),'reverse_time_s':root*tscale,'radius_center_m':.5*(rr[j]+rr[j+1]),'backward_velocity_center_m_s':.5*(vv[j]+vv[j+1]),'T_K':float(tmp[j]),'T_ratio_to_FEOS_contact':float(tmp[j]/contact['T_K']),'pressure_Pa':float(ad.pressure(den)[j]),'crossing_direction_density':-1 if residual[k]>0 else 1})
 unique={e['cell']:e for e in events};ew=np.array([e['mass_fraction'] for e in unique.values()]);et=np.array([e['reverse_time_s'] for e in unique.values()]);er=np.array([e['radius_center_m'] for e in unique.values()]);beyond=np.array([j in unique for j in range(n)])
 # This numerical reversal checks the reversible ODE only; it is not evidence
 # that a finite-grain source can deliver its boundary conditions.
 back=solve_ivp(rhs,(sol.t[-1],0),sol.y[:,-1],method='DOP853',rtol=a.rtol,atol=a.rtol*1e-2,max_step=a.max_step/tscale)
 reverse_error=back.y[:,-1]-z0
 np.savez_compressed(a.out.with_suffix('.npz'),time_s=ts*tscale,radius_m=rs,reverse_velocity_m_s=vs,density_kg_m3=ds,temperature_K=temps,pressure_Pa=ps,internal_J_kg=us,mass_kg=m,kinetic_J=np.array(kin),internal_J=np.array(intern),reverse_divergence_s=divs,initial_z=z0,final_z=sol.y[:,-1],initial_entropy_J_kg_K=np.array([x['s_J_kg_K'] for x in ad.initial]))
 result={'scope':__doc__,'stage':a.stage,'source_initial_sha256':sha(Path(inp['source_refinement']['source_initial'])),'source128_initial_sha256':sha(Path('/home/tkwa/src/protium-fusion/output/source128')/f'{a.stage}-H2-source128-carbon-median.initial.json'),'code_sha256':sha(__file__),'FEOS_library_sha256':sha(feos.library),'adiabats':ad.dump(),'profile_refinement':refinement,'numerics':{'zones':n,'rtol':a.rtol,'atol':a.rtol*1e-2,'max_step_s':a.max_step,'samples':len(ts),'requested_duration_s':a.duration,'elapsed_s':sol.t[-1]*tscale,'success':bool(sol.success),'message':sol.message,'nfev':sol.nfev,'viscosity':0.,'radiation':False,'nuclear':False},'energy':{'initial_total_J':initialE,'initial_internal_J':intern[0],'initial_kinetic_J':kin[0],'initial_FEOS_energy_relative_to_source_FastEOS':float(initialE/base['initial_energy_J']-1),'max_total_energy_residual_fraction':float(max(abs(ener/initialE-1))),'final_internal_J':intern[-1],'final_kinetic_J':kin[-1],'pressure_work_internal_to_kinetic_J':float(kin[-1]-kin[0])},'density_minmax_kg_m3':[float(ds.min()),float(ds.max())],'temperature_minmax_K':[float(temps.min()),float(temps.max())],'minimum_shell_width_m':float(np.diff(rs,axis=1).min()),'minimum_width_relative_to_initial':float((np.diff(rs,axis=1)/np.diff(r)).min()),'reverse_compression_max_mass_fraction':float(max(np.sum((divs<0)*w,axis=1))),'reverse_char_gradient_min_per_s':np.min(charmin,axis=0).tolist(),'return_to_initial_max_dimensionless_error':float(max(abs(reverse_error))),'return_to_initial_radius_max_m':float(max(abs(reverse_error[:n]))*Rscale),'return_to_initial_velocity_max_m_s':float(max(abs(reverse_error[n:]))*Vscale),'contact_crossed_mass_fraction':float(sum(ew)),'contact_active_mass_fraction':float(sum(w[active])),'contact_not_crossed_active_mass_fraction':float(sum(w[active&~beyond])),'contact_crossing_time_mass_quantiles_s':quant(et,ew) if len(et) else None,'contact_crossing_radius_mass_quantiles_m':quant(er,ew) if len(er) else None,'contact_crossing_max_relative_T_error':float(max(abs(e['T_ratio_to_FEOS_contact']-1) for e in events)) if events else None,'contact_events':events,'interpretation':'A continuous spherical fluid with a vacuum exterior is extended backward with fixed entropy. After a parcel reaches its declared postshock contact state, continuing it as fluid is an auxiliary boundary construction; the original incoming grains/contact geometry are not modeled. Convergence, nonnegative width and reversibility do not establish a globally shock-free continuum solution at all resolutions.','wall_seconds':time.monotonic()-start,'trace_sha256':sha(a.out.with_suffix('.npz'))}
 a.out.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['contact_events','adiabats']},indent=2),flush=True)
if __name__=='__main__':main()
