"""Exact-profile assembly kinematics and EOS matching; no 3D trajectory claimed."""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path.insert(0,str(ROOT));from fast_eos import FastEOS
from feos_contact import FEOS
Z={'p':1,'n14':7,'n15':7,'c12':6,'c13':6}
STEMS={s:Path('/home/tkwa/src/protium-fusion/output/source128')/f'{s}-H{h}-source128-carbon-median' for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quant(x,w):
 order=np.argsort(x);return np.interp([.1,.5,.9],np.cumsum(w[order])/sum(w)-.5*w[order]/sum(w),np.array(x)[order]).tolist()
def load(stage):
 stem=STEMS[stage];p=json.loads(stem.with_suffix('.initial.json').read_text());i=json.loads(stem.with_suffix('.inputs.json').read_text());b=json.loads(Path(i['source_refinement']['source_record']).read_text())
 r=np.r_[0,p['outer_radius_m']];v=np.r_[0,p['velocity_m_s']];m=np.array(p['mass_fraction_of_total'])*b['mass_kg'];T=np.array(p['T_K']);rho=m/(4*np.pi/3*np.diff(r**3));ye=sum(Z[n]*y for n,y in i['initial_Y'].items());yi=sum(i['initial_Y'].values())
 return p,i,b,r,v,m,T,rho,ye,yi

def main():
 eos=FastEOS();feos=FEOS();cases={}
 for s in STEMS:
  profile,i,b,r,v,m,T,rho,ye,yi=load(s);M=sum(m);w=m/M;nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
  st=eos.state(rho/1000,T,ye,yi);pressure=st.p*.1;u=st.e*1e-4;cs=np.sqrt(st.gamma1*pressure/rho)
  fe=[feos.state(s,rr,tt) for rr,tt in zip(rho,T)];fp=np.array([x['p_Pa'] for x in fe]);fu=np.array([x['u_J_kg'] for x in fe]);fS=np.array([x['s_J_kg_K'] for x in fe])
  div=4*np.pi*np.diff(r*r*v)*rho/m;accel=4*np.pi*r[1:]**2*(pressure-np.r_[pressure[1:],0])/nm
  reverse=-v;center_r=.5*(r[1:]+r[:-1]);center_v=.5*(reverse[1:]+reverse[:-1]);chars={}
  for sign in [-1,1]:
   lam=center_v+sign*cs;dt=np.diff(lam);converging=dt<0
   chars[str(sign)]={'initial_converging_interface_count':int(sum(converging)),'minimum_frozen_characteristic_intersection_s':float(min(-np.diff(center_r)[converging]/dt[converging])) if any(converging) else None,'scope':'Local planar characteristic slopes at the snapshot, frozen in time. Pressure evolution and spherical terms change trajectories; this is a steepening flag, not a predicted shock time.'}
  result={'source_initial':str(STEMS[s].with_suffix('.initial.json')),'source_initial_sha256':sha(STEMS[s].with_suffix('.initial.json')),'source_inputs':str(STEMS[s].with_suffix('.inputs.json')),'source_inputs_sha256':sha(STEMS[s].with_suffix('.inputs.json')),'mass_kg':M,'zones':len(m),'hot_outer_radius_m':r[-1],'hot_density_minmax_kg_m3':[float(min(rho)),float(max(rho))],'hot_temperature_minmax_K':[float(min(T)),float(max(T))],'hot_forward_velocity_minmax_m_s':[float(min(v)),float(max(v))],'hot_forward_divergence_minmax_per_s':[float(min(div)),float(max(div))],'hot_forward_expanding_mass_fraction':float(sum(w[div>0])),'hot_subsonic_cellmass_fraction':float(sum(w[np.abs(center_v)<cs])),'hot_U_fraction_prepared':float(m@u/b['initial_energy_J']),'FEOS_hot_U_relative_to_FastEOS':float((m@fu)/(m@u)-1),'FEOS_hot_energy_delta_fraction_prepared':float(m@(fu-u)/b['initial_energy_J']),'FEOS_hot_pressure_ratio_mass_quantiles':quant(fp/pressure,w),'FEOS_hot_temperature_scope':'Temperature, density and composition held fixed; FEOS has its own cold energy reference and partial ionization. Hot energy difference is an EOS comparison, not an extra source debit.','ballistic_backward_shell_crossing':bool(any(np.diff(reverse)<0)),'ballistic_backward_divergence_positive_all_cells':bool(all(-div>0)),'initial_sound_characteristics':chars,'cells':[{'index':j,'mass_fraction':float(w[j]),'rho_kg_m3':float(rho[j]),'T_K':float(T[j]),'r_center_m':float(center_r[j]),'v_forward_center_m_s':float(-center_v[j]),'sound_m_s':float(cs[j]),'pressure_Pa':float(pressure[j]),'u_J_kg':float(u[j]),'FEOS':fe[j],'forward_divergence_per_s':float(div[j]),'boundary_pressure_accel_m_s2':float(accel[j])} for j in range(len(w))]}
  if s in ['B1','B2']:
   contact=i['contact'];rc=contact['rho_kg_m3'];tc=contact['T_K'];contact_fast=eos.state(rc/1000,tc,ye,yi);contact_fe=feos.state(s,rc,tc)
   parent=np.array(i['source_refinement']['parent_cell_for_each_new_cell']);srcprofile=json.loads(Path(i['source_refinement']['source_initial']).read_text());srcmass=np.array(srcprofile['mass_fraction_of_total']);active=np.isin(parent,np.where(np.cumsum(srcmass)-.5*srcmass>=.001)[0]);above=rho>rc;selected=active&above
   times=np.full(len(m),np.nan);positions=np.full(len(m),np.nan)
   for j in np.where(selected)[0]:
    def f(t):rr=r+t*reverse;return m[j]/(4*np.pi/3*(rr[j+1]**3-rr[j]**3))-rc
    times[j]=brentq(f,0,.1,xtol=1e-15);rr=r+times[j]*reverse;positions[j]=.5*(rr[j]+rr[j+1])
   result['constant_contact']={'source_contact':contact,'reproduced_FEOS':contact_fe,'relative_FEOS_contact_pressure_error':contact_fe['p_Pa']/contact['p']-1,'FastEOS_over_FEOS_pressure_at_contact':float(contact_fast.p*.1/contact_fe['p_Pa']),'FastEOS_minus_FEOS_u_J_kg_at_contact':float(contact_fast.e*1e-4-contact_fe['u_J_kg']),'hot_mass_below_contact_density':float(sum(w[active&~above])),'hot_mass_above_contact_density':float(sum(w[selected])),'contact_entropy_relative_residual_active_max':float(max(abs(fS[active]/contact['s']-1))),'ballistic_contact_crossing_time_s_quantiles':quant(times[selected],w[selected]),'ballistic_contact_crossing_radius_m_quantiles':quant(positions[selected],w[selected]),'ballistic_contact_minmax_s':[float(min(times[selected])),float(max(times[selected]))],'ballistic_contact_minmax_m':[float(min(positions[selected])),float(max(positions[selected]))],'constant_velocity_compression_work_deficit_J':float(sum(m[selected]*(u[selected]-contact_fast.e*1e-4))),'constant_velocity_compression_work_deficit_fraction_prepared':float(sum(m[selected]*(u[selected]-contact_fast.e*1e-4))/b['initial_energy_J']),'scope':'Ballistic reverse geometry alone has no pressure forces or pdV energy transfer. Its contact crossings are geometric scales. Below-contact-density shells have no crossing on this monotone ballistic branch; droplet expansion, staged contact or a prior density turning point requires separate dynamics.'}
  cases[s]=result
  print(json.dumps({s:{k:v for k,v in result.items() if k not in ['cells','source_inputs','source_initial']}}),flush=True)
 out={'scope':__doc__,'code_sha256':sha(__file__),'FastEOS_code_sha256':sha(ROOT/'fast_eos.py'),'FastEOS_table_sha256':sha(ROOT/'electron-eos.npz'),'FEOS_library_sha256':sha(feos.library),'FEOS_helper_sha256':sha(P/'feos_contact.py'),'FEOS_material_database_sha256':hashlib.sha256(feos.db.encode()).hexdigest(),'cases':cases}
 (P/'initial-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
