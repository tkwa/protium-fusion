"""Pressure, expansion, radiation and isotope diagnostics from saved radial traces.

Run with the pinned NumPy/SciPy/pynucastro/Numba runtime used by the pulse model.
Default compares original64, optimized64, changed source64, and available128.
--manifest accepts [{label,family,result,trace?,reaction_flows?},...]. Paths can
be absolute or relative to docs/protium-fusion. New nuclear-rate variants should
set reaction_flows:false unless their exact rate arrays are supplied separately.
No hydro, archive edits, emissivity loss or new physics is introduced.
"""
from pathlib import Path
import argparse,json,gzip,hashlib,csv,sys,math
import numpy as np
from numba import njit
from scipy.integrate import trapezoid
import importlib.metadata
P=Path(__file__).resolve().parent;ROOT=P.parent;sys.path.insert(0,str(ROOT))
from fast_eos import FastEOS,KB_CGS,NA,ARAD_CGS
import pulse_model as pm
C=299792458.;MU=1.66053906660e-27;MEV=1.602176634e-13
A=pm.AA;Z=pm.ZZ;CO=pm.COEFF;RID=pm.RATE_INDEX;PREF=pm.PREF;DEXP=pm.DEXP;REACT=pm.REACT;NR=pm.NR;NN=pm.NN

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def array_sha(a):
 a=np.ascontiguousarray(a);h=hashlib.sha256();h.update(str(a.shape).encode());h.update(a.dtype.str.encode());h.update(a.tobytes());return h.hexdigest()
def resolve(p):
 p=Path(p);return p if p.is_absolute() else ROOT/p

def defaults():
 work=Path('/home/tkwa/src/protium-fusion/output/mixture-runs');rows=[]
 for stage,h in [('A','1'),('B1','1'),('B2','2'),('C13','1')]:rows.append(dict(label=f'{stage} original64',family='reference64',result=str(work/('control' if stage=='A' else 'grid')/f'{stage}-H{h}-fixed_energy-z64-k0p001-t1.json')))
 for stage,h,folder in [('A','3','refine'),('B1','2','grid'),('B2','2p5','refine'),('C13','4','grid')]:rows.append(dict(label=f'{stage} optimized64 H{h}',family='optimized64',result=str(work/folder/f'{stage}-H{h}-fixed_energy-z64-k0p001-t1.json')))
 for stage,name in [('A','AH3-N2000-local-entropy-fixed-energy-reference-floor'),('C13','C13H4-N2000-local-entropy-fixed-energy')]:rows.append(dict(label=f'{stage} source64',family='source64',result=f'component-design/launcher-mixture-yield/{name}.json',trace=f'/home/tkwa/src/protium-fusion/output/launcher-mixture-yield/{name}.npz'))
 for stage,h,fam in [('A','1','reference128'),('B1','1','reference128'),('B2','2','reference128'),('C13','1','reference128'),('A','3','optimized128'),('B1','2','optimized128'),('B1','1p5','candidate128'),('B2','2p5','optimized128'),('C13','3','candidate128'),('C13','4','optimized128')]:
  f=work/'mesh'/f'{stage}-H{h}-fixed_energy-z128-k0p001-t1.json'
  if f.exists():rows.append(dict(label=f'{stage} {fam} H{h}',family=fam,result=str(f)))
 return rows

@njit
def sampled_flows(temp,rho,y,weights):
 # Exact archived ReacLib arrays; unclipped signed polynomial matches solver.
 # Positive-inventory version only estimates physical gross reaction turnover.
 nt,nz=temp.shape;out=np.zeros((nt,NR));raw=np.zeros((nt,NR))
 for it in range(nt):
  for iz in range(nz):
   t=temp[it,iz]/1e9;tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)]);rates=np.zeros(NR)
   for s in range(len(RID)):
    ex=0.
    for k in range(7):ex+=CO[s,k]*tf[k]
    rates[RID[s]]+=np.exp(ex)
   for j in range(NR):
    f=rates[j]*PREF[j]*rho[it,iz]**DEXP[j];fr=f
    for k in range(3):
     n=REACT[j,k]
     if n>=0:f*=max(y[it,iz,n],0.);fr*=y[it,iz,n]
    out[it,j]+=weights[iz]*f;raw[it,j]+=weights[iz]*fr
 return out,raw

def event(values,time,kind='max'):
 i=int(np.argmax(values) if kind=='max' else np.argmin(values));return {'value':float(values[i]),'time_s':float(time[i]),'sample_index':i}
def wmean(x,w):return x@w

def reduce(entry,eos,do_flows):
 file=resolve(entry['result']);d=json.loads(file.read_text())
 if not d['success'] or d['termination_reason']!='final_time':raise ValueError(f'Incomplete case excluded: {file}')
 if d.get('reduced_three_species_A_network'):raise ValueError('Full27 network required')
 if d.get('coulomb_eos') or 'coulomb' in d.get('model','').lower():raise ValueError('Use the actual Coulomb EOS for this variant')
 trace=resolve(entry.get('trace',d.get('trace_path',str(file.with_suffix('.npz')))))
 z=np.load(trace);nt=len(z['time_s']);nz=int(z['zones']);nv=int(z['variables_per_zone']);names=z['species'].tolist();assert names==pm.SPECIES
 raw=z['state'].reshape(nz,nv,nt).transpose(2,0,1);y=raw[:,:,:NN];time=z['time_s'];m=z['cell_mass_kg'];M=float(sum(m));w=m/M;nodes=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
 T=raw[:,:,NN]*d['T0_K'];r=raw[:,:,NN+1]*d['R0_m'];v=raw[:,:,NN+2]*float(z['velocity_unit_cm_s'])/100.;ri=np.column_stack([np.zeros(nt),r[:,:-1]]);V=4*np.pi/3*(r**3-ri**3);rho=m[None,:]/V
 yp=np.maximum(y,0);norm=yp@A;ye=(yp@Z)/norm;yi=yp.sum(axis=2)/norm;state=eos.state(rho/1000,T,ye,yi);pressure=state.p*.1;pr=ARAD_CGS*.1*T**4/3;pi=rho*KB_CGS*NA*T*yi*1e-4;pe=pressure-pr-pi
 vi=np.column_stack([np.zeros(nt),v[:,:-1]]);Vdot=4*np.pi*(r*r*v-ri*ri*vi);compress=np.maximum(-(v-vi),0.);mode=d['viscosity_mode']
 if mode=='compression_cap':compress=np.minimum(compress,np.maximum(-Vdot*rho/m[None,:]*(r-ri),0.))
 elif mode=='volume_gate':compress=np.where(Vdot<0,compress,0.)
 elif mode=='divergence':compress=np.maximum(-Vdot*rho/m[None,:]*(r-ri)/3.,0.)
 else:raise ValueError('Unsupported or legacy artificial viscosity')
 sound_speed=np.sqrt(np.maximum(state.gamma1*pressure/rho,0.));qvis=rho*(d['artificial_viscosity_quadratic']*compress**2+d['artificial_viscosity_linear']*sound_speed*compress);viscous_heating=-(qvis*Vdot).sum(axis=1)
 U=(state.e*1e-4)@m;Ui=(1.5*KB_CGS*NA*T*yi*1e-4)@m;Ur=(ARAD_CGS*.1*T**4*V).sum(axis=1);Ue=U-Ui-Ur;K=.5*(v*v)@nodes
 factor=float(z['energy_unit_erg_g'])*1e-4;weak=raw[:,:,NN+3]@m*factor;escaped=raw[:,:,NN+4]@m*factor;nuclear=raw[:,:,NN+5]@m*factor;ynet=np.einsum('z,tzn->tn',w,y);cat={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}[d['stage']];cats=ynet[0,names.index(cat)]
 mass_q=np.r_[0.,np.cumsum(w)];radius50=np.array([np.cbrt(np.interp(.5,mass_q,np.r_[0.,x**3])) for x in r]);radius90=np.array([np.cbrt(np.interp(.9,mass_q,np.r_[0.,x**3])) for x in r]);volsum=V.sum(axis=1);Vweights=V/volsum[:,None]
 kappa=d['opacity_m2_kg'];tau=(kappa*rho*(r-ri)).sum(axis=1);L=4*np.pi*r[:,-1]**2*C*ARAD_CGS*.1*T[:,-1]**4/(4+3*kappa*rho[:,-1]*(r[:,-1]-ri[:,-1])/2)
 if not d['radiation_diffusion']:L=np.zeros_like(time)
 hist={'time_s':time,'temperature_max_K':T.max(axis=1),'temperature_massmean_K':T@w,'pressure_max_Pa':pressure.max(axis=1),'artificial_viscous_pressure_massmean_Pa':qvis@w,'artificial_viscous_pressure_max_Pa':qvis.max(axis=1),'pressure_plus_artificial_viscosity_massmean_Pa':(pressure+qvis)@w,'artificial_viscous_heating_W':viscous_heating,'pressure_massmean_Pa':pressure@w,'pressure_volumemean_Pa':(pressure*Vweights).sum(axis=1),'ion_pressure_massmean_Pa':pi@w,'electron_pressure_massmean_Pa':pe@w,'LTE_radiation_pressure_massmean_Pa':pr@w,'density_max_kg_m3':rho.max(axis=1),'density_massmean_kg_m3':rho@w,'density_volumemean_kg_m3':M/volsum,'surface_radius_m':r[:,-1],'outer_velocity_m_s':v[:,-1],'lumped_mass_rms_velocity_m_s':np.sqrt(2*K/M),'half_mass_radius_m':radius50,'ninety_percent_mass_radius_m':radius90,'kinetic_J':K,'internal_J':U,'ion_internal_J':Ui,'electron_positron_internal_J':Ue,'LTE_radiation_energy_J':Ur,'integrated_nuclear_J':nuclear,'hot_weak_escape_J':weak,'deposited_nuclear_J':nuclear-weak,'net_escaped_radiation_J':escaped,'surface_radiative_luminosity_W':L,'grey_center_to_surface_optical_depth':tau,'radiation_reservoir_over_luminosity_s':[float(u/l) if l>0 else None for u,l in zip(Ur,L)],'cumulative_alpha_per_initial_catalyst':(ynet[:,names.index('he4')]-ynet[0,names.index('he4')])/cats,'instantaneous_F_and_heavier_nuclei_per_initial_catalyst':ynet[:,Z>=9].sum(axis=1)/cats}
 extrema={k:event(hist[k],time,'min' if 'radius' in k else 'max') for k in ['pressure_max_Pa','pressure_massmean_Pa','pressure_volumemean_Pa','density_max_kg_m3','density_massmean_kg_m3','density_volumemean_kg_m3','surface_radius_m','half_mass_radius_m','ninety_percent_mass_radius_m','LTE_radiation_energy_J','surface_radiative_luminosity_W','grey_center_to_surface_optical_depth','artificial_viscous_pressure_massmean_Pa','artificial_viscous_pressure_max_Pa','artificial_viscous_heating_W']}
 peak_it,peak_iz=np.unravel_index(np.argmax(pressure),pressure.shape);massedges=np.r_[0.,np.cumsum(w)];pm_it=int(np.argmax(pressure@w))
 peak_context={'time_s':float(time[peak_it]),'cell_index':int(peak_iz),'cell_mass_fraction':float(w[peak_iz]),'enclosed_mass_fraction_inner':float(massedges[peak_iz]),'enclosed_mass_fraction_outer':float(massedges[peak_iz+1]),'radius_inner_m':float(ri[peak_it,peak_iz]),'radius_outer_m':float(r[peak_it,peak_iz]),'temperature_K':float(T[peak_it,peak_iz]),'density_kg_m3':float(rho[peak_it,peak_iz]),'pressure_Pa':float(pressure[peak_it,peak_iz]),'mass_fraction_above_half_this_pressure_at_same_time':float(w[pressure[peak_it]>.5*pressure[peak_it,peak_iz]].sum())}
 pressure_partition_at_peak_mean={'time_s':float(time[pm_it]),'ion_fraction':float((pi@w)[pm_it]/(pressure@w)[pm_it]),'electron_positron_fraction':float((pe@w)[pm_it]/(pressure@w)[pm_it]),'LTE_radiation_fraction':float((pr@w)[pm_it]/(pressure@w)[pm_it]),'grey_center_to_surface_optical_depth':float(tau[pm_it]),'artificial_viscosity_over_thermodynamic_pressure':float((qvis@w)[pm_it]/(pressure@w)[pm_it])}
 inventory=[]
 for j,n in enumerate(names):
  mm=ynet[:,j]/cats;hist['inventory_nuclei_per_initial_catalyst__'+n]=mm;ii=int(np.argmax(mm));inventory.append(dict(isotope=n,initial_nuclei_per_initial_catalyst=float(mm[0]),final_nuclei_per_initial_catalyst=float(mm[-1]),net_change_per_initial_catalyst=float(mm[-1]-mm[0]),sampled_peak_nuclei_per_initial_catalyst=float(mm[ii]),sampled_peak_time_s=float(time[ii]),final_mass_fraction=float(A[j]*ynet[-1,j])))
 flow=None
 if do_flows and entry.get('reaction_flows',True):
  if d.get('rate_update') or d.get('mixture_inputs',{}).get('rate_update'):raise ValueError('This case replaced nuclear rates in process; set reaction_flows:false or use its exact arrays')
  wanted=d.get('mixture_inputs',{}).get('network_source_sha256')
  if wanted and wanted!=sha(ROOT/'pulse_model.py'):raise ValueError('Network source mismatch; disable reaction flows or supply original arrays')
  f,fr=sampled_flows(T,rho/1000,y,w);turn=trapezoid(f,time,axis=0)/cats;turnraw=trapezoid(fr,time,axis=0)/cats;delta=(ynet[-1]-ynet[0])/cats;closure=pm.STOICH@turnraw-delta
  top=np.argsort(turn)[::-1][:20];table=[dict(reaction=str(pm.RATES[j]),turnover_per_initial_catalyst=float(turn[j]),signed_turnover=float(turnraw[j]),Q_MeV=float(pm.Q[j]),weak=bool(pm.RATES[j].weak_type)) for j in top]
  unseen=set(range(NR));pairs=[]
  while unseen:
   j=min(unseen);unseen.remove(j);rev=[k for k in unseen if np.all(pm.STOICH[:,k]==-pm.STOICH[:,j])];value=turn[j]
   if rev:
    for k in rev:unseen.remove(k);value-=turn[k]
   pairs.append(dict(reaction_orientation=str(pm.RATES[j]),opposite_reactions=[str(pm.RATES[k]) for k in rev],net_turnover_in_orientation_per_initial_catalyst=float(value)))
  pairs.sort(key=lambda p:abs(p['net_turnover_in_orientation_per_initial_catalyst']),reverse=True)
  signedE=float(np.dot(turnraw,pm.Q)*cats*M/MU*MEV);qclosure=(signedE-nuclear[-1])/max(abs(nuclear[-1]),1)
  aggregates={}
  for n in ['n','he4','o16']:
   coeff=pm.STOICH[names.index(n)];aggregates[n]={'gross_produced_per_initial_catalyst':float(np.maximum(coeff,0)@turn),'gross_consumed_per_initial_catalyst':float(np.maximum(-coeff,0)@turn),'net_inventory_change_per_initial_catalyst':float(delta[names.index(n)])}
  hc=pm.STOICH[Z>=9].sum(axis=0);aggregates['F_and_heavier']={'gross_entry_per_initial_catalyst':float(np.maximum(hc,0)@turn),'gross_exit_per_initial_catalyst':float(np.maximum(-hc,0)@turn),'net_inventory_change_per_initial_catalyst':float(delta[Z>=9].sum())}
  flow={'scope':'Trapezoidal quadrature on saved samples; gross turnover clips tiny negative inventories, signed comparison retains solver polynomial. Reversible capture turnover is not net production. No cold-decay/recovery flow included.','max_stoichiometric_quadrature_residual_nuclei_per_initial_catalyst':float(max(abs(closure))),'stoichiometric_L1_quadrature_residual_per_initial_catalyst':float(sum(abs(closure))),'nuclear_energy_quadrature_relative_error':qclosure,'max_clipped_minus_signed_turnover':float(max(abs(turn-turnraw))),'numerical_scope':'Tiny flows below the inventory/step-error scale, especially neutron production and cancellation, are not resolved physical yield predictions. A positive-inventory quadrature can consume clipped numerical trace neutrons.','top_gross_reactions':table,'top_net_reaction_pairs':pairs[:20],'side_channels':aggregates}
  for j in top[:6]:hist['reaction_rate_per_s_per_initial_catalyst__'+str(pm.RATES[j])]=f[:,j]/cats
 residual=(K+U+weak+escaped-d['initial_energy_J']-nuclear)/d['initial_energy_J'];netdelta=np.diff(escaped);checks=[np.all(V>0),np.all(np.isfinite(pressure)),abs(M/d['mass_kg']-1)<1e-12,abs(K[-1]/d['final_kinetic_J']-1)<1e-8,abs(U[-1]/d['final_internal_J']-1)<1e-8,abs(escaped[-1]-d['radiation_escaped_J'])/d['initial_energy_J']<1e-10,abs(weak[-1]-d['weak_energy_escaped_upper_J'])/d['initial_energy_J']<1e-10,np.max(abs(y@A-1))<1e-7,float(np.min(y*A))>-1e-7,np.all(netdelta>-1e-8*d['initial_energy_J']),abs(time[-1]/d['elapsed_s']-1)<1e-12,np.max(abs(residual))<3e-4,np.all(qvis>=0),np.all(qvis[Vdot>=0]==0)]
 summary={'label':entry['label'],'family':entry['family'],'stage':d['stage'],'zones':nz,'samples':nt,'elapsed_s':float(time[-1]),'initial_H_per_catalyst':float(ynet[0,names.index('p')]/cats),'initial_energy_J':d['initial_energy_J'],'nuclear_gain':d['deposited_nuclear_over_initial_energy'],'opacity_m2_kg':kappa,'solver_saved_step_peak_temperature_K':d['peak_temperature_K'],'sampled_peak_temperature_K':float(T.max()),'sampled_temperature_peak_fraction_of_solver_step_peak':float(T.max()/d['peak_temperature_K']),'extrema':extrema,'peak_cell_pressure_context':peak_context,'pressure_partition_at_peak_massmean':pressure_partition_at_peak_mean,'initial_surface_radius_m':float(r[0,-1]),'final_surface_radius_m':float(r[-1,-1]),'final_surface_velocity_m_s':float(v[-1,-1]),'final_lumped_mass_rms_velocity_m_s':float(np.sqrt(2*K[-1]/M)),'sampled_max_negative_mass_fraction':float(max(0,-np.min(y*A))),'initial_pressure_massmean_Pa':float(hist['pressure_massmean_Pa'][0]),'initial_pressure_volumemean_Pa':float(hist['pressure_volumemean_Pa'][0]),'net_escaped_radiation_J':float(escaped[-1]),'net_escaped_radiation_over_prepared_energy':float(escaped[-1]/d['initial_energy_J']),'net_escaped_radiation_over_prepared_plus_deposited_nuclear':float(escaped[-1]/(d['initial_energy_J']+nuclear[-1]-weak[-1])),'sampled_luminosity_integral_over_actual_escape':float(trapezoid(L,time)/escaped[-1]) if escaped[-1]>0 else None,'maximum_sampled_energy_residual_fraction_initial':float(max(abs(residual))),'inventory':inventory,'reaction_flows':flow,'provenance':{'result':str(file),'result_sha256':sha(file),'trace':str(trace),'trace_sha256':sha(trace),'FastEOS_sha256':sha(ROOT/'fast_eos.py'),'electron_table_sha256':sha(ROOT/'electron-eos.npz'),'nuclear_module_sha256':sha(ROOT/'pulse_model.py'),'pynucastro_version':importlib.metadata.version('pynucastro'),'nuclear_array_sha256':{k:array_sha(getattr(pm,k)) for k in ['COEFF','RATE_INDEX','PREF','DEXP','REACT','STOICH','Q','QNU','AA','ZZ']}},'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':bool(all(checks))}}
 return summary,{k:np.asarray(v).tolist() for k,v in hist.items()}

def plots(out,summaries,histories):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22,'savefig.dpi':180})
 colors={'reference64':'#455a64','optimized64':'#0069b4','source64':'#d66014'}
 fig,ax=plt.subplots(4,3,figsize=(14,13))
 for i,stage in enumerate(['A','B1','B2','C13']):
  for s in summaries:
   if s['stage']!=stage or s['family'] not in colors:continue
   h=histories[s['label']];t=np.array(h['time_s'])*1e6;c=colors[s['family']];label=s['family'].replace('64','')+' H:'+format(s['initial_H_per_catalyst'],'.3g')
   ax[i,0].semilogy(t,h['pressure_massmean_Pa'],color=c,label=label)
   ax[i,1].plot(t,h['half_mass_radius_m'],color=c,label=label)
   ax[i,2].semilogy(t,np.maximum(np.array(h['surface_radiative_luminosity_W']),1),color=c,label=label)
  ax[i,0].set_ylabel(stage+'\nMass-mean pressure (Pa)');ax[i,1].set_ylabel('Half-mass radius (m)');ax[i,2].set_ylabel('Net surface luminosity (W)')
  for x in ax[i]:x.set_xlabel('Time from prescribed inflow (µs)')
  ax[i,0].legend(fontsize=8)
 fig.suptitle('Saved-sample pulse dynamics: original, optimized and source-profile cases\nConstant-grey LTE radiation; 64 radial zones; sampled extrema',fontsize=13);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(out/'pulse-dynamics.png');fig.savefig(out/'pulse-dynamics.pdf');plt.close(fig)
 fig,axs=plt.subplots(2,2,figsize=(12,9))
 for ax,stage in zip(axs.ravel(),['A','B1','B2','C13']):
  selected=next((s for s in summaries if s['stage']==stage and s['family']=='source64'),next(s for s in summaries if s['stage']==stage and s['family']=='optimized64'));h=histories[selected['label']];t=np.array(h['time_s'])*1e6;Ein=selected['initial_energy_J']
  for key,label,c in [('kinetic_J','Coherent kinetic','#0069b4'),('ion_internal_J','Ion internal','#007f5f'),('electron_positron_internal_J','Electron/positron internal','#8e44ad'),('LTE_radiation_energy_J','LTE radiation reservoir','#d66014'),('net_escaped_radiation_J','Cumulative escaped radiation','#222222')]:ax.plot(t,np.array(h[key])/Ein,label=label,color=c)
  ax.set_title(selected['label']);ax.set_xlabel('Time (µs)');ax.set_ylabel('Energy / prepared energy');ax.legend(fontsize=8)
 fig.suptitle('Energy partition in selected source-response family\nElectron term includes degeneracy and pair-creation energy; LTE radiation is already part of internal energy',fontsize=12);fig.tight_layout(rect=[0,0,1,.93]);fig.savefig(out/'pulse-energy-partitions.png');fig.savefig(out/'pulse-energy-partitions.pdf');plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path);ap.add_argument('--output',type=Path,default=P/'pulse-diagnostics');ap.add_argument('--no-reaction-flows',action='store_true');ap.add_argument('--no-plots',action='store_true');a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True);entries=json.loads(a.manifest.read_text()) if a.manifest else defaults();eos=FastEOS();summaries=[];histories={}
 for ent in entries:
  s,h=reduce(ent,eos,not a.no_reaction_flows);summaries.append(s);histories[s['label']]=h;print(json.dumps({'case':s['label'],'peak_massmean_P_Pa':s['extrema']['pressure_massmean_Pa']['value'],'escaped_over_input':s['net_escaped_radiation_over_prepared_energy'],'checks':s['checks']}),flush=True)
 result={'scope':__doc__,'code_sha256':sha(__file__),'definitions':['All pressure/density/radius extrema are over archived samples, not all solver steps.','Solver temperature peak is max over accepted solver states in the archived result, not a continuous-time supremum.','Pressure fields are thermodynamic EOS pressure; the separate artificial-viscosity fields describe the numerical shock stress/heating already present in the trajectory.','LTE radiation reservoir and pressure are already included in FastEOS internal energy/pressure.','Actual cumulative net escaped radiation is the model loss. No gross free-free emission is added as another loss.','The grey LTE model imposes equilibrium radiation; it does not resolve photon production, spectrum, absorption lines or non-LTE thermalization.','Reaction turnovers are sampled quadrature diagnostics; isotope inventories and integrated model energy are more accurate conserved ledgers.'],'cases':summaries,'checks':{'count':sum(s['checks']['count'] for s in summaries),'passed':sum(s['checks']['passed'] for s in summaries),'all_passed':all(s['checks']['all_passed'] for s in summaries)}}
 (a.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');(a.output/'histories.json.gz').write_bytes(gzip.compress(json.dumps(histories,separators=(',',':'),allow_nan=False).encode(),mtime=0));(a.output/'manifest.json').write_text(json.dumps(entries,indent=2)+'\n')
 with (a.output/'comparison.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=['label','stage','zones','initial_H_per_catalyst','nuclear_gain','initial_pressure_massmean_Pa','sampled_peak_pressure_massmean_Pa','sampled_peak_pressure_max_Pa','sampled_peak_density_massmean_kg_m3','sampled_min_halfmass_radius_m','sampled_peak_LTE_radiation_J','net_escaped_radiation_over_prepared_energy','net_escaped_radiation_over_prepared_plus_deposited_nuclear']);writer.writeheader()
  for s in summaries:
   row={k:s[k] for k in writer.fieldnames if k in s};row.update(sampled_peak_pressure_massmean_Pa=s['extrema']['pressure_massmean_Pa']['value'],sampled_peak_pressure_max_Pa=s['extrema']['pressure_max_Pa']['value'],sampled_peak_density_massmean_kg_m3=s['extrema']['density_massmean_kg_m3']['value'],sampled_min_halfmass_radius_m=s['extrema']['half_mass_radius_m']['value'],sampled_peak_LTE_radiation_J=s['extrema']['LTE_radiation_energy_J']['value']);writer.writerow(row)
 if not a.no_plots:plots(a.output,summaries,histories)
 print(json.dumps(result['checks']),flush=True)
if __name__=='__main__':main()
