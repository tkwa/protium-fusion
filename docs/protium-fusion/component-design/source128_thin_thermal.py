"""Finite-backing screen on the coherent source128 family, no new hydrodynamics.
All source heat retained for a conservative thermal-storage comparison.
Uses actual kinetic+advected-U arrival spectra and separate prompt proxies.
"""
from pathlib import Path
import ast,json,math,hashlib,sys
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erfc,zeta
from scipy.linalg import eigh_tridiagonal
P=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/component-design');ROOT=P.parent
sys.path[:0]=[str(ROOT),str(ROOT/'mixture-design')]
from fast_eos import FastEOS
import pulse_model as pm
import material_cycle as mc
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=P/'source-family-carbon-median128/integration.json';packet=json.loads(source.read_text());family=packet['families']['source128'];plant=family['capacity_plants']['2000'];cycle=family['cycle']
R=plant['radius_m'];AREA=plant['area_m2'];RHO=2000.;C=299792458.;CUT=20000.;GUARD=86400.;DT=.25
rates={'n13':.001159113256783153,'o14':.009817978825996571,'o15':.005681556796002065};daughters={'n13':'c13','o14':'n14','o15':'n15'}
# Import only pure kernel definitions from the separately validated historical
# solver. Its old pulse-loading/main program is never executed or relabeled.
old=P/'receiver_thickness.py';tree=ast.parse(old.read_text());wanted={'pulse_bins','front_step','history','tail_bound'}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name in wanted],type_ignores=[]),str(old),'exec'))
hashes={str(p):sha(p)for p in[source,old,P/'receiver-thickness-validation.json',ROOT/'fast_eos.py',ROOT/'electron-eos.npz',ROOT/'pulse_model.py',Path(mc.__file__)]};stages={};eos=FastEOS()
for cat,item in family['stage_inputs'].items():
 selected=item['selected'];cache=selected['receiver'];stage=selected['stage_data'];d=selected['record'];name=stage['label'];rp=Path(item['record_path']);tp=Path(cache['source_paths']['trace'])
 assert math.isclose(cache['radius_m'],1e9,rel_tol=1e-12)and cache['completed_requested_duration']and cache['usable_ballistic']
 assert item['record_sha256']==sha(rp)
 for p in[rp,tp]:hashes[str(p)]=sha(p)
 tr=np.load(tp);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);st=tr['state'].reshape(nz,nv,-1);names=d['state_variables'];nn=names.index('T/T0');m=tr['cell_mass_kg']
 ro=np.r_[0.,st[:,names.index('r/R0'),-1]*d['R0_m']];v=np.r_[0.,st[:,names.index('v/cs0'),-1]*float(tr['velocity_unit_cm_s'])/100]
 Y=np.maximum(st[:,:nn,-1],0.);bar=Y@pm.AA;rho=m/(4*np.pi/3*np.diff(ro**3));u=eos.state(rho/1000,st[:,names.index('T/T0'),-1]*d['T0_K'],(Y@pm.ZZ)/bar,Y.sum(axis=1)/bar).e*1e-4
 drdv=np.diff(ro)/np.diff(v);b=R-ro[:-1]+drdv*v[:-1];off=float(tr['time_s'][-1])-drdv-R/C
 assert np.all(off<0)and np.all(u>0)and np.all(np.diff(v)>0)
 decay={n:(mc.nuc(n).mass-mc.nuc(daughters[n]).mass-mc.DEFAULT_COLD_NU[n])*stage['hot_output'].get(n,0)*stage['attempts_per_pulse']*mc.MEV/AREA for n in rates}
 cold=(stage['Q_cold_MeV']-stage['nu_cold_MeV'])*stage['attempts_per_pulse']*mc.MEV/AREA
 assert abs(sum(decay.values())/cold-1)<1e-6
 K=float(np.sum(m/np.diff(v)*np.diff(v**3)/6));U=float(m@u);prompt=float(d['radiation_escaped_J'])/AREA
 s={'m':m,'v':v,'b':b,'off':off,'slope':m/np.diff(v),'u':u,'decay_fluence':decay,'fluence':(K+U)/AREA+prompt+sum(decay.values()),'K_fluence':K/AREA,'K_nodal_fluence':float(d['final_kinetic_J'])/AREA,'U_fluence':U/AREA,'prompt_fluence':prompt,'period':plant['pulse_interval_s'][cat],'prompt_peak_reference_K':cache['prompt']['radial_ray_arrival_surface_screen']['peak_surface_rise_K']*(cache['radius_m']/R)**2,'prompt_duration_s':cache['prompt']['radial_ray_arrival_surface_screen']['duration_s'],'minimum_slot':cache['timing_at_radius']['mass999']['occupied_s']*R/cache['radius_m'],'impulse':cache['ballistic']['momentum_fluence_Pa_s']*(cache['radius_m']/R)**2,'peak_pressure':cache['ballistic']['peak_momentum_flux_Pa']*(cache['radius_m']/R)**3,'record':str(rp),'record_sha256':sha(rp),'trace':str(tp),'trace_sha256':sha(tp)}
 assert abs(K/cache['ballistic']['kinetic_J_smooth_mass_conserving']-1)<1e-11
 assert abs(U/cache['ballistic']['endpoint_internal_J']-1)<1e-10
 stages[name]=s
big=[n for n in stages if n!='B1'];qmean=sum(s['fluence']/s['period']for s in stages.values());guard_used=GUARD*sum(1/stages[n]['period']for n in big);assert guard_used<1
edges=np.arange(round(CUT/DT)+1)*DT;times=edges[1:]
for s in stages.values():
 parts=pulse_bins(s,edges,True);s['q_material']=parts[0]+parts[1];s['q_decay']=parts[2];s['q_total']=sum(parts[:3]);s['actual_early_fraction']=float(sum(sum(z)for z in pulse_bins(s,edges,False)[:3])*DT/(s['fluence']-s['prompt_fluence']))

def curvature_step(t,L,k,cp,alpha_expansion):
 # Positive magnitude for front-face heating, uniform bulk mode cancels.
 t=np.asarray(t);n=np.arange(1,600,2);lam=k/(RHO*cp)*(n*math.pi/L)**2
 out=np.zeros_like(t)
 for i in range(0,len(t),4096):
  a=t[i:i+4096];out[i:i+4096]=48*alpha_expansion/(RHO*cp*L*L)*np.sum((-np.expm1(-lam[:,None]*a))/((n*math.pi)**2*lam)[:,None],axis=0)
 return out
rows=[]
for k,cp in[(50.,1800.),(25.,1800.),(100.,1800.)]:
 for sigma in[500.,750.,1000.]:
  L=sigma/RHO;step=front_step(edges,L,k,cp,'N');ker=np.diff(step);assert min(ker)>-1e-12
  curvstep=curvature_step(edges,L,k,cp,5e-6);cg=np.diff(curvstep);peaks={}
  for name,s in stages.items():
   material=fftconvolve(s['q_material'],ker)[:len(times)];decay=fftconvolve(s['q_decay'],ker)[:len(times)];slow=material+decay
   # Prompt has finished before first0.25s slow bin; all-source surface tail
   # overstates finite-depth heating. Keep its flash maximum separate.
   age=np.maximum(times-s['prompt_duration_s'],1e-10);ptail=s['prompt_fluence']/np.sqrt(math.pi*k*RHO*cp*age)
   slow+=ptail
   flash=s['prompt_peak_reference_K']*math.sqrt(100*1000/(k*cp))
   curv=fftconvolve(s['q_total'],cg)[:len(times)]
   # Impulsive prompt-curvature contribution, with conservative finite first bin.
   curv+=s['prompt_fluence']*np.r_[cg[0]/DT,cg[1:]/DT]
   j=int(np.argmax(slow));jc=int(np.argmax(curv));curvmax=float(curv[jc]); fullslow=slow+np.cumsum(s['q_total'])*DT/(sigma*cp); residual_front=slow-curv*L/(2*5e-6)
   peaks[name]={'material_only_nonuniform_peak_K':float(max(material)),'material_plus_decay_nonuniform_peak_K':float(max(material+decay)),'slow_plus_prompt_tail_nonuniform_peak_K':float(slow[j]),'slow_peak_time_s':float(times[j]),'isolated_slow_full_front_rise_K_including_its_bulk_mode':float(max(fullslow)),'full_front_peak_time_s':float(times[int(np.argmax(fullslow))]),'free_plate_front_residual_thermal_strain_K_peak':float(max(abs(residual_front))),'front_stress_scale_Pa_at_E10GPa_nu0p2_alpha5em6':float(max(abs(residual_front)))*1e10*5e-6/.8,'material_contribution_at_slow_peak_K':float(material[j]),'decay_contribution_at_slow_peak_K':float(decay[j]),'all_surface_prompt_flash_reference_K':flash,'maximum_K_plus_U_flux_W_m2':float(max(s['q_material'])),'free_plate_curvature_peak_1_m_at_alpha5em6':curvmax,'curvature_peak_time_s':float(times[jc]),'free_plate_one_axis_sag_at_span_m':{str(span):curvmax*span**2/8 for span in[1.,3.,5.,10.,100.]},'span_for1cm_one_axis_sag_m':math.sqrt(.08/curvmax),'span_for10cm_one_axis_sag_m':math.sqrt(.8/curvmax)}
  # Mean power uses actual7230s cadence; old-tail bound permits the actual
  # smaller sequential mass999 slot when dispatch jitters.
  smallP=stages['B1']['minimum_slot'];smallcut=min(CUT,smallP/2)
  oldbig=max(stages[n]['fluence']for n in big)*history(L,k,cp,GUARD,CUT,'N');oldsmall=stages['B1']['fluence']*history(L,k,cp,smallP,smallcut,'N')
  tail=(sum(tail_bound(stages[n],GUARD,CUT)for n in big)+tail_bound(stages['B1'],smallP,smallcut))*L/(3*k)
  bulk=sum(s['fluence']for s in stages.values())/(sigma*cp);offset=qmean*L/(6*k);background=bulk+oldbig+oldsmall+tail+offset
  slowbound=1000+background+max(peaks[n]['slow_plus_prompt_tail_nonuniform_peak_K']for n in big)+peaks['B1']['slow_plus_prompt_tail_nonuniform_peak_K']
  # This intentionally loose reference preserves full all-surface photon case;
  # source-specific finite-depth alternatives are attached separately.
  fullbound=1000+background+max(max(peaks[n]['slow_plus_prompt_tail_nonuniform_peak_K'],peaks[n]['all_surface_prompt_flash_reference_K'])for n in big)+max(peaks['B1']['slow_plus_prompt_tail_nonuniform_peak_K'],peaks['B1']['all_surface_prompt_flash_reference_K'])
  rows.append({'sigma_kg_m2':sigma,'thickness_m':L,'k_W_m_K':k,'cp_J_kg_K':cp,'diffusion_time_s':L*L*RHO*cp/k,'nonuniform_fundamental_decay_s':L*L*RHO*cp/k/math.pi**2,'bulk_mode_total_swing_upper_K':bulk,'constant_mean_pickup_W_m2':qmean,'constant_mean_back_extraction_front_offset_K':offset,'bulk_floor_offset_to_keep_pickup_at_least1000K':qmean*L/(2*k),'additional_bulk_swing_for_bounded1day_dispatch_lateness_K':qmean*GUARD/(sigma*cp),'old_big_nonuniform_background_K':oldbig,'old_B1_nonuniform_background_K':oldsmall,'late_tails_nonuniform_bound_K':tail,'common_background_upper_K':background,'slow_front_temperature_all_phase_bound_K':slowbound,'all_surface_photon_front_temperature_reference_bound_K':fullbound,'isolated_peaks':peaks,'scope':'Constant mean extraction at back, finite uniform energy mode; all retained source heat at front for slow/decay upper screen. Volume-average perstage storage modes shifted so sum minimum is1000K. Actual pickup temperature need not be exactly1000K.'})
# Independent finite-volume operator check of nonuniform front response.
fv=[]
for N in[128,512,1024]:
 dx=1/N;diag=np.full(N,-2/dx**2);diag[0]=diag[-1]=-1/dx**2;lam,V=eigh_tridiagonal(diag,np.full(N-1,1/dx**2));zero=abs(lam)<1e-7;lam[zero]=0;tt=np.geomspace(.001,10,100);gain=np.empty((N,len(tt)));gain[zero]=tt;gain[~zero]=-np.expm1(lam[~zero,None]*tt)/(-lam[~zero,None]);front=(V[0]*(V[0]/dx))@gain+dx/2-tt;analytic=front_step(tt,1,RHO,1,'N')*RHO;err=float(max(abs(front/analytic-1)));fv.append({'cells':N,'max_relative_error':err});assert N<1024 or err<.0003
refine=[]
for name,s in stages.items():
 pp=[]
 for dt in[.5,.25,.125]:
  ee=np.arange(round(CUT/dt)+1)*dt;qq=sum(pulse_bins(s,ee,True)[:3]);kk=np.diff(front_step(ee,.25,50,1800,'N'));pp.append(float(max(fftconvolve(qq,kk)[:len(qq)])))
 refine.append({'stage':name,'dt_s':[.5,.25,.125],'peak_K':pp,'finest_relative_to_base':pp[-1]/pp[1]-1});assert abs(pp[-1]/pp[1]-1)<.002
# Independently reconstruct electrical/global mass/radiator/cadence ledgers.
ledger=[]
for fn,fam in packet['families'].items():
 for sig,p in fam['capacity_plants'].items():
  cc=fam['cycle'];F=p['helium_per_s'];n=p['gross_electric_W']-sum(p[x]for x in['driver_W','ordinary_auxiliary_W','pump_W','replacement_W','support_W'])
  checks={'mass':sum(p['mass_breakdown_kg'].values())/p['mass_kg']-1,'electrical':n/p['net_W']-1,'global_nonnu':(p['nuclear_non_neutrino_W']-p['escaping_source_W']-p['net_W'])/p['cold_rejection_W']-1,'nuclear':p['nuclear_non_neutrino_W']/(cc['H_nonnu_MeV_per_He']*mc.MEV*F)-1,'radiator':(p['cold_rejection_W']+p['environment_absorbed_W'])/(p['radiator_area_m2']*p['radiator_emission_W_m2'])-1}
  for cat,sd in cc['stage_pulse_data'].items():checks['period_'+cat]=p['pulse_interval_s'][cat]*sd['pulses_per_He']*F-1
  assert max(abs(x)for x in checks.values())<1e-10
  ledger.append({'family':fn,'sigma':sig,'relative_residuals':checks})
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':hashes,'reference':{'source_family':'source128 carbon median, exact prescribed source profile remap','selected_cases':plant['selected_cases'],'radius_m':R,'area_m2':AREA,'cadence_from_receiver_sigma':2000,'cadence_He_s':plant['helium_per_s'],'B1_period_s':stages['B1']['period'],'B1_history_minimum_spacing_bound_s':stages['B1']['minimum_slot'],'large_guard_s':GUARD,'large_guard_capacity_fraction':guard_used,'all_source_mean_W_m2':qmean,'mean_smooth_vs_exact_nodal_relative':qmean/(plant['captured_pulse_W']/.99/AREA)-1},'pulse_fluences':{n:{k:s[k]for k in['K_fluence','K_nodal_fluence','U_fluence','prompt_fluence','decay_fluence','fluence','period','actual_early_fraction','record','record_sha256','trace','trace_sha256']}for n,s in stages.items()},'thermal_rows':rows,'finite_volume_kernel_validation':fv,'actual_profile_timestep_refinement':refine,'independent_integration_ledger_checks':ledger,'qualification':['All-source retention is a conservative hot-buffer screen, not inherited99% capture in thin receivers.','Constant k,cp,rho are requirements/scenarios; no high-temperature damage property law or coupled radiative-cooling/vapor model.','Photon grey flash and finite-depth spectral proxies remain separately identified; a grey luminosity does not determine an escaping spectrum.','The bulk-mode bound is valid for periodic stage streams with constant mean pickup. One-day separation between large pulses is an explicit dispatch requirement; arbitrary burst schedules can require more storage.','A positive back-face delay kernel maps the bulk-mode history to the pickup temperature with constant offset minus qmean*L/(2*k). Raising the bulk floor by that offset keeps pickup at least1000K without fixing its temperature.','A one-day bounded dispatch lateness relative to nominal stage release times adds at most qmean*86400/(Sigma*cp) to the summed bulk storage range; this is an explicit schedule acceptance sensitivity.','Free plate curvature uses a selected constant thermal-expansion coefficient5e-6/K. This is a thermal-strain/segmentation screen, not a fatigue/joint or irradiated-material qualification.','No fixed1000K back-face or instantaneous peak heat-engine capacity is assumed.']}
(P/'source128-thin-thermal.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'reference':out['reference'],'summary':[{'sigma':r['sigma_kg_m2'],'k':r['k_W_m_K'],'bulk_K':r['bulk_mode_total_swing_upper_K'],'slow_bound_K':r['slow_front_temperature_all_phase_bound_K'],'all_surface_photon_bound_K':r['all_surface_photon_front_temperature_reference_bound_K'],'C13':r['isolated_peaks']['C13']}for r in rows],'FV':fv,'refine':refine,'ledger_max_relative':max(abs(x)for r in ledger for x in r['relative_residuals'].values())},indent=2))
