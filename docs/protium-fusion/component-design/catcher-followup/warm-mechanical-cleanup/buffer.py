"""Exact-dispatch, finite-conductance two-reservoir hot-gas buffer screen.
The selected operating initial condition is a precharged steady mean inventory.
All sources are retained; arrivals after the finite record remain inventory.
"""
from pathlib import Path
import json, math, hashlib, importlib.util
import numpy as np
from scipy.signal import lfilter
P=Path(__file__).resolve().parent
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website')
BASE=ROOT/'docs/protium-fusion/component-design/catcher-followup'
G=BASE/'gas'; W=BASE/'warm-mechanical-cleanup'
read=lambda x:json.loads(Path(x).read_text());sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('gas_chamber',G/'chamber.py');ch=importlib.util.module_from_spec(spec);spec.loader.exec_module(ch)
d=read(G/'inputs.json');pi=read(G/'pumping-interfaces.json');warm=read(W/'results.json')
KB=1.380649e-23;SB=5.670374419e-8;T=1000.;TC=300.;A=4*math.pi*d['R_m']**2
m=np.array([2*d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']]);vb=np.sqrt(8*KB*T/(math.pi*m));qmean=np.array(pi['mean_throughput_Pa_m3_s_per_wall_m2']);h=pi['geometry']['plenum_void_equivalent_m']+pi['geometry']['slit_void_equivalent_m'];Sch=np.array(pi['geometry']['species_channel_speed_m3_s_per_wall_m2']);Sp=.001;Se=1/(1/Sch+1/Sp)
# One physical molecular-flow outlet: speed scales as mean thermal speed.
# beta chosen for 1Pa summed hot partial pressure at the continuous mean flow.
beta=float(np.sum(qmean/(vb/4)));Sc=beta*vb/4;pcold=1e-4*qmean/sum(qmean);floor=pcold*math.sqrt(T/TC)
checks=[]
def check(name,ok):
 assert bool(ok),name
 checks.append(name)

def sources(dt,points):
 ts,HH,_=ch.emission(d,'p',dt,points);_,HE,_=ch.emission(d,'he4',dt,points)
 # chamber.emission assigns arrival at ts[i]. Apply the associated uniform
 # bin source on [ts[i],ts[i]+dt]; this adds 0..dt causal smoothing, no advance.
 q=np.vstack([HH,HE])*KB*T/(A*m[:,None])/dt
 total_dispatched=np.array([sum(d['stages'][e['cat']]['cold_'+ss+'_atomic_kg'] for e in d['events']) for ss in ['H','He']])*KB*T/(A*m)
 return ts,q,{'dt_s':dt,'arrival_quantiles_per_stage_species':points,'duration_s':len(ts)*dt,'emitted_PV_Pa_m_by_species':(q.sum(axis=1)*dt).tolist(),'dispatched_total_PV_Pa_m_by_species':total_dispatched.tolist(),'not_yet_arrived_PV_Pa_m_by_species':(total_dispatched-q.sum(axis=1)*dt).tolist(),'time_bin_scope':'ceil-causal arrival followed by uniform source over one dt bin; no source is shifted earlier. This is adequate for multi-day buffer pressure; exact local peak is reported separately.'}

def response(q,dt,V,K,initial='mean'):
 n=q.shape[1];Tend=n*dt;buf=[];src=[];stats=[]
 for j in range(2):
  M=np.array([[-Se[j]/h,Se[j]/(h*K)],[Se[j]/V,-(Se[j]/K+Sc[j])/V]])
  lam,U=np.linalg.eig(M);Ui=np.linalg.inv(U);B=np.array([1/h,0.]);c=(np.expm1(lam*dt)/lam)*(Ui@B);rr=np.exp(lam*dt)
  # subtract the exact mean-flow equilibrium, then filter each stable mode.
  eq=np.array([qmean[j]/Se[j]+(qmean[j]/Sc[j]+floor[j])/K,qmean[j]/Sc[j]+floor[j]])
  z0=np.zeros(2) if initial=='mean' else Ui@(-eq)
  mode=[]
  for k in range(2):
   z,_=lfilter([c[k]],[1.,-rr[k]],q[j]-qmean[j],zi=[rr[k]*z0[k]])
   mode.append(z)
  x=eq[:,None]+U@np.array(mode)
  x0=eq if initial=='mean' else np.zeros(2)
  # Independently integrate the exact modal trajectories over every bin.
  bmodal=Ui@B
  intm=(np.expm1(lam*dt)/lam)*(np.sum(mode,axis=1)-np.array(mode)[:,-1]+z0)+(bmodal/lam)*(np.expm1(lam*dt)/lam-dt)*sum(q[j]-qmean[j])
  integrated=eq*Tend+U@intm
  pintegral=float(integrated[1])
  balance_pintegral=(sum(q[j])*dt-h*(x[0,-1]-x0[0])-V*(x[1,-1]-x0[1]))/Sc[j]+floor[j]*Tend
  check(f'independent_integral_V{V}_K{K}_{initial}_{j}',abs(pintegral-balance_pintegral)<1e-7*max(1,abs(pintegral)))
  sintegral=(sum(q[j])*dt-h*(x[0,-1]-x0[0]))/Se[j]+pintegral/K
  output=Sc[j]*(pintegral-floor[j]*Tend)
  resid=sum(q[j])*dt-output-h*(x[0,-1]-x0[0])-V*(x[1,-1]-x0[1])
  check(f'mass_V{V}_K{K}_{initial}_{j}',abs(resid)<1e-11*max(1,sum(q[j])*dt))
  check(f'stable_V{V}_K{K}_{initial}_{j}',max(lam)<0)
  buf.append(x[1]);src.append(x[0]);stats.append({'species':['H2','He'][j],'buffer_pressure_min_Pa':float(min(x[1])),'buffer_pressure_max_Pa':float(max(x[1])),'buffer_pressure_exact_mean_Pa':float(pintegral/Tend),'source_pressure_discretized_max_Pa':float(max(x[0])),'source_pressure_exact_mean_Pa':float(sintegral/Tend),'source_relaxation_s':h/Se[j],'buffer_slow_eigenmode_s':float(-1/max(lam)),'initial_buffer_pressure_Pa':float(x0[1]),'initial_source_pressure_Pa':float(x0[0]),'final_buffer_pressure_Pa':float(x[1,-1]),'source_PV_Pa_m':float(sum(q[j])*dt),'cold_out_PV_Pa_m':float(output),'mass_balance_residual_Pa_m':float(resid)})
 buf=np.array(buf);src=np.array(src);total=buf.sum(axis=0);out=(Sc[:,None]*(buf-floor[:,None])).sum(axis=0)
 i=int(np.argmax(total));row={'buffer_volume_per_wall_area_m':V,'booster_zero_flow_K0_each_species':K,'initial_condition':initial,'total_buffer_pressure_min_Pa':float(min(total)),'total_buffer_pressure_max_Pa':float(max(total)),'total_buffer_pressure_exact_mean_Pa':sum(s['buffer_pressure_exact_mean_Pa'] for s in stats),'peak_time_s':float((i+1)*dt),'extrema_scope':'Sampled exact bin endpoints; within-bin extrema are not searched. The independent25/50/100s refinement quantifies the small buffer discretization change.','total_source_pressure_discretized_max_Pa':float(max(src.sum(axis=0))),'cold_outlet_PV_flow_min_Pa_m_s':float(min(out)),'cold_outlet_PV_flow_max_Pa_m_s':float(max(out)),'cold_outlet_peak_over_continuous_mean':float(max(out)/sum(qmean)),'species':stats}
 return row,total,buf,src

outrows=[];convergence=[];curves={}
for dt,points in [(100.,256),(50.,1024),(25.,2048)]:
 ts,q,meta=sources(dt,points)
 Vs=[.5] if dt!=50 else [.1,.5,1.]
 for V in Vs:
  r,tot,bf,sr=response(q,dt,V,1e6);r['source_discretization']=meta;outrows.append(r)
  if V==.5: convergence.append({'dt_s':dt,'points':points,'max_Pa':r['total_buffer_pressure_max_Pa'],'min_Pa':r['total_buffer_pressure_min_Pa'],'mean_Pa':r['total_buffer_pressure_exact_mean_Pa']})
  if dt==50:
   keep=np.unique(np.r_[np.arange(0,len(ts),20),np.argmax(tot),np.argmin(tot)]);curves[str(V)]={'time_s':((keep+1)*dt).tolist(),'total_pressure_Pa':tot[keep].tolist()}
 if dt==50:
  for K in [1e5]:
   r,*_=response(q,dt,.5,K);r['source_discretization']=meta;outrows.append(r)
  r,*_=response(q,dt,.5,1e6,'empty');r['source_discretization']=meta;outrows.append(r)
  q50=q;dt50=dt
check('pressure_resolution_max_under0p2percent',abs(convergence[-1]['max_Pa']/convergence[-2]['max_Pa']-1)<.002)
ref=next(r for r in outrows if r['buffer_volume_per_wall_area_m']==.5 and r['booster_zero_flow_K0_each_species']==1e6 and r['initial_condition']=='mean' and r['source_discretization']['dt_s']==50)
# Exact-cell source-peak bound from frozen gas packet plus the maximum
# buffer backpressure. The sampled bin peaks are not used for pump sizing.
local=next(x for x in pi['continuous_extraction_scenarios'] if x['downstream_pump_speed_m3_s_per_wall_m2_each_species']==.001)
peaks={s:p+ref['total_buffer_pressure_max_Pa']/1e6 for s,p in local['peak_pressure_by_stage_Pa'].items()}
# Physical outlet matches the older 1mm-radius,0.1m-long dogleg conduit.
co=warm['concrete_downstream_conduit'];beta_available=co['net_aperture_beta'];valve=beta/beta_available
# Actual pressure minimum determines the finite reverse-flow correction.
ratios=[]
for j,s in enumerate(ref['species']):
 ratios.append({'species':s['species'],'operating_hot_partial_Pmin_Pa':s['buffer_pressure_min_Pa'],'cold_partial_Pa':float(pcold[j]),'thermal_reverse_flux_fraction_at_operating_min':float(floor[j]/s['buffer_pressure_min_Pa']),'gross_outlet_speed_m3_s_per_m2':float(Sc[j]),'net_outlet_speed_at_operating_min_m_s':float(Sc[j]*(1-floor[j]/s['buffer_pressure_min_Pa']))})
# Small booster has measured-aperture-scale capture preserved at1000K.
alpha=np.array([warm['reference']['pump_capture_H2_requirement'],warm['reference']['pump_capture_He_comparison']]);fp_req=Sp/(alpha*vb/4);fp=float(max(fp_req));mA=warm['reference']['commercial_whole_pump_mass_per_nominal_aperture_kg_m2'];idle=warm['idle_loss_primary_comparison']['primary_measured_input_W_per_disk_m2']
booster={'specified_inlet_speed_each_species_m3_s_m2':Sp,'required_pump_face_fraction_from_benchmark_capture':fp_req.tolist(),'shared_face_selected_fraction':fp,'package_mass_at_benchmark_areal_density_kg_m2':fp*mA,'idle_W_m2_at_measured_small_pump_areal_input':fp*idle,'zero_flow_K0_reference':1e6,'zero_flow_K0_sensitivity':1e5,'source_peak_upper_Pa_by_stage':peaks,'scope':'Preserved hot capture, compression and package areal density are conditional engineering requirements. Actual small-stage counts/rotor dimensions are not derived. Source pressure rises during pulses; no instantaneous transfer or constant0.001Pa inlet is imposed.'}
# Nitrogen is a release sensitivity, not a selected chemical yield. Use fully
# decayed output, form N2, and display exact event PV/jump and steady flow.
Nrows=[];Nmean=0;N14rate=0.;N15rate=0.
for cat,s in d['stages'].items():
 N=(s['cold_output_counts'].get('n14',0)+s['cold_output_counts'].get('n15',0))*s['attempts_per_pulse']/2
 pv=N*KB*T/A;period=d['plant']['pulse_interval_s'][cat];Nmean+=pv/period
 N14rate+=s['cold_output_counts'].get('n14',0)*s['attempts_per_pulse']/period
 N15rate+=s['cold_output_counts'].get('n15',0)*s['attempts_per_pulse']/period
 Nrows.append({'stage':s['stage'],'N2_PV_Pa_m_if_all_nitrogen_mobile':pv,'mean_PV_Pa_m_s':pv/period,'instantaneous_pressure_jump_at_V0p5_Pa':pv/.5,'scope':'All stable N atoms hypothetically form N2 and enter this hot buffer. This reports an event-inventory bound; release timing/speciation are not inferred.'})
mN=2*(14*N14rate+15*N15rate)/(N14rate+N15rate)*1.66053906660e-27;vN=math.sqrt(8*KB*T/(math.pi*mN));SN=beta*vN/4
# Bursty N2 sensitivity: every stage's eventual nitrogen is released at
# the earliest resolved material arrival. This is a declared time-profile
# sensitivity, not a measured nitrogen release history or a strict maximum.
nq=np.zeros_like(q50[0]);pv_bycat={cat:row['N2_PV_Pa_m_if_all_nitrogen_mobile'] for cat,row in zip(d['stages'],Nrows)}
for e in d['events']:
 st=d['stages'][e['cat']];vmax=max(st['velocity_bounds_m_s']);rmax=max(st['radius_bounds_m']);tfirst=st['endpoint_time_s']+(d['R_m']-rmax)/vmax
 ix=int(math.ceil((e['focus_s']+tfirst)/dt50))
 if 0<=ix<len(nq):nq[ix]+=pv_bycat[e['cat']]/dt50
rN=math.exp(-SN/.5*dt50);Ncurve=lfilter([(-math.expm1(-SN/.5*dt50))/SN],[1,-rN],nq-Nmean)+Nmean/SN
# Reference H/He exactly on the same source grid.
_,base_curve,_,_=response(q50,dt50,.5,1e6)
Nscheduled=[{'volatile_N_fraction':f,'total_buffer_pressure_min_Pa':float(min(base_curve+f*Ncurve)),'total_buffer_pressure_max_Pa':float(max(base_curve+f*Ncurve)),'N2_outflow_peak_PV_Pa_m_s':float(max(SN*f*Ncurve)),'scope':'Finite0.5m buffer and physicalN2 outlet. N feed is intentionally impulsive at earliest resolved material arrival; its finite booster delay is unresolved and can smooth this sensitivity.'} for f in [0,.1,1.]]
Ncases=[{'volatile_N_fraction':f,'added_continuous_mean_N2_pressure_Pa_same_passive_outlet':f*Nmean/SN,'added_mean_N2_PV_Pa_m_s':f*Nmean,'largest_single_event_pressure_jump_Pa_V0p5':f*max(x['N2_PV_Pa_m_if_all_nitrogen_mobile'] for x in Nrows)/.5} for f in [0,.1,1.]]
# An explicit small-cassette hot storage box:0.5m depth,1m face span,0.5mm
# dense carbon-composite skins,posts/ties. Not a freestanding global vessel.
rho=1800.;t=.0005;pitch=1.;V=.5;allow=100e6;E=10e9;pdesign=10.;panelstress=.75*pdesign*pitch**2/t**2;panelbow=5*pdesign*pitch**4/(32*E*t**3)
massskins=2*rho*t;massside=2*V*rho*t/pitch;ties=rho*pdesign*V/allow
vessel={'void_depth_m':V,'cassette_span_m':pitch,'two_skin_thickness_each_m':t,'sidewall_thickness_m':t,'selected_density_kg_m3':rho,'design_gauge_pressure_Pa':pdesign,'skin_and_sidewall_areal_mass_kg_m2':massskins+massside,'ideal_tensile_tie_mass_floor_kg_m2':ties,'stress_allowance_Pa_requirement':allow,'simply_supported_one_way_strip_stress_Pa':panelstress,'linear_one_way_strip_bow_m':panelbow,'E_Pa_scenario':E,'mechanical_scope':'Flat skins require closer ribs, corrugation or tension to limit bow; this is an explicit scale screen. Both reservoir walls remain near1000K, with the station cold insulation outside. No isolated global pressure-shell stress is assigned. Seals, valves, vibration, access and actual refractory laminate qualification require the named service allowance.'}
# Choose0.1m stiffening pitch to make a directly usable flat-skin screen.
vessel['ribs_0p1m_pitch_strip_stress_Pa']=panelstress*.1**2;vessel['ribs_0p1m_pitch_strip_bow_m']=panelbow*.1**4
vessel['two_faces_0p1m_grid_ribs_10mm_deep_0p5mm_thick_mass_kg_m2']=2*2/.1*.01*t*rho
vessel['explicit_ties_cross_section_m2']=1e-6
vessel['explicit_tie_pitch_m']=.1
vessel['explicit_tie_mass_kg_m2']=rho*V*1e-6/.1**2
vessel['tie_stress_at_design_pressure_Pa']=pdesign*.1**2/1e-6
vessel['explicit_hot_box_mass_with_ribs_kg_m2']=vessel['skin_and_sidewall_areal_mass_kg_m2']+vessel['two_faces_0p1m_grid_ribs_10mm_deep_0p5mm_thick_mass_kg_m2']+vessel['explicit_tie_mass_kg_m2']
vessel['tied_grid_scope']='Tension ties connect both hot skins on the0.1m grid. Face ribs span between ties, not across an unsupported1m cassette. Internal positive gauge pressure loads ties in tension. Off-design handling and thermal mismatch remain separate.'
# Separate global-cleanup foreline is boosted into the same reservoir only
# after compression, so primary pulses never set its low pressure.
global_q=[]
for j,ss in enumerate(['H2','He']):
 row=next(z for z in warm['species'] if z['species']==ss)
 global_q.append(row['g_if_K1e9_for_this_species']*vb[j]/4*1e-10)
global_q=np.array(global_q)
global_booster_S=global_q/(1e-3-ref['total_buffer_pressure_max_Pa']/1e6)
global_extra={'simultaneous_reference_partial_Pin_each_Pa':1e-10,'realized_g_H2_He':[next(z for z in warm['species'] if z['species']==ss)['g_if_K1e9_for_this_species'] for ss in ['H2','He']],'gas_PV_Pa_m_s':global_q.tolist(),'ratio_to_primary_mean_PV':float(sum(global_q)/sum(qmean)),'separate_global_foreline_pressure_each_Pa':.001,'required_booster_speed_each_m_s_at_buffer_upper_bound':global_booster_S.tolist(),'added_buffer_mean_partial_pressure_Pa':(global_q/Sc).tolist(),'scope':'These are simultaneous species flows through the fixed H2-sized global geometry at the specified partial pressures. They enter through a distinct small booster; no direct pressure connection joins the primary pulse plenum to the global-cleanup foreline.'}
ref['initial_hot_buffer_gas_mass_kg_m2']=float(sum(np.array([v['initial_buffer_pressure_Pa'] for v in ref['species']])*.5*m/(KB*T)))
ref['initial_hot_buffer_gas_mass_kg']=ref['initial_hot_buffer_gas_mass_kg_m2']*A
seal=[{'fraction_of_primary_mean_flow_to_guard':f,'allowed_total_leak_conductance_per_wall_area_m_s_at_reference_Pmax':f*sum(qmean)/ref['total_buffer_pressure_max_Pa'],'equivalent_H2_effusive_open_area_fraction':f*sum(qmean)/ref['total_buffer_pressure_max_Pa']/(vb[0]/4)} for f in [1e-4,.001,.01]]
booster['idle_W_m2_at_modern_next85H20W_per_DN63_area']=fp*20/(math.pi*.063**2/4)
compression={'steady_primary_1000K_to1Pa_from_mean_source_50percent_work_W_m2':float(2*sum(qmean*np.log(np.maximum(1.,(qmean/Sc)/(qmean/Se))))),'gas_work_scope':'Mean-state ideal-gas free-energy scale only. Actual pulsed pump work, windage and idle require the separate apparatus budget; no compression work is credited when the source exceeds the buffer pressure.'}
month=365.25*86400/12
regen=[]
for retained in [.01,.1,.5,1.]:
 pv=retained*sum(qmean)*month
 regen.append({'retained_fraction_same_H_He':retained,'hold_time_s':month,'retained_cohort_PV_Pa_m':pv,'instant_release_pressure_jump_Pa_V0p5':pv/.5,'minimum_release_duration_s_if_total_outlet_capped2mean':retained/(1+retained)*month,'He_cohort_pressure_jump_Pa_V0p5':retained*qmean[1]*month/.5})
out={'title':'Exact-dispatch hot foreline buffer','family':d['family'],'gas_order':['H2','He'],'source_sha256':{str(x.relative_to(ROOT)):sha(x) for x in [G/'inputs.json',G/'pumping-interfaces.json',G/'chamber.py',W/'results.json']},'code_sha256':sha(__file__),'parameters':{'hot_T_K':T,'cold_intercept_T_K':TC,'source_plenum_void_depth_m':h,'channel_speeds_m_s':Sch.tolist(),'booster_inlet_speed_m_s':Sp,'source_series_speeds_m_s':Se.tolist(),'continuous_mean_PV_Pa_m_s':qmean.tolist(),'outlet_effective_area_beta':beta,'passive_species_outlet_speeds_m_s':Sc.tolist(),'cold_total_pressure_Pa':sum(pcold),'cold_partial_pressure_Pa':pcold.tolist(),'buffer_initialization':'Mean-flow equilibrium before the finite actual-dispatch record; actual precharge stock is explicitly recorded. Empty startup is separate.'},'cases':outrows,'convergence':convergence,'reference':ref,'cold_conduit':{'benchmark_pipe':co,'effective_transmission_setting_fraction':valve,'beta_required':beta,'beta_available_full_open':beta_available,'actual_operating_min_backstream':ratios,'all_admitted_hot_radiation_W_m2_upper':co['all_admitted_hot_radiation_to_cold_W_m2_upper'],'solid_conduction_W_m2':co['pipe_wall_conduction_W_m2'],'scope':'Same geometrical mouth/pipe as existing isolated foreline example, with valve/restriction reducing molecular transmission. Full mouth radiation is still charged; flow throttling alone receives no optical credit. Output is intercepted at300K before retained He cryosorption. Cold capture capacity and refrigeration remain separately charged.'},'primary_warm_booster':booster,'hot_storage_box':vessel,'hermetic_envelope_requirements':{'leak_conductance_sensitivities':seal,'scope':'Dense liner/seals or a separately collected guard envelope are required. Raw porous graphite is not assigned gas-tightness. The fraction is diverted internally to the guard, not permanent material loss. Connection to the central cavity would require its own stricter gas-column budget; this packet does not assign that routing.'},'compression_scale':compression,'separate_global_foreline_booster':global_extra,'regeneration':{'reference_hold_s':month,'total_cold_outlet_cap_over_continuous_mean_scenario':2.,'cases':regen,'all_mobile_N2_one_month_extra_pressure_jump_Pa_V0p5':Nmean*month/.5,'scope':'Retained cohorts replace their prompt source, rather than being double-counted. The release-duration formula assumes backgroundprompt=(1-f)mean and totalregulatedoutlet<=2mean. Localvalves must limit discharge, or multipletiles mustshare a capacity-qualifiedheader/sorberbank before scheduling can average their loads. The storedcohort maystay in traps/sealedholdingcells until admitted. The hotbuffer pressure vessel and seals mustaccept any chosenburst inventory; no automaticcoldcapacity credit is assigned.'},'nitrogen_sensitivity':{'N2_mean_mass_kg_baryonic_isotope_proxy':mN,'N15_fraction_by_continuous_nitrogen_atoms':N15rate/(N14rate+N15rate),'N2_speed_same_outlet_m_s':SN,'stages':Nrows,'cases':Ncases,'scheduled_impulsive_release_sensitivity':Nscheduled,'outlet_beta_for1Pa_total_mean_if_allN2':beta*(1+Nmean/SN),'outlet_beta_full_open_margin_for1Pa_allN2':beta_available/(beta*(1+Nmean/SN)),'chemical_separator':'Retain a separate staged H/N/O chemical recovery branch ahead of the helium sorber. Nitrogen-bearing molecules and oxygen-bearing intermediates require removal before He cryosorption. No present separator efficiency is assigned; the root account carries its positive-work and apparatus allowances.','scope':'No mobile fraction or molecular form is selected; these bound additional inventory and mean throughput. Gas speed uses the continuous actualN14/N15 inventory-weighted mean molecular mass with baryonic14/15mass approximation. Stage-dependent isotopologue velocity and realN/O molecular chemistry remain separate. Actual N time-history is not substituted for H/He.'},'global_booster_and_buffer_scope':'The tabulated primary H/He dynamic curves exclude the separate global-cleanup stream; its specified simultaneous-flow increment is reported and adds about0.1percent to inventory/outlet demand. The physical conduit retains substantial unthrottled margin.','checks_passed':len(checks),'checks':checks,'limits':['No qualified1000K turbo/drag pump, seals or hot rotor are asserted.','Source relaxation is explicitly finite. The narrow gas arrival peaks use the frozen cellwise analytic upper screen, rather than the buffer time bins.','Gas is instantaneously thermalized/released at physical wall arrival for this maximum-throughput response; measured trapping and chemistry can delay it.','An empty buffer initially has vanishing driving pressure to the cold outlet; it must fill or be precharged before the quoted operating minimum applies.','The tiny cold gas conduit does not erase cold cryosorber view, cycling, working capacity or refrigerator loads.','Local pressure/temperature and mean free path may invalidate molecular flow in a chosen booster geometry; its high-pressure stages require a physical drag/channel design.']}
(P/'buffer.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');(P/'buffer-curves.json').write_text(json.dumps(curves,separators=(',',':'))+'\n')
print(json.dumps({'reference':ref,'convergence':convergence,'booster':booster,'conduit':out['cold_conduit'],'N':Ncases,'box':vessel,'checks':len(checks)},indent=2))
