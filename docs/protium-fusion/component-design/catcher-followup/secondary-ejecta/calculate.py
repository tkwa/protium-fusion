"""Bounded carbon-ejecta source and ballistic recapture screen.
The yield fits and RES expression are extrapolation diagnostics, not a reactor yield model.
"""
from pathlib import Path
import csv,json,math,hashlib,sys
import numpy as np
P=Path(__file__).resolve().parent; B=P.parent
sys.path.insert(0,str(B/'stopping'));import spectra as sp
MU=1.66053906660e-27;EV=1.602176634e-19;KB=1.380649e-23;MC=12*MU;R=1.1e9;AREA=4*math.pi*R*R;YEAR=365.25*86400
sha=sp.sha
sources={}
def read(p):
 p=Path(p);sources[str(p)]=sha(p);return json.loads(p.read_text())
surf=read(B/'stopping/surface-spectra.json');temp=read(B/'geometry/temperature-histories.json');gas=read(B/'gas/source-leakage-followup/results.json')
sources[str(B/'stopping/spectra.py')]=sha(B/'stopping/spectra.py')
rs,pm,st,decay,weak,hashes=sp.setup();sources.update(hashes)
checks=[]
def check(name,condition,**data):
 checks.append(dict(name=name,passed=bool(condition),**data));assert condition,(name,data)
def vb(T,m=MC):return math.sqrt(8*KB*T/(math.pi*m))
def group(z):return {1:'H',2:'He',6:'C',7:'N',8:'O'}.get(z,'other')
groups=['H','He','C','N','O']
# IPP9/82 pp42,46 calculated individual fits; N is explicitly interpolated sensitivity.
fits={1:(27.3,.0271),2:(25.4,.169),6:(44.,.620),7:(72.,.805),8:(100.,.990)}
fit_exp={1:(27.3,.0271),2:(30.2,.387),6:(35.,1.20),7:(67.5,1.095),8:(100.,.990)}
def yield_fit(E,A,Z,exp=False):
 Eth,Q=(fit_exp if exp else fits)[Z]
 ETF=30.74*(A+12)/12*Z*6*math.sqrt(Z**(2/3)+6**(2/3))
 eps=np.maximum(E/ETF,1e-100)
 sn=.5*np.log1p(1.2288*eps)/(eps+.1728*np.sqrt(eps)+.008*eps**.1504)
 ratio=np.minimum(Eth/np.maximum(E,1e-100),1.)
 return Q*sn*(1-ratio**(2/3))*(1-ratio)**2
stages={};rate={g:0. for g in groups};model_rates={'physical_fit_calculated':0.,'physical_fit_experimental_CO_proxy':0.,'RES_temperature_flux_extrapolation':0.,'RES_only_T_ge1200_extrapolation':0.};crosschecks=[]
for stage,row in surf['stages'].items():
 x=sp.packets(Path(row['record_path']),Path(row['trace_path']),pm,decay,order=32)
 fine=sp.packets(Path(row['record_path']),Path(row['trace_path']),pm,decay,order=64)
 sources[row['record_path']]=sha(row['record_path']);sources[row['trace_path']]=sha(row['trace_path'])
 profiles=[t for t in temp['profiles'] if t['stage']==stage];tp=max(profiles,key=lambda z:z['maximum_150um_skin_front_K'])
 arr=np.genfromtxt(tp['csv'],delimiter=',',names=True);sources[tp['csv']]=sha(tp['csv'])
 rel=x['arrival']-R/pm.C
 T=np.interp(rel,arr['relative_to_prompt_arrival_s'],arr['front_with150um_skin_k4_K'])
 allcounts=x['n'][:,pm.ZZ>0].sum(axis=1);flux=allcounts/x['dvweight']*x['arrival_jacobian'];gammacm=flux/1e4
 out={'period_s':row['period_s'],'temperature_profile':tp['csv'],'groups':{},'models':{},'total_ion_fluence_m2':float(allcounts.sum()),'RES_ion_fluence_fraction_within_calibration_flux':float(allcounts[(gammacm>=1e13)&(gammacm<=1e17)].sum()/allcounts.sum())}
 contrib={key:np.zeros(len(x['eu'])) for key in model_rates};finecon=np.zeros(len(fine['eu']))
 for G in groups:
  js=[j for j,z in enumerate(pm.ZZ) if group(z)==G];n=x['n'][:,js].sum(axis=1);N=float(n.sum());period=row['period_s'];rate[G]+=N/period
  E=np.concatenate([x['eu']*pm.AA[j]*1e6 for j in js]);weights=np.concatenate([x['n'][:,j] for j in js]);lo=float(weights[E>1e4].sum()/N) if N else 0;mev=float(weights[E>1e6].sum()/N) if N else 0
  out['groups'][G]={'incident_nuclei_m2':N,'carbon_mass_kg_m2_per_Y1':N*MC,'number_rate_m2_s':N/period,'incident_energy_MeV_quantiles':sp.quantile(E/1e6,weights) if N else {},'fraction_above_10keV':lo,'fraction_above_1MeV':mev}
  for j in js:
   E=x['eu']*pm.AA[j]*1e6;counts=x['n'][:,j];Z=int(pm.ZZ[j]);A=float(pm.AA[j]);yc=yield_fit(E,A,Z);ye=yield_fit(E,A,Z,True)
   contrib['physical_fit_calculated']+=counts*yc
   contrib['physical_fit_experimental_CO_proxy']+=counts*ye
   # IAEA1991 p67 formula: Gamma ions/cm2/s. Applied outside energy validation explicitly.
   yres=54*A**1.18*np.exp(-.78*EV/(KB*T))*(np.maximum(gammacm,1e-99)/1e16)**(-.1)
   contrib['RES_temperature_flux_extrapolation']+=counts*yres
   contrib['RES_only_T_ge1200_extrapolation']+=counts*yres*(T>=1200)
   finecon+=fine['n'][:,j]*yield_fit(fine['eu']*A*1e6,A,Z)
  cc=sum(float((x['n'][:,j]*yield_fit(x['eu']*pm.AA[j]*1e6,float(pm.AA[j]),int(pm.ZZ[j]))).sum()) for j in js)
  out['groups'][G]['physical_fit_calculated_mean_C_atoms_per_ion']=cc/N if N else 0
  out['groups'][G]['physical_fit_calculated_C_mass_kg_m2']=cc*MC
 for key,vals in contrib.items():
  N=float(vals.sum());emtime=float(vals@x['arrival']/N) if N else 0
  fluxC=vals/x['dvweight']*x['arrival_jacobian'];peak_index=int(np.argmax(fluxC))
  out['models'][key]={'peak_sampled_C_atoms_m2_s':float(fluxC[peak_index]),'peak_sampled_time_after_ignition_s':float(x['arrival'][peak_index]),'peak_C_atoms_per_internal_nm2_s_at_area40':float(fluxC[peak_index]/40/1e18),'emitted_C_atoms_m2':N,'emitted_C_kg_m2':N*MC,'emitted_C_mass_ug_m2':N*MC*1e9,'emission_mean_time_after_ignition_s':emtime,'number_weighted_effective_Y_all_ions':N/allcounts.sum()}
  model_rates[key]+=N*MC/period
 c=out['models']['physical_fit_calculated']['emitted_C_atoms_m2'];delta=float(finecon.sum()/c-1);check(stage+' physical-yield 32/64 quadrature',abs(delta)<1e-5,relative_change=delta)
 expected=sum(row['whole_pulse_species_wall_nuclei_m2'][name] for name,z in zip(surf['nuclear_species'],surf['charges']) if z>0)
 check(stage+' charged number closure',abs(allcounts.sum()/expected-1)<1e-12)
 out['physical_fit_32to64_relative_change']=delta
 cno=x['n'][:,np.isin(pm.ZZ,[6,7,8])].sum(axis=1);cnoFlux=cno/x['dvweight']*x['arrival_jacobian'];out['Y_CNO1_peak_sampled_C_atoms_m2_s']=float(cnoFlux.max());out['Y_CNO1_peak_C_atoms_per_internal_nm2_s_at_area40']=float(cnoFlux.max()/40/1e18);out['CNO_mean_arrival_s']=float(cno@x['arrival']/cno.sum());stages[stage]=out

# Exact spherical single-flight radial column. Cosine emission, absorbing next wall.
def hx(x):
 return ((1+x*x)*(-np.log(np.abs(1-x)))-(2*x-x*x)/2)/(x*x)
def quadrature(n):
 z,w=np.polynomial.legendre.leggauss(n);return np.r_[(z+1)/2,(z+3)/2],np.r_[w/2,w/2]
def response(t,mean_eV,n=256):
 t=np.atleast_1d(np.asarray(t,float));ans=np.ones(len(t));ix=np.where(t>0)[0];x,w=quadrature(n);h=hx(x);kt=mean_eV*EV/2
 for part in np.array_split(ix,max(1,math.ceil(len(ix)/256))):
  a=MC/(2*kt)*(R/t[part])**2
  ans[part]=np.sum(w*h*2*a[:,None]**2*x**3*np.exp(-a[:,None]*x*x),axis=1)
 return ans
x,w=quadrature(1024);norm=float(w@hx(x));check('mono-speed impulse integrates to 4R/v',abs(norm/4-1)<1e-6,integral=norm)
# Flux-Maxwell speed law is a scenario kernel; mean emitted E=2kT, no measured graphite energy distribution inferred.
kernels={}
for E in [2*KB*1000/EV,2*KB*1500/EV,1.,3.,10.]:
 T=E*EV/(2*KB);tt=np.geomspace(.002,100,1600)*R/math.sqrt(2*KB*T/MC);hh=response(tt,E,256);hh2=response(tt,E,512);j=int(np.argmax(hh2));err=float(np.max(abs(hh-hh2)));check('kernel quadrature E='+str(E),err<2e-4,max_absolute_difference=err)
 kernels[str(E)]={'mean_ejecta_energy_eV':E,'equivalent_temperature_K':T,'vbar_m_s':vb(T),'steady_column_seconds_per_emission_flux':4*R/vb(T),'single_pulse_peak_column_per_emitted_wall_mass':float(hh2[j]),'peak_time_s':float(tt[j]),'mean_single_flight_lifetime_s':4*R/(3*vb(T))}
# Constant-Y response coefficients are the most reusable result.
coeff={}
for G,Ndot in rate.items():
 q=Ndot*MC;coeff[G]={'incident_nuclei_s':Ndot*AREA,'incident_nuclei_m2_s':Ndot,'C_kg_s_per_Y1':q*AREA,'C_kg_m2_year_per_Y1':q*YEAR,'erosion_um_year_per_Y1_at_rho2000':q*YEAR/2000*1e6,'atomic_C_column_kg_m2_per_Y1_by_mean_energy_eV':{e:q*k['steady_column_seconds_per_emission_flux'] for e,k in kernels.items()}}
q_cno=sum(rate[g] for g in ['C','N','O'])*MC
for name,q in model_rates.items():
 model_rates[name]={'C_kg_s':q*AREA,'C_kg_m2_year':q*YEAR,'erosion_um_year_at_rho2000':q*YEAR/2000*1e6,'atomic_column_by_mean_energy_eV':{e:q*k['steady_column_seconds_per_emission_flux'] for e,k in kernels.items()}}
fast=next(r for r in gas['enthalpy_limits'] if r['mass_speed_quantile']==.999);ideal=gas['ideal_mobile100nm_routing_g0p01']
# Preserve exact baseline object; use the published ideal column only.
base_column=ideal['total_column_kg_m2']
check('ideal baseline species sum',abs(sum(ideal['species_column_kg_m2'])/base_column-1)<1e-14)
thresholds={}
for E,k in kernels.items():
 col=q_cno*k['steady_column_seconds_per_emission_flux'];thresholds[E]={'Y_CNO_for_10pct_ideal_HHe_column':.1*base_column/col,'Y_CNO_for_10pct_fast_total_column_limit':.1*fast['maximum_total_gas_column_kg_m2'][0]/col,'Y_CNO_using_entire_remaining_fast_column_above_ideal_HHe':fast['maximum_added_column_above_ideal100nm_g0p01_kg_m2'][0]/col}

# Finite 370-day dispatch, impulsive spatially-uniform emission; first capture sticks.
events_path=B.parent/'final-family-dispatch/sigma350/events.csv';sources[str(events_path)]=sha(events_path);events=list(csv.DictReader(events_path.open()))
def chronology(dt,E,model):
 end=max(float(r['main_capture_end_s']) for r in events);times=np.arange(0,end+dt,dt);em=np.zeros(len(times));mean_arr={s:(stages[s]['models'][model]['emission_mean_time_after_ignition_s'] if model!='Y_CNO_1' else stages[s]['CNO_mean_arrival_s']) for s in stages}
 for r in events:
  stage=r['stage'];t=float(r['focus_s'])+mean_arr[stage];i=min(int(t/dt+.5),len(times)-1)
  mass=stages[stage]['models'][model]['emitted_C_kg_m2'] if model!='Y_CNO_1' else sum(stages[stage]['groups'][g]['carbon_mass_kg_m2_per_Y1'] for g in ['C','N','O'])
  em[i]+=mass
 kernel=response(times,E,256);nfft=1<<(2*len(times)-1).bit_length();col=np.fft.irfft(np.fft.rfft(em,nfft)*np.fft.rfft(kernel,nfft),nfft)[:len(times)];check('chronology nonnegative '+model+str(dt)+str(E),float(col.min())>-1e-18)
 j=int(col.argmax());return times,col,{'peak_column_kg_m2':float(col[j]),'peak_day':float(times[j]/86400),'mean_column_kg_m2':float(np.mean(col)),'grid_spacing_s':dt,'start':'empty cavity at first dispatch; no earlier ejecta history'}
chron={}
for E in [2*KB*1000/EV,2*KB*1500/EV,1.,3.]:
 for model in (['Y_CNO_1','RES_temperature_flux_extrapolation','RES_only_T_ge1200_extrapolation'] if E<1 else ['physical_fit_calculated','Y_CNO_1']):
  t,c,r=chronology(1800,E,model);_,_,fine=chronology(900,E,model);err=fine['peak_column_kg_m2']/r['peak_column_kg_m2']-1;check('chronology step '+model+str(E),abs(err)<.01,peak_relative_change=err);r['900s_peak_relative_change']=err;r['900s_peak_column_kg_m2']=fine['peak_column_kg_m2'];chron[f'{model}_E{E}']=r
  if E==3. or E==2*KB*1000/EV:np.savetxt(P/(model+'-E'+str(round(E,6))+'-column.csv'),np.c_[t,c],delimiter=',',header='time_s,radial_column_kg_m2',comments='')

# Stable volatile species are separate molecules, with their whole mass and warm-pump capture.
volatile=[]
for species,mass,carbons,G in [('CH4',16,1,'H'),('CO',28,1,'O'),('CO2',44,1,'O'),('HCN',27,1,'N'),('C2N2',52,2,'N')]:
 # Y here means C atoms leaving into cavity per ion of the named group; routing is included.
 Ndot=rate[G];q=Ndot/carbons*mass*MU
 for capture in [.01,.1,1.]:
  col=4*R*q/(capture*vb(1000,mass*MU))
  volatile.append({'molecule':species,'incident_group':G,'capture':capture,'carbon_atoms_per_molecule':carbons,'column_per_C_emitted_per_incident_ion_kg_m2':col,'yield_for_10pct_ideal_HHe_column':.1*base_column/col,'yield_for_entire_remaining_fast_column_above_ideal_HHe':fast['maximum_added_column_above_ideal100nm_g0p01_kg_m2'][0]/col,'scope':'stable molecule diagnostic at1000K; emitted-carbon yield includes liberation and cavity routing. Radioactive O daughters may alter molecule before recapture.'})

# Granularity: representative 2cm C13 pellet and carbon grains; no fracture/dust production model.
ap=.02;rho_core=1011.6428571428571;fill=.8;mcore=4*math.pi/3*ap**3*rho_core*fill;v=fast['speed_m_s'];dust=[]
for ad in [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]:
 md=4*math.pi/3*2000*ad**3;cross=math.pi*(ap+ad)**2;d={'grain_radius_m':ad,'grain_mass_kg':md,'energy_per_impact_J_kg_core_at_fast_speed':.5*md*v*v/mcore,'column_kg_m2_for_expected_1_hit':md/cross,'column_kg_m2_for_1pct_probability_at_least_one_hit':-math.log(.99)*md/cross}
 d['steady_column_kg_m2_for_1e-6_of_incoming_mass_as_dust_by_speed_m_s']={}
 qin=sum(surf['stages'][s]['total_mass_fluence_kg_m2_baryon']/stages[s]['period_s'] for s in stages)
 for vd in [1,10,100,1000]:d['steady_column_kg_m2_for_1e-6_of_incoming_mass_as_dust_by_speed_m_s'][str(vd)]=4*R*qin*1e-6/vd
 dust.append(d)
# RES reference values stay separate from incident-energy dependent physical yields.
res_reference=[{'ion_mass_u':M,'T_K':T,'Gamma_cm2_s':1e14,'C_atoms_per_ion':54*M**1.18*math.exp(-.78*EV/(KB*T))*(1e14/1e16)**(-.1)} for M in [1,4,12,14,16] for T in [1000,1200,1500]]
result={'scope':'Actual128-zone deliveredB1/C13H2 incident spectra and exact350 cadence; incoming normal radial incidence. Physical-yield fits extrapolated to MeV energies and N interpolation are explicit model diagnostics. RES extrapolation uses pulse temperature/flux but is unvalidated in incident energy, isotope decay, damage history and mixed chemistry. AtomicC transport is cosine/flux-Maxwell speed law with absorbing next wall; volatile molecules and dust are separate. No confidence interval assigned to uncalibrated yields.','constants':{'R_m':R,'area_m2':AREA,'carbon_atom_kg':MC,'graphite_rho_kg_m3':2000.,'core_radius_m':ap,'core_mass_kg':mcore,'core_fill':fill,'fast_speed_m_s':v,'ideal_HHe_column_kg_m2':base_column,'fast_enthalpy_cases':fast},'yield_fit_parameters':{'calculated':fits,'experimental_He_C_with_calculated_O_N_interpolation':fit_exp,'N_note':'N Q/Eth are midpoint interpolation between C/O and have no fit measurement; isotope mass enters ETF.','formula':'Y=Q Sn(E/ETF)[1-(Eth/E)^(2/3)](1-Eth/E)^2; Sn=.5ln(1+1.2288eps)/(eps+.1728sqrt(eps)+.008eps^.1504); ETF=30.74(A+12)/12*Z*6*sqrt(Z^(2/3)+6^(2/3)) eV'},'stages':stages,'constant_yield_response_coefficients':coeff,'model_diagnostics':model_rates,'atomic_transport_kernels':kernels,'constant_CNO_yield_thresholds_unit_sticking':thresholds,'finite_dispatch_atomic_column_scenarios':chron,'volatile_carbon_scenarios':volatile,'dust_granularity_scenarios':dust,'dust_incoming_baryon_mass_flux_kg_m2_s':qin,'RES_reference_diagnostic':res_reference,'checks':checks,'source_sha256':sources,'code_sha256':sha(__file__),'RES_temperature_scope':'Each stage reuses the warmest exported150um-skin profile; most B1 pulses are cooler. All-T expression and optional T>=1200K mask are separate unsupported extrapolations, not bounds or predictions.'}
sp.write(P/'results.json',result)
print(json.dumps({k:result[k] for k in ['constant_yield_response_coefficients','model_diagnostics','constant_CNO_yield_thresholds_unit_sticking','finite_dispatch_atomic_column_scenarios']},indent=2));print(len(checks),'checks passed')
