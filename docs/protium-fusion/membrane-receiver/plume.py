"""Thermal evaporation and first-flight plume screen; no sputtering model.

Graphite equilibrium partial pressures come from NIST-JANAF C1/C2/C3 tables.
Source temperature histories are reconstructed from the saved receiver inputs,
keeping microsecond prompt heating separate from slow debris/decay heating.
"""
import json,math,hashlib,sys,argparse
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import roots_genlaguerre
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OLD=ROOT/'open-loop-redesign'
KB=1.380649e-23;NA=6.02214076e23;RG=8.314462618;MU=1.66053906660e-27
R=1.1e9;AREA=4*math.pi*R*R;K=100.;RHO=2000.;CP=1000.;ALPHA=K/(RHO*CP);L=2.5;YEAR=365.25*86400
# Manually transcribed primary table values, 1000..2100K at100K intervals.
# Formation reaction nC(graphite)->Cn(g): pressure/bar =10**logKf.
TABLE_T=np.arange(1000.,2101.,100.)
LOGK=np.array([
[-29.286,-25.869,-23.023,-20.615,-18.552,-16.765,-15.202,-13.823,-12.598,-11.503,-10.517,-9.626],
[-33.640,-29.652,-26.332,-23.526,-21.124,-19.044,-17.227,-15.625,-14.203,-12.933,-11.791,-10.760],
[-31.400,-27.544,-24.338,-21.633,-19.320,-17.322,-15.578,-14.044,-12.685,-11.472,-10.384,-9.402]])
HF=np.array([
[719.475,719.360,719.188,718.968,718.709,718.415,718.092,717.744,717.373,716.984,716.577,716.156],
[840.132,839.353,838.484,837.548,836.559,835.528,834.465,833.375,832.266,831.140,830.000,828.849],
[813.156,811.105,808.974,806.781,804.538,802.255,799.940,797.595,795.226,792.834,790.421,787.989]])
CMASS=np.arange(1,4)*12.0107e-3/NA
SOURCES=[
{'url':'https://janaf.nist.gov/tables/C-003.html','use':'C gas logKf and formation enthalpy versus graphite; standardpressure1bar.'},
{'url':'https://janaf.nist.gov/tables/C-113.html','use':'C2 gas logKf and formation enthalpy versus graphite.'},
{'url':'https://janaf.nist.gov/tables/C-138.html','use':'C3 gas logKf and formation enthalpy versus graphite.'},
{'url':'https://doi.org/10.1063/1.1743248','use':'Thorn/Winslow1957 equilibrium vapor and evaporation coefficients; individual C/C2/C3 coefficients .37/.34/.08 around2400K are a sensitivity here, not measured at1000–2000K.'},
{'url':'https://arxiv.org/abs/2002.07835','use':'Long etal2020 measured evolving nanoparticle sublimation initially about400times bulk; motivates a rough/damaged-surface sensitivity, not an upperbound for irradiated receiver.'},
{'url':'https://pubs.rsc.org/en/content/articlelanding/1973/f1/f19736901628','use':'Nagai etal1973 silica effusion1823–1983K: Gibbsreaction182050−58.29T cal/mol; vaporizationcoefficient.022±.008. SimpleSiO/O2 approximation here omits atomicO and other species.'},
{'url':'https://doi.org/10.3847/0004-637X/824/2/103','use':'Fegley etal2016 thermochemical vapor calculation gives1.54Pa total at2000K, predominantlySiO/O2/O/SiO2.'},
{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6716270/','use':'Plante/Sessoms1973 tungsten Langmuir data2600–3100K,log10P(atm)=7.871−45385/T; lowerT extrapolation labeled.'}]

def graphite(T,coeff=None):
 T=np.asarray(T);assert T.min()>=1000 and T.max()<=2100,(T.min(),T.max())
 # Linear interpolation in1/T of logP approximates local Clapeyron behavior.
 pressures=np.array([1e5*10**np.interp(1/T,1/TABLE_T[::-1],row[::-1]) for row in LOGK])
 j=pressures*np.sqrt(CMASS.reshape((-1,)+(1,)*T.ndim)/(2*math.pi*KB*T))
 if coeff is not None:j=j*np.asarray(coeff).reshape((-1,)+(1,)*T.ndim)
 h=np.array([np.interp(T,TABLE_T,row) for row in HF])*1000/(CMASS.reshape((-1,)+(1,)*T.ndim)*NA)
 return pressures,j,h

def slab_step(t):
 t=np.asarray(t);out=2*np.sqrt(ALPHA*t)/(K*np.sqrt(math.pi));late=t>2000
 if np.any(late):
  lam=ALPHA*((np.arange(80)+.5)*math.pi/L)**2
  out[late]=L/K-(2*ALPHA/(K*L))*np.sum(np.exp(-lam[:,None]*t[late])/lam[:,None],axis=0)
 return out

def thermal_flux(item,dt=.25,end=200000):
 d=item['record'];n=int(end/dt);t=(np.arange(n)+.5)*dt
 v=np.r_[0.,d['profile']['velocity_m_s']];dm=d['mass_kg']*np.array(d['profile']['mass_fraction_of_total']);slope=dm/np.diff(v)
 ke=np.sum(slope*np.diff(v**3)/6);slope*=d['final_kinetic_J']/ke
 speed=R/t;i=np.searchsorted(v,speed,side='right')-1;ok=(i>=0)&(i<len(dm));q=np.zeros(n)
 q[ok]=.5*speed[ok]**4/R*slope[i[ok]]/AREA*(1+d['final_internal_J']/d['final_kinetic_J'])
 # Decay lifetimes and emitted non-neutrino energies already frozen in receiver.
 stage=SCREENS[d['stage_cat']];tlight=R/299792458.;delay_t=np.maximum(t-tlight,0)
 rates={'n13':math.log(2)/597.9,'o14':math.log(2)/70.606,'o15':math.log(2)/122.24}
 # Use exact same library rates from frozen energy/flux metadata when available.
 for parent,energy in stage['delayed_decay_energy_J'].items():
  rate=RATES[parent];q+=np.where(t>=tlight,energy/AREA*rate*np.exp(-rate*delay_t),0)
 weights=np.diff(slab_step(np.arange(n+1)*dt))
 rise=fftconvolve(q,weights)[:n]
 # Late prompt heat response: first0.1s treated on a separate fine grid below.
 age=np.maximum(t-tlight,1e-30);tail=d['radiation_escaped_J']/AREA*np.sqrt(ALPHA)/(K*np.sqrt(math.pi*age));tail[t<tlight+.1]=0
 rise+=tail
 return t,rise,q

def prompt_history(item):
 tt=np.array(item['prompt_time_s']);ee=np.array(item['prompt_cumulative_J'])/AREA
 end=tt[-1];fine=np.linspace(0,end,12001);delta=fine[1];flux=np.diff(np.interp(fine,tt,ee))/delta
 n=np.arange(len(flux));weights=2*np.sqrt(ALPHA*delta)/(K*np.sqrt(math.pi))*(np.sqrt(n+1)-np.sqrt(n))
 rise=fftconvolve(flux,weights)[:len(flux)];t=fine[1:]
 # Exact piecewiseconstant original trace convolution after prompt, logarithmically sampled.
 tailt=np.geomspace(end*(1+1e-9),.1,1600)
 slopes=np.diff(ee)/np.diff(tt);later=np.zeros_like(tailt)
 for a,b,s in zip(tt[:-1],tt[1:],slopes):later+=s*(np.sqrt(np.maximum(tailt-a,0))-np.sqrt(np.maximum(tailt-b,0)))
 later*=2*np.sqrt(ALPHA)/(K*np.sqrt(math.pi))
 return np.r_[0,t,tailt],np.r_[0,rise,later]

def integrated_emission(t,temp):
 p,j,h=graphite(temp);fluence=np.trapezoid(j,t,axis=1);latent=np.trapezoid(j*h,t,axis=1)
 weighted=np.trapezoid(j*temp,t,axis=1)/np.maximum(fluence,1e-300)
 return {'mass_kg_m2_by_Cn':fluence.tolist(),'total_mass_kg_m2':float(sum(fluence)),
 'total_mass_kg':float(sum(fluence)*AREA),'latent_J_m2':float(sum(latent)),
 'mass_weighted_T_by_Cn_K':weighted.tolist(),'peak_T_K':float(max(temp)),
 'peak_equilibrium_P_Pa':float(max(p.sum(axis=0))),'peak_mass_flux_kg_m2_s':float(max(j.sum(axis=0)))}

def sphere_plume(T,m):
 x,w=roots_genlaguerre(160,1);v=np.sqrt(2*KB*T*x/m)
 def survival(t):a=v*t/R;return float(w@np.maximum(1-a*a/4,0))
 def column(t):
  a=v*t/R;ok=a<2;a=a[ok];ww=w[ok];b=1+a*a;u=a/2
  integral=-(1-u)/a-b/(2*a*a)*np.log((b-2*a)/(b-2*a*u))
  # Floating cancellation atveryearlytimes notused inreportedrange.
  return float(ww@integral)
 vbar=math.sqrt(8*KB*T/(math.pi*m));mean=4*R/(3*vbar)
 return {'thermal_gas_mean_speed_m_s':vbar,'mean_first_wall_return_s':mean,
 'return_time_quantiles_s':{str(f):brentq(lambda t:survival(t)-(1-f),1,1e8) for f in [.01,.5,.9,.99]},
 'survival_at_delays':{str(t):survival(t) for t in [6774.31305795,86400.,604800.,2592000.]},
 'radial_column_per_emitted_surface_mass':{str(t):column(t) for t in [6774.31305795,86400.,604800.,2592000.]},
 'scope':'Diffusecosinelaw emission and fluxMaxwell speeds, nointermolecularcollisions, firstwallhitsticks. Spheregravity cancels inside; magneticfieldsdo notremove neutrals.'}

def alternatives(T):
 # SiO2->SiO+.5O2, pSiO=2pO2 equilibriumstoichiometry approximation.
 dg=(182050-58.29*T)*4.184;k=math.exp(-dg/(RG*T));psi=1.01325e5*(math.sqrt(2)*k)**(2/3);po=psi/2
 jsi=.022*(psi*math.sqrt(.0440849/(2*math.pi*RG*T))+po*math.sqrt(.0319988/(2*math.pi*RG*T)))
 pw=101325*10**(7.871-45385/T);jw=pw*math.sqrt(.18384/(2*math.pi*RG*T))
 return {'T_K':T,'silica_approx_total_P_Pa':psi+po,'silica_approx_flux_kg_m2_s_alpha022':jsi,
 'silica_scope':'1823–1983Kmeasurement; outsidethisrangeextrapolation; SiO/O2only andequilibriumstoichiometry approximateopenvacuumcomposition.',
 'tungsten_P_Pa':pw,'tungsten_flux_kg_m2_s_alpha1':jw,'tungsten_scope':'Extrapolates2600–3100Kprimaryfit; noirradiationorchemicalenhancement.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--dt',type=float,default=.25);ap.add_argument('--output',default='results.json');args=ap.parse_args()
 inputs=json.loads((OLD/'collector-inputs.json').read_text());station=json.loads((OLD/'station-inputs.json').read_text());periodic=json.loads((OLD/'collector-periodic.json').read_text());SCREENS=station['receiver']['stages']
 sys.path.insert(0,str(ROOT));from nuclear_screen import LIB
 RATES={p:LIB.get_rate_by_name(p+'(,)'+d).eval(1e8) for p,d in [('n13','c13'),('o14','n14'),('o15','n15')]}
 rate=.8/(station['receiver']['occupied_s_per_He_at_R1e9']['mass999']*R/1e9);intervals={k:1/(rate*p['pulses_per_He']) for k,p in station['cycle']['stage_pulse_data'].items()}
 history=periodic['accumulated_background_upper_K'];b1=periodic['isolated_stage_peak_bounds_K']['n15'];rows={}
 for cat,item in inputs['stages'].items():
  item['record']['stage_cat']=cat;t,rise,q=thermal_flux(item,dt=args.dt);pt,pr=prompt_history(item)
  # Constant background+worst B1maximum throughoutlargepulse is conservative.
  background=1000+history+(b1 if cat!='n15' else 0)
  slow=integrated_emission(t,background+rise)
  initial_beta_flux=sum(energy/AREA*RATES[parent] for parent,energy in SCREENS[cat]['delayed_decay_energy_J'].items())
  fast=integrated_emission(pt,background+pr+2*initial_beta_flux*np.sqrt(ALPHA*pt)/(K*np.sqrt(math.pi)))
  late_T=background+rise[-1]+q[-1]*L/K
  late_duration=max(intervals[cat]-t[-1],0)
  late_mass=graphite(np.asarray(late_T))[1]*late_duration
  total=np.array(slow['mass_kg_m2_by_Cn'])+np.array(fast['mass_kg_m2_by_Cn'])+late_mass
  rows[cat]={'mean_interval_s':intervals[cat],'background_K':background,'slow':slow,'prompt':fast,
   'total_mass_kg_m2_by_Cn':total.tolist(),'total_mass_kg_m2':float(total.sum()),'total_mass_kg':float(total.sum()*AREA),
   'annual_gross_evaporation_kg_m2':float(total.sum()*YEAR/intervals[cat]),'annual_gross_evaporation_kg':float(total.sum()*AREA*YEAR/intervals[cat]),
   'assumed_bulk_coefficients_sensitivity_mass_kg':float(total@np.array([.37,.34,.08])*AREA),
   'all_species_400x_surface_enhancement_mass_kg':float(total.sum()*400*AREA),
   'latent_J_m2':slow['latent_J_m2']+fast['latent_J_m2'],
   'late_window_constant_upper_T_K':float(late_T),'late_window_s':float(late_duration),'late_window_mass_upper_kg':float(late_mass.sum()*AREA),
   'late_window_logic':'Pastpositiveheatcontributionsdecay becausekerneldecreases; futuret^-4/decayforcing boundedbycurrentflux timesL/k. ConstantworstB1backgroundovercounts. Latentfigure excludes negligiblelatewindow.',
   'thermal_grid_dt_s':args.dt,'thermal_end_s':float(t[-1]),
   'thermal_scope':'Finitecooledslabslowresponse; exactfastpromptconvolution; separatelyintegrated0-.1spromptandslowhistorydoublecountsnegligiblebackground. Assumedpropertiesanddiffusionmodelretained.'}
 annual=sum(z['annual_gross_evaporation_kg'] for z in rows.values());allow=station['cycle']['fresh_hydrogen_kg_per_He']*rate*.0014
 propertyrows=[]
 for T in [1000.,1200.,1400.,1600.,1800.,1854.,1900.,1965.326,2000.,2100.]:
  p,j,h=graphite(T);propertyrows.append({'T_K':T,'partial_P_Pa':p.tolist(),'total_P_Pa':float(p.sum()),'mass_flux_kg_m2_s_alpha1':float(j.sum()),'mass_flux_kg_m2_s_bulk_coefficients':float(j@np.array([.37,.34,.08])),'latent_W_m2_alpha1':float(sum(j*h))})
 plumes={str(i+1):sphere_plume(1700,CMASS[i]) for i in range(3)}
 # Fresh-C loss and pollution requirements differ. Escape fractionsareexplicitvariables.
 opening=[]
 for f in [0.,1e-6,1e-3,.01,1.]:opening.append({'irreversible_escape_fraction':f,'fresh_C_kg_s':annual/YEAR*f,'fraction_material_budget':annual/YEAR*f/allow,'annual_irreversible_kg_m2':annual/AREA*f})
 contamin=[]
 maxsig=max(r['total_mass_kg_m2'] for r in rows.values())
 for sigma in [1e-20,1e-19,1e-18]:contamin.append({'assumed_effective_interaction_cross_section_m2':sigma,'max_single_pulse_atom_column_m2':3*maxsig/(12.0107*MU),'Poisson_encounter_optical_depth_screen':3*maxsig/(12.0107*MU)*sigma,'scope':'Factor3spatialcolumnscreen, allcarboncountedasatoms; crosssections areparametersnot measuredvalues; nofuel-grainshielding credit.'})
 # For a well-mixed gas with cold sticking surface fractionf, S=fAvbar/4.
 cond=[]
 for f in [1.,.1,.01,.001]:
  v=plumes['1']['thermal_gas_mean_speed_m_s'];S=f*AREA*v/4;tau=(4*math.pi*R**3/3)/S
  cond.append({'cold_condensing_area_fraction':f,'pumping_speed_m3_s':S,'well_mixed_clearing_time_s':tau,'average_condensation_mass_kg_s':annual/YEAR,'scope':'Noactiveionization; molecularflow andunitsticking. Actualnearwallcapturecanbequicker withshortgapgeometry.'})
 out={'scope':__doc__,'radius_m':R,'area_m2':AREA,'source_tables':{'T_K':TABLE_T.tolist(),'log10_Kf_C1_C2_C3':LOGK.tolist(),'formation_enthalpy_kJ_mol_C1_C2_C3':HF.tolist(),'interpolation':'logPlinear1/T; hlinearT; p°1bar; idealgasactivity; puregraphiteactivity1'},
 'primary_sources':SOURCES,'source_access_note':'NISTHTMLtablesretrievedthroughwebtool; directtxtdownloadreturned403 andwasnotretried. Numericcolumnsaboveareexplicittranscription.',
 'graphite_properties':propertyrows,'pulse_cases':rows,'annual_gross_evaporation_kg':annual,'average_gross_evaporation_kg_s':annual/YEAR,'annual_gross_evaporation_kg_m2':annual/AREA,'material_C_allowance_kg_s':allow,
 'loss_fraction_cases':opening,'diffuse_firstflight_plumes_at1700K':plumes,'interaction_column_screens':contamin,'condensation_area_screens':cond,
 'alternative_membrane_materials': [alternatives(t) for t in [1000.,1500.,1800.,1854.,1965.326,2000.]],
 'limitations':['Eachstageis integratedfor200ks; its longermeaninterval is boundedbyconstantterminaluppertemperature. OverlappinglowB1backgroundwindows deliberatelyovercount ratherthan resetthewall.',
 'Noimpactdrivensputtering,desorption,fragmentation,thermalstress,dustorchemicalerosionmodel. These can exceed equilibrium evaporation.',
 'Bulk JANAFphaseequilibriumisnotanupperboundfor a metastable irradiated/nanostructured surface;400xsensitivityisillustrative.',
 'Uniformwallthermalhistoryandconstantpropertiesare retained; sourcegeometrycouldcreatehotspots.',
 'Promptsurfaceabsorptionandenergyarrivalmodelremainreducedapproximations; realphotondepthcanreducepeakskinheating.',
 'Condensation coefficient, grain collisioncrosssections and pulse contamination tolerance remainengineeringinputs.',
 'Silica/Wcomparisonsrequire actualmembranetemperature; athinfilmcanbe muchhotterthanmassivegraphite.',
 'The factor3singlepulsecolumnis ascreen,notarigorousalltimesbound; analyticdiffusespherecolumnfactorsareprovided separately.',
 'Microscopicvacuumcollisionratesdo notdirectlyequate to bulkfuelmixing or ignitionfailure.'],
 'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD/'collector-inputs.json',OLD/'station-inputs.json',OLD/'collector-periodic.json']}}
 (HERE/args.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'pulses':{k:{q:z[q] for q in ['total_mass_kg','annual_gross_evaporation_kg','latent_J_m2']}|{'slow_peak':z['slow']['peak_T_K'],'prompt_peak':z['prompt']['peak_T_K']} for k,z in rows.items()},'annual_gross_kg':annual,'mean_kg_s':annual/YEAR,'material_allowance_kg_s':allow,'plume_C1':plumes['1']},indent=2))
