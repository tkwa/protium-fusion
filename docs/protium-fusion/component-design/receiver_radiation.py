#!/usr/bin/env python3
"""Bounded photon-deposition, metal-vapor and oxide-equilibrium interface screen.
No new hydrodynamics; reads frozen corrected64 pulse histories.
Run: python receiver_radiation.py --root ../ --output receiver-radiation.json
"""
import argparse,importlib.util,json,hashlib,math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.signal import fftconvolve
RG=8.314462618;NA=6.02214076e23;R=1.1e9;AREA=4*math.pi*R**2
# Primary NIST XCOM, keV and cm2/g. Repeated energy rows retain actual edges.
TABLES={
'Al':'''1 1185 1183
1.5 402.2 400.1
1.5596 362.1 360
1.5596 3957 3829
2 2263 2204
3 788 773.2
4 360.5 354.5
5 193.4 190.2
6 115.3 113.3
8 50.33 49.18
10 26.23 25.43
15 7.955 7.487
20 3.441 3.094
30 1.128 .8778
40 .5685 .3601
50 .3681 .184
60 .2778 .1099
80 .2018 .05511
100 .1704 .03794
150 .1378 .02827
200 .1223 .02745
300 .1042 .02816
400 .09276 .02862
500 .08445 .02868
600 .07802 .02851
800 .06841 .02778
1000 .06146 .02686''',
'Sn':'''1 8157 8144
1.5 3296 3287
2 1665 1657
3 614.3 607.7
3.9288 311.4 305.8
3.9288 928.5 884.5
4 939.3 894.5
4.1561 846.9 807.4
4.1561 1145 1084
4.30764 1060 1002
4.4647 971.2 922.1
4.4647 1117 1059
5 847.1 806.4
6 529.4 506.5
8 250 240
10 138.4 132.7
15 46.64 44.13
20 21.46 19.90
29.2001 7.760 6.876
29.2001 43.60 15.14
30 41.21 14.90
40 19.42 9.835
50 10.70 6.314
60 6.564 4.211
80 3.029 2.101
100 1.676 1.189
150 .6091 .4119
200 .3260 .1990
300 .1639 .08179
400 .1156 .05138
500 .09374 .0398
600 .08113 .03417
800 .06662 .02878
1000 .058 .02601'''}
TABLES={k:np.array([list(map(float,r.split())) for r in v.splitlines()]) for k,v in TABLES.items()}
AL_T=np.arange(500,2101,100.)
AL_LOGK=np.array([-27.348,-21.629,-17.552,-14.502,-12.137,-10.292,-8.805,-7.570,-6.529,-5.640,-4.872,-4.203,-3.614,-3.093,-2.628,-2.211,-1.836])
# JANAF formation Gibbs energies kJ/mol, at listed temperatures.
GT=np.array([1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000.])
G={
'Al':[197.027,185.422,173.917,162.503,151.172,139.921,128.742,117.631,106.585,95.600,84.673],
'AlO':[-15.315,-21.704,-27.988,-34.182,-40.296,-46.343,-52.331,-58.269,-64.164,-70.023,-75.850],
'AlO2':[-102.421,-102.681,-102.885,-103.034,-103.133,-103.185,-103.191,-103.153,-103.074,-102.953,-102.794],
'Al2O':[-230.375,-235.573,-240.587,-245.433,-250.124,-254.670,-259.082,-263.368,-267.534,-271.587,-275.532],
'Al2O2':[-407.314,-405.518,-403.565,-401.467,-399.235,-396.877,-394.400,-391.811,-389.117,-386.321,-383.428],
'O':[187.681,181.165,174.619,168.047,161.453,154.840,148.210,141.564,134.905,128.234,121.552],
'O2':[0]*11,
'Al2O3cr':[-1361.437,-1328.286,-1295.228,-1262.264,-1229.393,-1196.617,-1163.934,-1131.342,-1098.841,-1066.426,-1034.096]}
STO={'Al':(1,0),'AlO':(1,1),'AlO2':(1,2),'Al2O':(2,1),'Al2O2':(2,2),'O':(0,1),'O2':(0,2)}
M={s:(a*26.9815385+b*15.999)*1e-3 for s,(a,b) in STO.items()}
SOURCES={'Al_XCOM':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html','Sn_XCOM':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z50.html','Al_gas':'https://janaf.nist.gov/tables/Al-005.html','Al_liquid':'https://janaf.nist.gov/tables/Al-004.html','AlO_gas':'https://janaf.nist.gov/tables/Al-074.html','AlO2_gas':'https://janaf.nist.gov/tables/Al-077.html','Al2O_gas':'https://janaf.nist.gov/tables/Al-092.html','Al2O2_gas':'https://janaf.nist.gov/tables/Al-094.html','Al2O3_alpha':'https://janaf.nist.gov/tables/Al-096.html','O_gas':'https://janaf.nist.gov/tables/O-001.html'}
def loadmod(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def al_vapor(T):
 T=np.asarray(T);assert T.min()>=500 and T.max()<=2100
 p=1e5*10**np.interp(1/T,1/AL_T[::-1],AL_LOGK[::-1]);j=p*np.sqrt(M['Al']/(2*math.pi*RG*T))
 return p,j

def oxide_equilibrium(T):
 # Surface chemical equilibrium plus congruent effusive atomic flux O/Al=3/2.
 # x and y are elemental chemical potentials/RT, relative to pure Al and half O2.
 g={s:np.interp(T,GT,v)*1000/(RG*T) for s,v in G.items()}
 def pressures(y):
  x=(g['Al2O3cr']-3*y)/2
  return {s:math.exp(a*x+b*y-g[s])*1e5 for s,(a,b) in STO.items()}
 def f(y):
  p=pressures(y);al=sum(STO[s][0]*p[s]/math.sqrt(M[s]) for s in STO);o=sum(STO[s][1]*p[s]/math.sqrt(M[s]) for s in STO)
  return math.log(o/al)-math.log(1.5)
 y=brentq(f,-200,30);p=pressures(y);j={s:v*math.sqrt(M[s]/(2*math.pi*RG*T)) for s,v in p.items()}
 ratio=sum(STO[s][1]*j[s]/M[s] for s in STO)/sum(STO[s][0]*j[s]/M[s] for s in STO)
 return {'temperature_K':T,'partial_pressure_Pa':p,'total_pressure_Pa':sum(p.values()),'mass_flux_kg_m2_s':sum(j.values()),'species_mass_flux':j,'congruent_flux_O_per_Al':ratio}

def coeff(mat,e):
 t=TABLES[mat].copy()
 for i in range(len(t)-1):
  if t[i,0]==t[i+1,0]:t[i,0]*=(1-1e-10)
 return np.array([np.exp(np.interp(np.log(e),np.log(t[:,0]),np.log(t[:,i]))) for i in [1,2]])*.1

def main(root):
 root=Path(root);mem=root/'membrane-receiver';old=root/'open-loop-redesign'
 sp=loadmod(mem/'spectral_proxy.py','saved_spectral');pl=loadmod(mem/'plume.py','saved_plume')
 pulse=json.loads((mem/'pulse-loading.json').read_text());inputs=json.loads((old/'collector-inputs.json').read_text());station=json.loads((old/'station-inputs.json').read_text());periodic=json.loads((old/'collector-periodic.json').read_text());oldplume=json.loads((mem/'results.json').read_text())
 pl.SCREENS=station['receiver']['stages'];pl.RATES={p:math.log(2)/t for p,t in [('n13',597.9),('o14',70.606),('o15',122.24)]}
 matprops={'Al':{'rho':2375.,'cp':31.751/.0269815385,'k':80.,'scope':'JANAF liquid cp; rho,k are explicit engineering values'},'Al2O3':{'rho':3900.,'cp':124.771/.101960077,'k':10.,'scope':'JANAF alpha crystal1000K cp; density and conductivity are engineering values'},'Sn':{'rho':6700.,'cp':250.,'k':30.,'scope':'All thermophysical values are illustrative parameters pending qualification'}}
 edges=np.concatenate([t[:,0] for t in TABLES.values()]);e=np.unique(np.r_[np.geomspace(1,1000,6500),edges,edges*(1-1e-10)])
 mc={m:coeff(m,e) for m in ['Al','Sn']};mc['Al2O3']=mc['Al']*(2*.0269815385/.101960077)+sp.coeff('O',e)*(3*.015999/.101960077)
 checks=[];flash={}
 for stage,s in pulse['stages'].items():
  pr=s['prompt'];df=np.diff(pr['cumulative_fluence_J_m2']);temp=np.array(pr['outer_temperature_K']);kt=.5*(temp[:-1]+temp[1:])*sp.KB_KEV
  cases={}
  for fac in [.5,1.,2.]:
   spec,low,hi=sp.planck_parts(kt*fac,df,e);rows={}
   # Unresolved sub-keV energy placed entirely at the surface with its original time bins.
   lowgrid=np.geomspace(1e-7,1,1800);lowbins=np.trapz(sp.spectrum(lowgrid,kt*fac),lowgrid,axis=1)*df
   times=np.array(pr['emission_time_s']);times=times-times[0];fine=np.linspace(0,times[-1],12001);step=fine[1]
   lowflux=np.diff(np.interp(fine,times,np.r_[0,np.cumsum(lowbins)]))/step
   squareweights=np.sqrt(np.arange(len(lowflux))+1)-np.sqrt(np.arange(len(lowflux)))
   for mat,p in matprops.items():
    mu,men=mc[mat];sig=p['rho']*1e-4;loss=-np.expm1(-sig*mu);dep=men/mu*loss
    trans=float(np.trapz(spec*np.exp(-sig*mu),e));heat=float(np.trapz(spec*dep,e));surf=float(np.trapz(spec*men,e)/p['cp'])
    checks.extend([0<=heat<=sum(df),0<=trans<=sum(df),surf>=heat/(sig*p['cp'])])
    rows[mat]={'thickness_m':1e-4,'areal_mass_kg_m2':sig,'resolved_transfer_J_m2':heat,'resolved_uncollided_transmission_J_m2':trans,'resolved_surface_sensible_heat_equivalent_K_no_conduction':surf,'layer_mean_sensible_heat_equivalent_K':heat/(sig*p['cp']),'unresolved_low_high_full_retention_J_m2':low+hi,'below1keV_all_surface_conduction_peak_rise_K':float(max(fftconvolve(lowflux,squareweights)[:len(lowflux)])*2*math.sqrt(step)/(math.sqrt(math.pi*p['k']*p['rho']*p['cp']))),'scope':'mu_en first-interaction energy transfer with bulk local retention. Sensible equivalent precedes phase change; no secondary electron or repeated-fluorescence transport.'}
   cases[str(fac)]={'temperature_multiplier':fac,'below1keV_J_m2':low,'above1MeV_J_m2':hi,'rows':rows}
  flash[stage]={'total_prompt_J_m2':float(sum(df)),'Planck_shape_cases':cases}
 # Preserve old separated temperature histories, no constant peak-time replacement.
 vaporrows={};dt=.5
 oxide_grid=np.linspace(1000,2000,201);oxide_j=np.array([oxide_equilibrium(float(T))['mass_flux_kg_m2_s'] for T in oxide_grid])
 def oxflux(T):return np.exp(np.interp(1/np.asarray(T),1/oxide_grid[::-1],np.log(oxide_j[::-1])))
 for cat,item in inputs['stages'].items():
  item['record']['stage_cat']=cat;t,rise,q=pl.thermal_flux(item,dt=dt);pt,pr=pl.prompt_history(item)
  background=oldplume['pulse_cases'][cat]['background_K'];slowT=background+rise
  # Conditional backed metal follows carbon temperature; adds no contact jump here.
  pp,j=al_vapor(slowT);basej=float(al_vapor(np.array(background))[1]);dM=float(np.trapz(j-basej,t))
  fp,fj=al_vapor(background+pr);fast=float(np.trapz(fj-basej,pt))
  # Baseline and excess separate; no peak-to-duration extrapolation.
  interval=oldplume['pulse_cases'][cat]['mean_interval_s']
  vaporrows[cat]={'interval_s':interval,'background_K':background,'slow_peak_K':float(max(slowT)),'slow_excess_Al_kg_m2_per_pulse':dM,'prompt_excess_using_old_all_surface_carbon_history_kg_m2_per_pulse':fast,'mean_excess_Al_source_kg_s':(dM+fast)*AREA/interval,'oxide_excess_mass_kg_m2_same_thermal_history_plus10K':float(np.trapz(oxflux(slowT+10)-oxflux(background+10),t)+np.trapz(oxflux(background+pr+10)-oxflux(background+10),pt)),'peak_Al_vapor_Pa':float(max(pp)),'peak_Al_mass_flux_kg_m2_s':float(max(j)),'slow_end_s':float(t[-1]),'scope':'200ks excess integral, small late tail omitted; unchanged carbon temperature history is conditional, ignores evaporation heat feedback and altered material heat transfer. Prompt is a deliberately old all-surface upper reference, not the Planck flash model.'}
 gasrows=[]
 for T in [700.,750.,800.,900.,933.45,1000.,1100.,1200.,1400.]:
  p,j=map(float,al_vapor(np.array(T)));vbar=math.sqrt(8*RG*T/(math.pi*M['Al']));sig=p*M['Al']*R/(RG*T)
  checks.append(abs(4*j*R/vbar/sig-1)<1e-12)
  gasrows.append({'T_K':T,'P_Pa':p,'flux_kg_m2_s':j,'whole_sphere_gross_source_kg_s':AREA*j,'equilibrium_radial_mass_column_kg_m2':sig,'first_wall_return_mean_s':4*R/(3*vbar)})
 sat1000=next(z for z in gasrows if z['T_K']==1000)['equilibrium_radial_mass_column_kg_m2'];budgets=[]
 for a in [1e-4,1e-3,1e-2]:
  for u in [1000.,10000.,100000.]:
   sig=8*1000*a*u/(3*(2.4e6)**2)
   # Congruent equilibrium of cold Al deposit gives a conservative residual floor.
   trap=brentq(lambda T:float(al_vapor(np.array(T))[0])*M['Al']*R/(RG*T)-sig,500,1000)
   budgets.append({'grain_radius_m':a,'allocated_gas_heating_J_kg':u,'energy_retention_fraction':1.,'allowed_radial_gas_column_kg_m2':sig,'maximum_hot_area_escape_over_effective_cold_area_at1000K':sig/sat1000,'cold_Al_deposit_equilibrium_temperature_limit_K':trap,'ideal_single_aperture_depth_over_radius_if_coldwalls_unitstick':math.sqrt(sat1000/sig-1)})
 ox=[oxide_equilibrium(float(T)) for T in GT];checks.extend(abs(z['congruent_flux_O_per_Al']-1.5)<1e-8 for z in ox)
 hper=station['cycle']['returned_fuel_atoms_per_He']['p']+station['cycle']['ash_atoms_per_He']['p']
 hflow=hper*periodic['He_rate_s']*1.67262192369e-27;hgas=[]
 for label,molar in [('H',1.00784e-3),('H2',2.01568e-3)]:
  speed=math.sqrt(8*RG*1000/(math.pi*molar));column=4*hflow*R/(AREA*speed)
  hgas.append({'species':label,'T_K':1000.,'postburn_protons_per_He':hper,'hydrogen_mass_flow_kg_s':hflow,'mean_speed_m_s':speed,'unit_first_return_capture_radial_column_kg_m2_if_all_released':column,'escape_fraction_over_global_capture_fraction_limit_for_1cm_10kJkg':8*1000*.01*10000/(3*(2.4e6)**2)/column,'scope':'Neutral gas molecular-flow mean reference; effective global capture must actually exist. Multiply column by escape fraction / capture fraction. H2 does not condense at1000K.'})
 dust=[]
 for a in [1e-3,1e-2]:
  energy_budget=10000.;retention=1.;grainrho=1000.;dustrho=3000.;vrel=2.4e6;pmax=.001
  b=a*(2*energy_budget*grainrho/(retention*vrel*vrel*dustrho))**(1/3)
  number=-math.log1p(-pmax)/(math.pi*(a+b)**2);mass=number*4*math.pi*b**3*dustrho/3
  dust.append({'fuel_grain_radius_m':a,'dust_radius_m_for10kJkg_single_hit':b,'allowed_fraction_grains_above_energy_budget':pmax,'critical_size_dust_radial_column_limit_kg_m2':mass,'uniform_dust_source_kg_s_at100m_s_unit_first_return_capture':mass*AREA*100/(4*R),'scope':'Poisson projected-area interception, one catastrophic-equivalent hit; a grain-damage fraction, not the probability an entire macropulse has any hit. Dust size/speed/energy-retention are assumptions; no source yield is inferred.'})
 files=[mem/'pulse-loading.json',mem/'spectral_proxy.py',mem/'plume.py',mem/'results.json',old/'collector-inputs.json',old/'station-inputs.json',old/'collector-periodic.json']
 return {'scope':'Conditional material-interface screen. Corrected64 pulse inputs; frozen data untouched. Primary atomic attenuation and equilibrium thermochemistry plus explicit reduced transport assumptions.','primary_sources':SOURCES,'material_parameters':matprops,'flash':flash,'Al_gas_properties':gasrows,'Al_backed_surface_pulse_excess':vaporrows,'Al1000K_baseline_source_kg_s':next(z for z in gasrows if z['T_K']==1000)['whole_sphere_gross_source_kg_s'],'Al_diagnostic_baseline1000_plus_excess_on_amplified_background_kg_s':next(z for z in gasrows if z['T_K']==1000)['whole_sphere_gross_source_kg_s']+sum(z['mean_excess_Al_source_kg_s'] for z in vaporrows.values()),'Al_aggregate_scope':'The baseline1000 plus excess on conservative amplified backgrounds is a diagnostic, not an exact operating mean or a rigorously complete time-envelope integral. The1000K baseline alone establishes the open-cavity vapor concern.','packet_gas_acceptance':budgets,'hydrogen_gas_return':hgas,'dust_acceptance':dust,'oxide_equilibrium':ox,'oxide_scope':'Seven ideal-gas species Al,AlO,AlO2,Al2O,Al2O2,O,O2; alpha-Al2O3 activity1; congruent effusive O/Al flux3/2 and unit evaporation coefficients. No hydrogen/carbon reduction, radiation chemistry, sputtering, dissociation kinetics or defect enhancement. This is a thermal equilibrium reference, not irradiated service prediction.','source_tables':{'XCOM_keV_cm2_g':{k:v.tolist() for k,v in TABLES.items()},'Al_vapor_T':AL_T.tolist(),'Al_vapor_logKf':AL_LOGK.tolist(),'oxide_T':GT.tolist(),'formation_G_kJ_mol':G},'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':all(checks)},'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',default=str(Path(__file__).resolve().parent.parent));p.add_argument('--output',default=str(Path(__file__).with_name('receiver-radiation.json')));a=p.parse_args();o=main(a.root);Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({'checks':o['checks'],'Al_base_kg_s':o['Al1000K_baseline_source_kg_s'],'Al_hybrid_diagnostic_kg_s':o['Al_diagnostic_baseline1000_plus_excess_on_amplified_background_kg_s'],'Al_pulses':o['Al_backed_surface_pulse_excess'],'oxide':[z for z in o['oxide_equilibrium'] if z['temperature_K'] in [1000,1400,1700,2000]],'flash':{s:d['Planck_shape_cases']['1.0'] for s,d in o['flash'].items()}},indent=2))
