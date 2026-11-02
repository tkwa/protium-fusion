"""Atomic threshold inventories and charge-conditional expansion timescales.

No charge-state network, no temperature trajectory, no new hydrodynamics.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--output',type=Path,default=Path(__file__).with_name('charge-state.json'));args=ap.parse_args();root=args.root
loadpath=root/'membrane-receiver/pulse-loading.json';cyclepath=root/'open-loop-redesign/station-inputs.json';opspath=root/'open-loop-redesign/operating-cases.json'
p=json.loads(loadpath.read_text());cycle=json.loads(cyclepath.read_text())['cycle'];ops=json.loads(opspath.read_text())['cases']['high_source_rating_surface']
MU=p['atomic_mass_unit_kg'];EV=1.602176634e-19;KB=1.380649e-23;Z=np.array(p['Z']);A=np.array(p['A']);R=p['radius_m'];area=4*math.pi*R*R
# Rounded source values. N/O are NIST table values reproduced in Horny thesis.
# Trace elements are bounded separately, rather than filled with unverified data.
ION={1:[13.5984],2:[24.5874,54.4178],6:[11.2603,24.3833,47.8878,64.4939,392.087,489.993],7:[14.53414,29.6013,47.44924,77.4735,97.8902,552.0718,667.046],8:[13.61806,35.11730,54.9355,77.41353,113.8990,138.1197,739.29,871.4101]}
# Verner & Ferland table1: recombination from bare ion to H-like ion.
RR={'H+':[7.982e-11,.7480,3.148,7.036e5],'He2+':[1.891e-10,.7524,9.370,2.774e6],'C6+':[6.556e-10,.7567,65.23,2.446e7],'N7+':[7.586e-10,.7563,90.15,3.338e7],'O8+':[8.616e-10,.7563,119.1,4.352e7]}
def rr(name,T):
 a,b,t0,t1=RR[name];x=math.sqrt(T/t0);y=math.sqrt(T/t1);return a/(x*(1+x)**(1-b)*(1+y)**(1+b))
def quant(x,w):
 idx=np.argsort(x);cum=np.cumsum(w[idx])/sum(w);return np.interp([.01,.5,.99],cum,np.array(x)[idx]).tolist()
first=np.array([ION[int(z)][0] if z in ION else 0 for z in Z]);full=np.array([sum(ION[int(z)]) if z in ION else 0 for z in Z]);trace=np.array([14*z**3 if z not in ION else 0 for z in Z]);temps=[10.,197.,1000.,1e4,1e5,1e6]
checks=[]
def check(name,v):
 checks.append({'name':name,'passed':bool(v)});assert v,name
stages={}
for name,s in p['stages'].items():
 c=s['cells'];m=np.array(c['mass_kg']);r=np.array(c['radius_bounds_m']);v=np.array(c['velocity_bounds_m_s']);Y=np.array(c['transport_Y_endpoint']);ui=np.array(c['specific_ion_internal_J_kg']);T=ui*MU/(1.5*KB*Y.sum(1));rho=m/(4*math.pi/3*np.diff(r**3));pack=s['spectral_packets'];idx=np.array(pack['cell_index']);t=np.array(pack['arrival_time_s'])-s['endpoint']['time_s'];rrin=r[idx]+v[idx]*t;rrout=r[idx+1]+v[idx+1]*t;rho_wall=m[idx]/(4*math.pi/3*(rrout**3-rrin**3));yw=np.array(pack['Y_at_wall']);ne=rho_wall/MU*(yw@Z)/1e6;ni=rho_wall/MU*yw[:,Z>0].sum(1)/1e6;te=(rrout**3-rrin**3)/(3*(rrout**2*v[idx+1]-rrin**2*v[idx]));mf=np.array(pack['mass_fluence_kg_m2']);vel=np.array(pack['velocity_m_s']);kw=.5*mf*vel**2
 counts=np.array(s['totals']['species_at_wall_m2']);knownN=counts[Z>0].sum();atomic1=counts@first*EV;atomicfull=counts@full*EV;traceupper=counts@trace*EV;K=s['totals']['kinetic_fluence_exact_nodal_J_m2'];rates={}
 for temp in temps:
  one={}
  for ion in RR:
   alpha=rr(ion,temp);chi=ne*alpha*te;one[ion]={'alpha_cm3_s':alpha,'n_alpha_texp_energy_quantiles_1_50_99':quant(chi,kw),'n_alpha_texp_max_all_packets':float(chi.max())}
  # NRL singly-ionized diagnostic uses one free electron per nucleus, not Z electrons.
  a3=8.75e-27*(KB*temp/EV)**(-4.5);ch3=ni*ni*a3*te
  ae=(3/(4*math.pi*ne*1e6))**(1/3);ai=(3/(4*math.pi*ni*1e6))**(1/3);gamma_e=1.43996448e-9/(ae*(KB*temp/EV));gamma_n7=49*1.43996448e-9/(ai*(KB*temp/EV))
  rates[str(temp)]={'co_temperature_coupling_diagnostics':{'electron_Gamma_max':float(gamma_e.max()),'bare_N7_Gamma_max':float(gamma_n7.max())},'radiative_bare_ion_capture':one,'three_body_singly_ionized_diagnostic':{'alpha3_cm6_s':a3,'n2_alpha3_texp_energy_quantiles_1_50_99':quant(ch3,kw),'maximum_all_packets':float(ch3.max())}}
 check(name+' ballisticdensitypositive',np.all(rho_wall>0));check(name+' masspacketclosure',abs(mf.sum()/s['totals']['mass_fluence_kg_m2']-1)<1e-12);check(name+' chargeupper',np.all(ne>=ni*(1-1e-12)));check(name+' knownspeciesdominance',counts[np.isin(Z,list(ION))].sum()/counts.sum()>.999)
 stages[name]={'endpoint_temperature_K_range':[float(T.min()),float(T.max())],'endpoint_density_kg_m3_range':[float(rho.min()),float(rho.max())],'wall_ballistic_density_kg_m3_energy_quantiles_1_50_99':quant(rho_wall,kw),'fully_stripped_electron_density_upper_cm3_energy_quantiles_1_50_99':quant(ne,kw),'fully_stripped_electron_density_upper_cm3_max':float(ne.max()),'expansion_time_s_energy_quantiles_1_50_99':quant(te,kw),'atomic_threshold_inventory':{'single_charge_known_J_m2':atomic1,'fully_strip_known_J_m2':atomicfull,'omitted_trace_fullstrip_upper_J_m2':traceupper,'single_charge_50pct_production_J_m2':2*(atomic1+traceupper),'fully_strip_50pct_production_J_m2':2*(atomicfull+traceupper),'single_charge_50pct_cost_over_kinetic':2*(atomic1+traceupper)/K,'fully_strip_50pct_cost_over_kinetic':2*(atomicfull+traceupper)/K,'known_elements_fraction_nuclei':counts[np.isin(Z,list(ION))].sum()/counts.sum()},'charge_fluence_C_m2':{'single_all_nuclei':knownN*EV,'fully_stripped':float(counts@Z*EV)},'temperature_conditional_rates':rates}
# Apply the repeating cycle's stage frequencies, retaining actual finite-flight daughters.
catmap={'c12':'B2','c13':'C13','n14':'A','n15':'B1'};averages={k:0. for k in ['single_charge_threshold_W_m2','fully_strip_threshold_W_m2','trace_upper_W_m2','single_charge_current_A_m2','fully_stripped_current_A_m2']}
for cat,s in cycle['stage_pulse_data'].items():
 rate=ops['helium_per_s']*s['pulses_per_He'];data=stages[catmap[cat]];a=data['atomic_threshold_inventory'];q=data['charge_fluence_C_m2'];averages['single_charge_threshold_W_m2']+=rate*a['single_charge_known_J_m2'];averages['fully_strip_threshold_W_m2']+=rate*a['fully_strip_known_J_m2'];averages['trace_upper_W_m2']+=rate*a['omitted_trace_fullstrip_upper_J_m2'];averages['single_charge_current_A_m2']+=rate*q['single_all_nuclei'];averages['fully_stripped_current_A_m2']+=rate*q['fully_stripped']
for target in ['single_charge','fully_strip']:
 a=averages[target+'_threshold_W_m2']+averages['trace_upper_W_m2'];averages[target+'_50pct_production_W_m2']=2*a;averages[target+'_50pct_share_of_outgoing_kinetic']=2*a/(cycle['bulk_kinetic_MeV_per_He']*ops['helium_per_s']*EV*1e6/area);averages[target+'_50pct_share_of_frozen_net']=2*a/(ops['net_W']/area)
# Planck tail bound: at x>>1, Bose factor differs negligibly from exp(-x).
x=ION[6][0]/(KB*1000/EV);tail=15/math.pi**4*math.exp(-x)*(x**3+3*x*x+6*x+6)
check('Plancktail_wien_valid',x>100)
for ion,z in [('H+',1),('C6+',6),('N7+',7),('O8+',8)]:
 te_ev=KB*1e4/EV;beta=13.5984*z*z/te_ev;seaton=5.2e-14*z*math.sqrt(beta)*(.43+.5*math.log(beta)+.469*beta**(-1/3))
 check('RR_vs_NRL_Seaton_15pct_'+ion,abs(rr(ion,1e4)/seaton-1)<.15)
out={'scope':'Charge-conditional rate diagnostics and neutral-atom threshold inventory; no collisional-radiative evolution or predicted wall charge fraction. Actual temperatures require a separate expansion/heating model.','receiver_radius_m':R,'ionization_energies_eV':ION,'rr_fit_coefficients':[{'ion':k,'a_cm3_s':v[0],'b':v[1],'T0_K':v[2],'T1_K':v[3]} for k,v in RR.items()],'formulas':{'ballistic_density':'m/[4pi/3*((r_outer+v_outer*t)^3-(r_inner+v_inner*t)^3)]','expansion_time':'rho/abs(d rho/dt)=V/(dV/dt)','radiative_test':'ne_upper*alpha_RR(Te)*t_exp','threebody_test':'n_ion^2*alpha3(Te)*t_exp, singly-ionized diagnostic only','atomic_energy':'sum(number_fluence_species*ground-state ionization threshold); not plasma thermodynamic free energy','50pct_production':'2*atomic threshold, no charge separation, electron collection, heating, beam transport or hardware power included'},'stages':stages,'selected_station_averages':averages,'blackbody_1000K':{'carbon_first_ionization_x':x,'energy_fraction_above_carbon_ground_threshold':tail,'hemisphere_energy_flux_above_threshold_W_m2':5.670374419e-8*1000**4*tail},'sources':[{'url':'https://physics.nist.gov/cuu/pdf/CODATA_RMP2000.pdf','use':'NIST static TableIII ground-state H/He/C successive and total binding energies; rounded atomic threshold arithmetic.'},{'url':'https://ft.nephy.chalmers.se/~vojtech/Theses/Doctoral_state_exam_Thesis_Vojtech_Horny.pdf','use':'Primary author-hosted plasma thesis,Table5.1,N/O successive thresholds reproduced fromNIST; rounded0.1eV totals sufficient.'},{'url':'https://www.pa.uky.edu/~verner/rec.html','use':'VernerFerland1996 author data and fit validity3–1e9K;bare-to-H-like RR subset.'},{'url':'https://www.pa.uky.edu/~verner/dima/rec//table1.dat','use':'Downloaded RR coefficients;electron countN=1 means recombinedH-likeproduct.'},{'url':'https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf','use':'Previously retained primaryNRL reprint,p54–56 time-dependent charge kinetics,hydrogenicRR,singly-ionizedthreebody,conditionsforLTE.'},{'url':'https://arxiv.org/abs/astro-ph/0604144','use':'Badnell2006:level/metastable-resolvedRR and complementaryDR rates needed for a fullnetwork.'}], 'retained_atomic_source_sha256':{'verner-table1.dat':'2382f6896fe58ee1da58c9bd563fd7d02c0d94410dbd4150767eed26a54f5091','NRL_Plasma_Formulary_2023.pdf':'1a13dc8f0a1e493806adde235b22ccdd80eed939936a4622b36d0b563bea2596'},'access_notes':['NISTASD CGI request returned403; no further attempt was made to that endpoint. Static NISTpaper and author-hosted NIST-reproduced table provide threshold values.','OfficialNRL publicindex returned403; the already retained authoritativereprint was used.'],'source_hashes':{str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [loadpath,cyclepath,opspath]},'checks':{'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks}}
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'averages':averages,'stages':{k:{'inventory':v['atomic_threshold_inventory'],'ne':v['fully_stripped_electron_density_upper_cm3_energy_quantiles_1_50_99'],'RR_O197_max':v['temperature_conditional_rates']['197.0']['radiative_bare_ion_capture']['O8+']['n_alpha_texp_max_all_packets'],'3body197_max':v['temperature_conditional_rates']['197.0']['three_body_singly_ionized_diagnostic']['maximum_all_packets']} for k,v in stages.items()},'checks':out['checks']['passed']},indent=2))
