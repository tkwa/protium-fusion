#!/usr/bin/env python3
"""C13 receiver local-velocity-moment/temperature requirements, no hydro or PIC.
Uses the frozen shared endpoint packet. Run python3 debris_temperature.py.
Temperatures are conditional second moments; a Maxwellian is not guaranteed.
"""
import json,hashlib,math
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent;SRC=P.parent/'membrane-receiver/pulse-loading.json'
KB=1.380649e-23;EV=1.602176634e-19;MU=1.66053906660e-27;ME=9.1093837139e-31;C=299792458.;EPS=8.8541878128e-12
SOURCES={
'beta_Q_convention':'https://www.nndc.bnl.gov/nudat3/help/betadecay.jsp',
'O15_endpoint_primary_review':'https://lss.fnal.gov/archive/2024/pub/fermilab-pub-24-0489-sqms-v.pdf',
'nuclear_levels_N15':'https://nucldata.tunl.duke.edu/nucldata/HTML/A=15/15_04_1991.pdf',
'collision_formulary':'https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Plasma_Formulary_2023.pdf',
'particle_stopping':'https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf',
'lepton_escape_and_thermalization':'https://arxiv.org/abs/1605.07218',
'collective_self_confinement':'https://academic.oup.com/mnras/article/541/3/2197/8203299'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def quant(x,w,ps=(.01,.1,.5,.9,.99)):
    order=np.argsort(x);xx=np.array(x)[order];ww=np.array(w)[order];cdf=np.cumsum(ww)/sum(ww)
    return {str(p):float(np.interp(p,cdf,xx)) for p in ps}
def i2(x):
    x=np.asarray(x,dtype=float);small=x<.01;out=np.empty_like(x);out[small]=x[small]/3-x[small]**2/4+x[small]**3/10-x[small]**4/36+x[small]**5/168
    z=x[~small];out[~small]=(2-np.exp(-z)*(z*z+2*z+2))/(z*z)
    return out

def beta(q=1.732,A=15,Z=7,a=0.,n=24000):
    # Allowed spectrum, point-charge Fermi factor for positrons; no radiative/finite-size correction.
    t=np.linspace(1e-9,q-1e-9,n);me=.510998950;p=np.sqrt(t*(t+2*me));b=p/(t+me);nu=q-t;farg=2*math.pi/137.035999084*Z/b;fermi=farg/np.expm1(np.minimum(farg,700));w=p*(t+me)*nu**2*fermi
    norm=np.trapz(w,t);mean=lambda y:float(np.trapz(w*y,t)/norm)
    p2=p*p+nu*nu+2*a*b*p*nu/3
    return {'beta_endpoint_MeV':q,'daughter_A':A,'daughter_Z':Z,'angular_correlation_coefficient_scenario':a,'mean_positron_kinetic_MeV':mean(t),'mean_neutrino_MeV':mean(nu),'mean_daughter_recoil_eV':mean(p2)/(2*A*931.49410242)*1e6,'positron_momentum_endpoint_MeV_c':math.sqrt(q*(q+2*me)),'daughter_recoil_max_eV':q*(q+2*me)/(2*A*931.49410242)*1e6,'scope':'Allowed beta phase-space approximation with repulsive point Fermi factor. Angular coefficient is a sensitivity, not an evaluated O15 mixing-ratio fit; endpoint and lifetime are evaluated-source inputs.'}

def main():
    d=json.loads(SRC.read_text());s=d['stages']['C13'];cells=s['cells'];pa=s['spectral_packets'];sp=d['species'];aa=np.array(d['A']);zz=np.array(d['Z']);m=np.array(cells['mass_kg']);r=np.array(cells['radius_bounds_m']);v=np.array(cells['velocity_bounds_m_s']);y=np.array(cells['transport_Y_endpoint']);ci=np.array(pa['cell_index']);vp=np.array(pa['velocity_m_s']);mp=np.array(pa['mass_fluence_kg_m2']);t=np.array(pa['arrival_time_s'])-s['endpoint']['time_s'];gamma=1/np.sqrt(1-(vp/C)**2);yi=y.sum(1);ui=np.array(cells['specific_ion_internal_J_kg']);T0=ui*MU/(1.5*KB*yi);rho0=m/(4*np.pi/3*np.diff(r**3));rw0=r[ci]+v[ci]*t;rw1=r[ci+1]+v[ci+1]*t;rho=m[ci]/(4*np.pi/3*(rw1**3-rw0**3));texp=(rw1**3-rw0**3)/(3*(rw1**2*v[ci+1]-rw0**2*v[ci]));Ywall=np.array(pa['Y_at_wall']);ni=rho/MU*Ywall.sum(1);ne=rho/MU*(Ywall@zz);number=mp/MU*Ywall.sum(1);K=.5*mp*vp**2;N=sum(number);KF=sum(K)
    tg=T0[ci]*(rho/rho0[ci])**(2/3);abar=(Ywall@aa)/Ywall.sum(1);tgeo=abar*MU*r[-1]**2/(3*KB*t*t)
    par=sp.index('o15');o15=mp/MU*y[ci,par];lam=next(x['rate_s_inverse'] for x in d['weak_decay_channels'] if x['parent']=='o15');x=lam*t/gamma;dec=-np.expm1(-x);weight=i2(x)
    betas={str(a):beta(a=a) for a in [-1/3,0.,1.]};b=betas['0.0'];recoil=o15*weight*b['mean_daughter_recoil_eV']*EV;recoil_raw=o15*dec*b['mean_daughter_recoil_eV']*EV;tr=2*recoil/(3*KB*number)
    betaJ=o15*weight*b['mean_positron_kinetic_MeV']*1e6*EV;betaJ_raw=o15*dec*b['mean_positron_kinetic_MeV']*1e6*EV
    # Thermal deposition scenarios use the same t_decay^2/t_arrival^2 adiabatic dilution.
    # Fully stripped electron count is only a charge-conditional denominator.
    electrons=mp/MU*(Ywall@zz);branches=[]
    for f in [0.,1e-6,.001,.01,.1,1.]:
        heat=f*betaJ
        added_bulk=sum(recoil_raw-recoil+f*(betaJ_raw-betaJ))
        branches.append({'instantaneous_positron_energy_deposition_fraction':f,'late_thermal_energy_J_m2':float(sum(heat)),'electron_only_temperature_equivalent_full_strip_number_mean_K':float(2*sum(heat)/(3*KB*sum(electrons))),'common_particle_temperature_if_equilibrated_K':float(2*sum(heat+recoil)/(3*KB*sum(number+electrons))),'random_energy_fraction_of_bulk_plus_random':float(sum(heat+recoil)/(KF+sum(heat+recoil))),'bulk_fraction_of_bulk_plus_random':float(KF/(KF+sum(heat+recoil))),'additional_coherent_energy_in_ideal_closed_branch_J_m2':float(added_bulk),'energy_closed_bulk_fraction':float((KF+added_bulk)/(KF+sum(recoil_raw)+f*sum(betaJ_raw))),'energy_closed_total_increment_J_m2':float(sum(recoil_raw)+f*sum(betaJ_raw)),'scope':'Conditional instantaneous local beta kinetic-energy deposition followed by nonrelativistic adiabatic dilution. Fixed-KF fractions are diagnostic ratios against the archived trajectory. Energy-closed branch assigns raw-minus-filtered recoil/beta energy to additional coherent expansion work, with no new radiation loss. This is a global ideal ledger, not a recalculated parcel trajectory. No annihilation photons or neutrinos added. Actual confinement/delayed thermalization changes the kernel; no retention probability assigned.'})
    coll=[]
    for Ti in [100.,1000.,3000.,10000.,100000.]:
        tev=Ti*KB/EV
        for z in [1,3,7]:
            freq=4.8e-8*z**4*ni/1e6*10/(math.sqrt(14)*tev**1.5)
            # Equal-charge nitrogen proxy, not a full multispecies Landau collision operator.
            ai=(3/(4*np.pi*ni))**(1/3);coupling=z*z*EV*EV/(4*np.pi*EPS*ai*KB*Ti)
            coll.append({'trial_Ti_K':Ti,'fixed_ion_charge_scenario':z,'Coulomb_logarithm_assumed':10,'nu_ii_times_local_expansion_time_energy_quantiles':quant(freq*texp,K),'Gamma_ii_energy_quantiles':quant(coupling,K),'scope':'NRL weak-coupling like-ion relaxation estimate, N14 mass proxy, n_i from exact shell expansion; charge is assumed, not predicted. Invalid when strong coupling or neutral fraction dominates.'})
    equil=[]
    for temp in [400.,1000.,3000.,10000.,100000.]:
        tev=temp*KB/EV
        for z in [1,3,7]:
            nue=3.2e-9*z*z*(ni/1e6)*10/(14*tev**1.5)
            nui=z*nue
            equil.append({'trial_common_temperature_K':temp,'nitrogen_proxy_ion_charge':z,'Coulomb_logarithm_assumed':10,'electron_relaxation_times_local_expansion_time_energy_quantiles':quant(nue*texp,K),'ion_relaxation_times_local_expansion_time_energy_quantiles':quant(nui*texp,K),'temperature_difference_decay_time_s_energy_quantiles':quant(1/(nue+nui),K),'scope':'NRL near-equal-temperature two-species energy exchange. N14 proxy with n_e=Z*n_i; not actual charge prediction or integration from unequal temperatures. High charge can share recoil energy with electrons within flight time; warm beta-heated electrons can couple much more slowly.'})
    mapping=[]
    for Ti,Te in [(3000.,100.),(3000.,10000.),(10000.,10000.),(1e5,1e5),(1e7,1e7)]:
        thermal=1.5*KB*(N*Ti+sum(electrons)*Te)
        mapping.append({'ion_temperature_K':Ti,'electron_temperature_K':Te,'electrons_full_strip_scenario':True,'thermal_J_m2':float(thermal),'bulk_percent_of_bulk_plus_thermal':float(100*KF/(KF+thermal))})
    totals=s['totals'];energy={k:float(m@np.array(cells[n])/s['area_m2']) for k,n in [('ions','specific_ion_internal_J_kg'),('electrons_EOS','specific_electron_internal_J_kg'),('radiation','specific_radiation_internal_J_kg')]}
    beta_kin_mean=b['mean_positron_kinetic_MeV'];column=totals['mass_fluence_kg_m2']
    direct=[{'mass_stopping_power_MeV_cm2_g_scenario':st,'single_column_energy_lost_MeV':st*column*.1,'fraction_of_mean_positron_kinetic':st*column*.1/beta_kin_mean,'scope':'Order-of-magnitude direct-crossing stopping diagnostic;1kg/m²=.1g/cm². Column is one receiver-normal pulse fluence, not a complete angular/time-dependent path. Magnetic/collective confinement and early dense deposition excluded.'} for st in [1.,3.,10.]]
    # Independent numerical integration verifies recoil/filter kernel including small-x limit.
    checks=[]
    for a in [1e-6,.001,.1,1.,10.,100.]:
        z=np.linspace(0,1,20001);num=np.trapz(a*np.exp(-a*z)*z*z,z);checks.append(abs(float(i2(np.array([a]))[0])-num)<1e-8)
    checks +=[abs(sum(number)/totals['nuclei_fluence_m2']-1)<1e-10,abs(KF/totals['kinetic_fluence_smooth_J_m2']-1)<1e-10,abs(sum(energy.values())/totals['residual_internal_fluence_J_m2']-1)<1e-8,all(weight<=dec),all(weight>=0)]
    checks += [abs((row['additional_coherent_energy_in_ideal_closed_branch_J_m2']+row['late_thermal_energy_J_m2']+sum(recoil))/row['energy_closed_total_increment_J_m2']-1)<1e-12 for row in branches]
    out={'scope':'C13 local random-motion and conditional temperature screen at receiver. Initial radial velocity spread and successive arrival speeds are coherent flow; only peculiar velocity relative to local mean contributes to temperature. No new hydro, collisional-radiative, kinetic-electron or magnetic evolution.','source_sha256':sha(SRC),'code_sha256':sha(__file__),'primary_sources':SOURCES,'constants':{'kB_J_K':KB,'eV_J':EV,'atomic_mass_kg':MU},'radius_m':d['radius_m'],'endpoint':{'time_s':s['endpoint']['time_s'],'outer_radius_m':r[-1],'gas_temperature_min_max_K':[float(min(T0)),float(max(T0))],'gas_density_min_max_kg_m3':[float(min(rho0)),float(max(rho0))],'energy_fluences_J_m2':energy,'bulk_kinetic_J_m2':KF,'internal_to_bulk_ratio':sum(energy.values())/KF},'wall':{'energy_arrival_1_50_99_s':s['smooth']['arrival_1_50_99_percent_energy_s'],'ion_density_energy_quantiles_m3':quant(ni,K),'electron_density_full_strip_energy_quantiles_m3':quant(ne,K),'expansion_time_energy_quantiles_s':quant(texp,K),'gas_only_adiabatic_temperature_energy_quantiles_K':quant(tg,K),'force_free_from_endpoint_geometric_temperature_bound_energy_quantiles_K':quant(tgeo,K),'electron_per_nucleus_full_strip':float(sum(electrons)/N),'kinetic_MeV_per_nucleus':float(KF/N/EV/1e6),'initial_number_fractions':{sp[j]:float(sum(mp/MU*y[ci,j])/N) for j in range(len(sp)) if sum(mp/MU*y[ci,j])/N>1e-7}},'recoil':{'dominant_parent':'o15','lifetime_half_s':math.log(2)/lam,'beta_spectrum_scenarios':betas,'decayed_O15_per_all_nuclei':float(sum(o15*dec)/N),'recoil_without_spatial_filter_J_m2':float(sum(recoil_raw)),'recoil_with_decaytime_squared_filter_J_m2':float(sum(recoil)),'filtered_ion_temperature_number_mean_K':float(2*sum(recoil)/(3*KB*N)),'filtered_ion_temperature_energy_quantiles_K':quant(tr,K),'unfiltered_ion_temperature_number_mean_K':float(2*sum(recoil_raw)/(3*KB*N)),'bulk_fraction_with_filtered_recoil':float(KF/(KF+sum(recoil))),'O15_positron_kinetic_energy_unfiltered_J_m2':float(sum(betaJ_raw)),'O15_positron_kinetic_energy_filtered_J_m2':float(sum(betaJ)),'kernel':'For x=lambda*t/gamma,I2=[2-exp(-x)*(x²+2x+2)]/x² = integral0^t lambda/gamma*exp(-lambda*tau/gamma)*(tau/t)^2dtau. A kick at tau becomes local peculiar velocity approximately(tau/t)*delta_v in homologous expansion. Finite source size/early acceleration and perturbed arrival sampling remain unresolved.','scope':'O15 dominates. O14,N13 initial abundances <3e-8 each, so their gamma/beta recoils are negligible here; other stages differ. Effective Ti can describe a cold core plus minority fast daughters instead of Maxwellian equilibrium.'},'beta_retention_branches':branches,'direct_crossing_positron_stopping_scenarios':direct,'fixed_charge_collision_checks':coll,'electron_ion_equal_temperature_relaxation_checks':equil,'temperature_to_bulk_fraction_mapping':mapping,'energy_upper_scope':{'all_inflight_decay_total_Q_including_unrecoverable_neutrino_J_m2':totals['inflight_decay_total_Q_upper_J_m2'],'all_decay_Q_over_bulk':totals['inflight_decay_total_Q_upper_J_m2']/KF,'note':'This total mass-difference Q is an intentionally loose ceiling, not heat deposited in the plume.'},'judgment':{'conditional_ion_temperature_central_K':3000,'conditional_subjective80pct_interval_K':[1000,10000],'conditions':'Pulse-number-averaged recoil-equivalent ion second moment, or characteristic O15-bearing stream; monotone expansion, no strong shocks/external braking, no significant beta energy transferred to ions beyond daughter recoil, and negligible net recoil-energy transfer to electrons/radiation. Interval is subjective judgment around the computed recoil scale, including momentum correlations, phase-space filtering and unresolved early thermal tail. It is neither a spatial percentile band nor a Monte Carlo confidence interval. Cold O15-poor outer parcels fall far below it.','fully_stripped_recoil_shared_common_temperature_K':float(2*sum(recoil)/(3*KB*sum(number+electrons))),'unconditional_ion_temperature':'No preferred unconditional value: sharing with electrons, radiative losses, and beta deposition are not solved. Full-strip equilibration can give about400K instead of the ion-retained3000K scale.','electron_temperature':'No defensible unconditional central or80%interval from saved state; separate beta deposition/charge/thermalization branches supplied. A common Ti=Te is an extra assumption.'},'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':all(checks)}}
    (P/'debris-temperature.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'wall':out['wall'],'recoil':out['recoil'],'beta_branches':branches},indent=2))
if __name__=='__main__':main()
