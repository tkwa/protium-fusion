"""Postprocess saved opaque-pulse trajectories; no hydrodynamic evolution."""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from microphysics import *

ap=argparse.ArgumentParser()
ap.add_argument('--repo',default=str(Path(__file__).resolve().parents[2]))
ap.add_argument('--trace-dir',default='/tmp/protium-research')
ap.add_argument('--out',default=str(Path(__file__).parent))
ap.add_argument('--stems',nargs='+')
a=ap.parse_args();ROOT=Path(a.repo);OUT=Path(a.out);OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT))
import pulse_model as pm
from fast_eos import FastEOS
from pynucastro.screening import PlasmaState,make_screen_factors,chugunov_2009
from numba import njit

@njit
def screening_all(T,rho,Y,Z,factor):
    out=np.empty(len(T))
    for i in range(len(T)):
        out[i]=chugunov_2009(PlasmaState(T[i],rho[i]/1000,Y[i],Z),factor)
    return out

stems=['compression-A-128-full-burning','compression-B1-128-burning',
       'compression-B2-128-burning','compression-C13-64-burning','compression-He-64-r30-burning']
if a.stems:stems=a.stems
eos=FastEOS();results=[]
for stem in stems:
    meta=ROOT/(stem+'.json');trace=Path(a.trace_dir)/(stem+'.npz')
    d=json.loads(meta.read_text());raw=np.load(trace)
    nz=d['zones'];nv=int(raw['variables_per_zone']);nn=nv-5
    z=raw['state'].reshape(nz,nv,-1);nt=z.shape[-1];t=raw['time_s']
    m=d['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3)
    r=np.vstack([np.zeros(nt),z[:,nn+1]*d['R0_m']]);dr=np.diff(r,axis=0)
    rc=.5*(r[1:]+r[:-1]);rho=m[:,None]/(4*np.pi/3*np.diff(r**3,axis=0))
    T=z[:,nn]*d['T0_K'];Y=z[:,:nn].transpose(0,2,1)
    assert raw['species'].tolist()==pm.SPECIES
    Y=np.maximum(Y,0);Y=Y/np.sum(Y*pm.AA,axis=-1)[...,None]
    ye=np.sum(Y*pm.ZZ,axis=-1);yi=np.sum(Y,axis=-1)
    ne=rho/MU*ye;z2ni=rho/MU*np.sum(Y*pm.ZZ**2,axis=-1)
    state=eos.state(rho/1000,T,ye,yi)
    cv=state.c_v*1e-4;u=state.e*1e-4;p=state.p*.1
    cs=np.sqrt(state.gamma1*p/rho)
    nuclear=-m[:,None]*np.sum((Y[:,1:]-Y[:,:-1])*pm.MASS,axis=-1)*pm.MEV/pm.MU
    inc=np.maximum(nuclear,0)
    w=np.zeros_like(T);w[:,:-1]+=.5*inc;w[:,1:]+=.5*inc
    wt=(m[:,None]/d['mass_kg'])*np.ones_like(T)
    eta,nminus,nplus=fermions(T,ne)
    kappa_s,valid_sc=poutanen(T,eta,nminus+nplus,rho)
    # Outside its published fit domain retain only an explicitly named
    # diagnostic Thomson/phase-space estimate, not a qualified opacity.
    xf=HB*(3*np.pi*np.pi*ne)**(1/3)/(ME*C)
    ef=ME*C*C*(np.sqrt(1+xf*xf)-1)
    fallback=SIGT*(nminus+nplus)/rho*np.minimum(1,1.5*KB*T/ef)
    kappa_diag=np.where(valid_sc,kappa_s,fallback)
    alpha_ff,emission,alpha_planck=freefree(T,ne,z2ni,x=5.)
    alpha_s=kappa_diag*rho
    therm_time=1/(C*alpha_ff)
    tau_eff=np.sqrt(3*alpha_ff*(alpha_ff+alpha_s))*dr
    wp=np.sqrt(ne*E*E/(EPS*ME))
    plasma_over_kt=HB*wp/(KB*T)
    # If this very severe common suppression still yields a tiny thermal
    # relaxation time/large cell optical depth, coefficient refinement has
    # substantial margin. It is a stress test, not a certified lower bound.
    suppressed=1e-4
    tvar=np.minimum(d['elapsed_s'],np.minimum(
      T/np.maximum(np.abs(np.gradient(T,t,axis=1)),1e-99),
      rho/np.maximum(np.abs(np.gradient(rho,t,axis=1)),1e-99)))
    ui,pi,gdominant,ge=coulomb(T,rho,Y,pm.ZZ)
    # Classical electron/ion temperature equilibration scale, summed over ions.
    # NRL equal-temperature Maxwellian exchange rate in SI via the standard
    # Spitzer expression. lnLambda=1 deliberately slows exchange; degeneracy
    # and mildly relativistic corrections remain explicit limitations.
    invtei=np.zeros_like(T)
    for j,(massnumber,charge) in enumerate(zip(pm.AA,pm.ZZ)):
        if charge==0:continue
        ni=rho/MU*Y[...,j];mi=massnumber*MU
        invtei+=8*np.sqrt(2*np.pi)*ni*charge**2*EC**2/(3*ME*mi)*(KB*T/ME+KB*T/mi)**(-1.5)
    tei=1/invtei
    kcond,kclass,kdeg=conductivity_envelope(T,ne,z2ni)
    chi=kcond/(rho*cv); lcond=np.sqrt(chi*tvar)
    facek=np.maximum(kcond[:-1],kcond[1:])
    gradT=np.diff(T,axis=0)/np.diff(rc,axis=0)
    luminosity=-4*np.pi*r[1:-1]**2*facek*gradT
    transported=np.trapezoid(np.abs(luminosity),t,axis=1)
    # Summing interfaces intentionally overcounts energy transported repeatedly.
    # It is not the radiated/lost energy and assumes the frozen trajectories.
    faceenergy=float(transported.sum())
    boundary_lum=np.vstack([np.zeros(nt),luminosity,np.zeros(nt)])
    conductive_cell_power=-np.diff(boundary_lum,axis=0)
    redistributed=float(np.trapezoid(.5*np.abs(conductive_cell_power).sum(axis=0),t))
    # Critical nuclear reactants, including both sequential alpha captures
    # for the triple-alpha factory. No automatic reverse-rate adjustment here.
    pairs={'A':[('p','n14')],'B1':[('p','n15')],
           'B2':[('p','c12'),('p','n13')],'C13':[('p','c13')],
           'He':[('he4','he4'),('he4','be8')]}[d['stage']]
    screens={}
    for first,second in pairs:
        sf=make_screen_factors(first,second)
        logs=screening_all(T.ravel(),rho.ravel(),Y.reshape(-1,nn),pm.ZZ,sf).reshape(T.shape)
        screens[first+'+'+second]={
            'log_enhancement_burn_weighted_q01_q50_q99':weighted_quantile(logs,w),
            'enhancement_burn_weighted_q01_q50_q99':weighted_quantile(np.exp(np.minimum(logs,700)),w),
            'initial_enhancement_mass_q01_q50_q99':weighted_quantile(np.exp(np.minimum(logs[:,0],700)),m),
            'max_log_enhancement':float(logs.max())}
    active=w>0
    q=lambda values:weighted_quantile(values,w)
    record={
       'stem':stem,'stage':d['stage'],'trace_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),
       'metadata_sha256':hashlib.sha256(meta.read_bytes()).hexdigest(),
       'zones':nz,'samples':nt,'duration_s':d['elapsed_s'],
       'input_energy_J':d['initial_energy_J'],'net_nuclear_J':d['nuclear_energy_J'],
       'positive_sampled_nuclear_increments_J':float(inc.sum()),
       'T_K_range':[float(T.min()),float(T.max())],
       'rho_kg_m3_range':[float(rho.min()),float(rho.max())],
       'eta_range':[float(eta.min()),float(eta.max())],
       'burn_weighted_q01_q50_q99':{
          'T_K':q(T),'rho_kg_m3':q(rho),'eta':q(eta),
          'Poutanen_kappa_m2_kg_within_domain':q(kappa_s),
          'freefree_alpha_x5_per_m':q(alpha_ff),
          'freefree_Planck_kappa_m2_kg':q(alpha_planck/rho),
          'freefree_x5_relaxation_s':q(therm_time),
          'effective_optical_depth_per_cell_x5':q(tau_eff),
          'relaxation_over_local_variation_time_suppressed_1e4':q(therm_time/tvar/suppressed),
          'effective_optical_depth_per_cell_suppressed_1e4':q(tau_eff*suppressed),
          'pair_fraction_of_lepton_density':q(2*nplus/(nminus+nplus)),
          'equilibrium_pair_rest_energy_fraction':q(2*nplus*ME*C*C/(rho*u)),
          'classical_electron_ion_equilibration_s_lnLambda1':q(tei),
          'classical_electron_ion_equilibration_over_variation_time':q(tei/tvar),
          'plasma_frequency_energy_over_kT':q(plasma_over_kt),
          'electron_fermi_energy_over_kT':q(ef/(KB*T)),
          'dominant_ion_Gamma':q(gdominant),
          'ion_ion_pressure_fraction':q(-pi/p),
          'ion_ion_internal_energy_fraction':q(-ui/u),
          'radiation_internal_energy_fraction':q(ARAD*T**4/(rho*u)),
          'conductivity_envelope_W_m_K':q(kcond),
          'conductive_diffusion_length_m':q(lcond),
          'conductive_diffusion_length_over_cell_width':q(lcond/dr),
          'conductive_diffusion_length_over_temperature_scale':q(lcond*np.abs(np.gradient(T,axis=0)/np.gradient(rc,axis=0))/T)},
       'Poutanen_outside_domain_burn_weight_fraction':float(w[~valid_sc].sum()/w.sum()),
       'initial_ion_ion_energy_fraction_input':float(-m@ui[:,0]/d['initial_energy_J']),
       'initial_ion_ion_pressure_fraction_mass_q01_q50_q99':weighted_quantile(-pi[:,0]/p[:,0],m),
       'initial_ion_ion_Gamma_mass_q01_q50_q99':weighted_quantile(gdominant[:,0],m),
       'conductive_sum_absolute_interface_transfer_over_nuclear':faceenergy/d['nuclear_energy_J'],
       'conductive_redistributed_energy_over_nuclear':redistributed/d['nuclear_energy_J'],
       'conductive_max_interface_transfer_over_nuclear':float(transported.max()/d['nuclear_energy_J']),
       'screening':screens,
       'diagnostic_flags':{
         'freefree_nonrelativistic_mildly_relativistic_Tmax_MeV':float(KB*T.max()/E/1e6),
         'freefree_electron_degeneracy_or_plasma_dispersion_requires_correction':True,
         'conductivity_100x_sensitivity_allowed_but_not_rigorous_upper_bound':True,
         'Coulomb_estimate_ion_ion_only_electron_polarization_and_mixing_omitted':True,
         'no_coupled_hydro_or_netgain_prediction':True}}
    results.append(record)
    np.savez_compressed(OUT/(stem+'-diagnostics.npz'),time_s=t,T_K=T,rho_kg_m3=rho,
      mass_kg=m,radius_m=r,normalized_molar=Y,nuclear_weight_J=w,eta=eta,
      kappa_scatter=kappa_s,alpha_ff_x5_per_m=alpha_ff,tau_eff_x5=tau_eff,
      radiation_fraction=ARAD*T**4/(rho*u),ion_ion_pressure_fraction=-pi/p,
      conductivity_envelope=kcond,conductive_length_m=lcond)
    print(stem,json.dumps({k:record[k] for k in ['Poutanen_outside_domain_burn_weight_fraction',
       'initial_ion_ion_energy_fraction_input','conductive_sum_absolute_interface_transfer_over_nuclear']}),flush=True)
out={'scope':'Frozen-trajectory microphysical diagnostics; no altered pulse yield, gain or plant closure claimed.',
     'quantile_weights':'Positive nuclear mass-defect changes between saved checkpoints, half assigned to each endpoint. Not a reaction-resolved flow integral.',
     'rows':results}
(OUT/'pulse-transport-audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
