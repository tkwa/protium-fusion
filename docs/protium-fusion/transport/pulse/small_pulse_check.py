"""Independent microphysics extraction for the driver's scaled A trajectory.

Uses the already verified diagnostic coefficients; it does not evolve fluid,
calculate merger losses, or predict the selected fictional measured gain.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from microphysics import *
ap=argparse.ArgumentParser()
ap.add_argument('--repo',default=str(Path(__file__).resolve().parents[2]))
ap.add_argument('--out',default=str(Path(__file__).parent))
a=ap.parse_args();ROOT=Path(a.repo);OUT=Path(a.out)
sys.path.insert(0,str(ROOT))
import pulse_model as pm
p=ROOT/'launcher-reference/driver-physics/scaled-hot-A-64-R30'
d=json.loads(p.with_suffix('.json').read_text());raw=np.load(p.with_suffix('.npz'))
nz=d['zones'];nv=int(raw['variables_per_zone']);nn=nv-5
z=raw['state'].reshape(nz,nv,-1);nt=z.shape[-1];t=raw['time_s']
ind=[pm.SPECIES.index(x) for x in raw['species'].tolist()]
A=pm.AA[ind];Z=pm.ZZ[ind]
m=d['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3)
r=np.vstack([np.zeros(nt),z[:,nn+1]*d['R0_m']]);dr=np.diff(r,axis=0)
rho=m[:,None]/(4*np.pi/3*np.diff(r**3,axis=0))
T=z[:,nn]*d['T0_K'];Y=z[:,:nn].transpose(0,2,1)
Y=np.maximum(Y,0);Y/=np.sum(Y*A,axis=-1)[...,None]
inc=np.maximum(-m[:,None]*np.sum(np.diff(Y,axis=1)*pm.MASS[ind],axis=-1)*pm.MEV/pm.MU,0)
w=np.zeros_like(T);w[:,:-1]+=.5*inc;w[:,1:]+=.5*inc
ne=rho/MU*np.sum(Y*Z,axis=-1);z2=rho/MU*np.sum(Y*Z**2,axis=-1)
eta,nm,np_=fermions(T,ne);ks,valid=poutanen(T,eta,nm+np_,rho)
af,em,aplanck=freefree(T,ne,z2)
tau=np.sqrt(3*af*(af+ks*rho))*dr
var=np.minimum(d['elapsed_s'],np.minimum(
 T/np.maximum(abs(np.gradient(T,t,axis=1)),1e-99),
 rho/np.maximum(abs(np.gradient(rho,t,axis=1)),1e-99)))
q=lambda x:weighted_quantile(x,w)
QMeV=pm.MASS[ind][0]+pm.MASS[ind][1]-pm.MASS[ind][2]
x=QMeV/(ME*C*C/pm.MEV)
# Total Klein-Nishina cross section on a stationary free electron.
kn=.75*((1+x)/x**3*(2*x*(1+x)/(1+2*x)-np.log(1+2*x))
       +np.log(1+2*x)/(2*x)-(1+3*x)/(1+2*x)**2)
required={'selected_fictional_measured_capture_stage_gain':2.6,
 'qualification':'Chosen future experimental outcome; no computed 3D loss or engineering forecast.',
 'delivered_input_J':d['initial_energy_J'],
 'implied_nuclear_energy_J':2.6*d['initial_energy_J'],
 'input_plus_nuclear_energy_J':3.6*d['initial_energy_J'],
 'implied_capture_conversion_at_same_input':2.6*d['initial_energy_J']/(d['nuclear_energy_J']/d['desired_product_per_initial_catalyst']),
 'source_electric_energy_if_driver_efficiency_0p3_J':d['initial_energy_J']/.3,
 'returned_electric_over_source_if_full_collection_and_heat_efficiency_0p4':.3*.4*3.6,
 'isotropic_equivalent_total_energy_fluence_at_1e9m_J_m2':3.6*d['initial_energy_J']/(4*np.pi*1e18),
 'isotropic_equivalent_total_energy_fluence_at_1e10m_J_m2':3.6*d['initial_energy_J']/(4*np.pi*1e20)}
checks={'all_cells_positive_density':bool(np.all(rho>0)),
 'ordered_radii':bool(np.all(dr>0)),'ordered_saved_times':bool(np.all(np.diff(t)>0)),
 'nuclear_weight_sum_relative_error':abs(float(w.sum()/inc.sum()-1)),
 'trace_end_equals_metadata_duration_relative_error':abs(float(t[-1]/d['elapsed_s']-1))}
assert checks['all_cells_positive_density'] and checks['ordered_radii'] and checks['ordered_saved_times']
assert checks['nuclear_weight_sum_relative_error']<1e-12 and checks['trace_end_equals_metadata_duration_relative_error']<1e-12
record={'scope':'Frozen actual scaled prescribed-hot A trajectory, reduced three-species network. No cold-material, driver, merger or multidimensional validation.',
 'trace_sha256':hashlib.sha256(p.with_suffix('.npz').read_bytes()).hexdigest(),
 'metadata_sha256':hashlib.sha256(p.with_suffix('.json').read_bytes()).hexdigest(),
 'diagnostic_microphysics_sha256':hashlib.sha256(Path(__file__).with_name('microphysics.py').read_bytes()).hexdigest(),
 'extractor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'zones':nz,'saved_samples':nt,'mass_kg':d['mass_kg'],'reference_radius_m':d['R0_m'],
 'initial_outer_radius_m':d['initial_surface_radius_m'],
 'minimum_sampled_outer_radius_m':float(r[-1].min()),'duration_s':d['elapsed_s'],
 'computed_idealized_capture_gain':d['deposited_nuclear_over_initial_energy'],
 'computed_capture_conversion':d['desired_product_per_initial_catalyst'],
 'computed_peak_temperature_K':d['peak_temperature_K'],
 'sampled_peak_density_kg_m3':float(rho.max()),
 'energy_residual_fraction_initial':d['energy_residual_fraction_initial'],
 'fit_outside_domain_nuclear_weight_fraction':float(w[~valid].sum()/w.sum()),
 'quantile_weight':'Positive nuclear mass-defect increments, half assigned to each sampled endpoint.',
 'burn_weighted_q01_q50_q99':{'rho_kg_m3':q(rho),'T_K':q(T),'Compton_opacity_m2_kg':q(ks),
   'effective_thermalization_depth_per_cell_x5':q(tau),
   'effective_thermalization_depth_both_coefficients_divided100':q(tau*.01),
   'absorption_time_over_local_variation_both_coefficients_divided100':q(100/(C*af*var))},
 'prompt_capture_gamma_diagnostic':{'capture_Q_MeV':float(QMeV),
  'stationary_electron_Klein_Nishina_over_Thomson':float(kn),
  'initial_mixture_stationary_electron_Compton_m2_kg':float(kn*SIGT*(8/15)/MU),
  'qualification':'Single-scattering scale only. Dense-plasma spectral deposition requires a transport calculation; pair production is omitted.'},
 'selected_fiction_and_conservation_requirements':required,'input_and_extraction_checks':checks}
(OUT/'small-pulse-check.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
print(json.dumps(record,indent=2))
