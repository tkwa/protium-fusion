"""Independent invariants and quadrature checks for the C13 beam screen."""
from pathlib import Path
import json
from math import pi
import numpy as np
from scipy.integrate import quad
from nonneutral_screen import (Q,MC,MP,COULOMB,DATA,nuclear_and_drag,energy_case,
 confinement,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N)

checks={}
# Close the isotope cycle at arbitrary p,d substitution for neutron capture.
reactions=[({'c13':-1,'p':-1,'n13':1,'n':1},1),
           ({'n13':-1,'c13':1},1),({'d':-2,'he4':1},.5)]
for y in [0,.2,.61,1.]:
 seq=reactions+[({'c13':-1,'p':-1,'c12':1,'d':1},y),
    ({'c12':-1,'n':-1,'c13':1},y),({'n':-1,'p':-1,'d':1},1-y)]
 net={}
 for d,n in seq:
  for key,v in d.items():net[key]=net.get(key,0)+v*n
 assert abs(net['p']+2)<1e-14 and abs(net['he4']-.5)<1e-14
 assert max(abs(v) for k,v in net.items() if k not in ['p','he4'])<1e-14
checks['isotope_cycle_closed_for_pd_substitution']=[0,.2,.61,1.]
checks['pd_and_regeneration_Q_closure_MeV']=Q['pd']+Q['c12_capture']-Q['np_capture']
assert abs(checks['pd_and_regeneration_Q_closure_MeV'])<1e-8
assert abs(DD_CHARGED_PER_NEW_N+DD_THERMAL_PER_NEW_N-.5*Q['net_D'])<1e-12

# Numerical impact-parameter integral, independent of the Coulomb-log formula.
E=4.5; ratio=MC/MP; ecm=ratio/(1+ratio)*E; b0=6*COULOMB/(2*ecm)
f=4*ratio/(1+ratio)**2; logL=25.
exact=quad(lambda x:2*pi*f*b0*b0*np.exp(2*x)/(1+np.exp(2*x)),-30,logL,
           epsabs=1e-9,epsrel=1e-11)[0]/100
approx=2*pi*f*b0*b0*logL/100
checks['rutherford_energy_quadrature_relative_error']=abs(exact/approx-1)
assert checks['rutherford_energy_quadrature_relative_error']<1e-9

def theta_lab(x):
 theta=2*np.arctan(np.exp(-x))
 return np.arctan2(ratio*np.sin(theta),1+ratio*np.cos(theta))
ang=quad(lambda x:2*pi*b0*b0*np.exp(2*x)*theta_lab(x)**2,-25,logL,
         epsabs=1e-8,epsrel=1e-10)[0]/100
ang_approx=8*pi*(ratio/(1+ratio))**2*b0*b0*logL/100
checks['angular_variance_small_angle_relative_error']=abs(ang_approx/ang-1)
assert checks['angular_variance_small_angle_relative_error']<.05

# Independent polynomial root check of the conducting-cylinder cold-fluid
# instability threshold. Normalize frequencies by omega_p and speeds by
# omega_p/k_perp so only ratios remain in the root calculation.
a,b=1.,1.2
threshold=(a**(2/3)+b**(2/3))**1.5
roots={}
for factor in [.8,1.2]:
 V=factor*threshold; ap=a*a;ac=b*b
 rr=np.roots([1,-2*V,V*V-ap-ac,2*ac*V,-ac*V*V])
 roots[str(factor)]=float(max(abs(rr.imag)))
assert roots['0.8']>.1 and roots['1.2']<1e-8
checks['cold_pipe_polynomial_max_imaginary_root']=roots

# Recover the monoenergetic result as the Gaussian width tends to zero,
# away from knots in the evaluated pn curve.
zero=nuclear_and_drag(4.501,0);narrow=nuclear_and_drag(4.501,1e-5)
checks['narrow_spectrum_drag_relative_error']=abs(narrow['ion_recoil_drag_MeV_per_pn']/zero['ion_recoil_drag_MeV_per_pn']-1)
assert checks['narrow_spectrum_drag_relative_error']<1e-7
checks['maximum_energy_closure_MeV']=0.
for E in np.linspace(3.4,7.,20):
 for sd in [0,.1]:
  d=nuclear_and_drag(float(E),sd)
  checks['maximum_energy_closure_MeV']=max(checks['maximum_energy_closure_MeV'],abs(d['energy_closure_MeV']))
  # A perfect converter and accelerator recover nuclear Q irrespective of
  # scattering or pd substitution. Side-cost and auxiliary inputs are zero.
  ideal=energy_case(d,1,1,1,0,0)['net_MeV_per_pn_before_unmodeled_plant_loads']
  assert abs(ideal-Q['heat_per_new_neutron'])<1e-8
assert checks['maximum_energy_closure_MeV']<1e-8

# Report rather than suppress the source's independent-production mismatch.
xs=DATA['channels']
def sig(za,e):
 d=xs[f'{za}_0'];return float(np.interp(e,d['energy_MeV'],d['sigma_barn']))
checks['source_independent_production_mismatch_at_5_MeV_barn']={
 'D':sig(1002,5.),'C12':sig(6012,5.),'difference':sig(1002,5.)-sig(6012,5.)}
assert checks['source_independent_production_mismatch_at_5_MeV_barn']['difference']>.01
out=Path(__file__).with_name('nonneutral-validation.json');out.write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
