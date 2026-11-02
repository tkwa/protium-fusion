"""Frozen-charge-moment ion-sphere stress added consistently to the saved EOS.
All state inputs/outputs are cgs, matching FastEOS. No production file is changed.
F_C = u_C = -C rho^(1/3) Ye^(4/3), independent of temperature.
The charge-weighted <Z^(2/3)> is frozen from the initial fuel mixture.
This over-attractive strong-coupling scale is a sensitivity, not a precision EOS.
"""
from pathlib import Path
import sys
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from fast_eos import FastEOS,NA
E2_CGS=2.307077552e-19 # e^2/(4 pi eps0), erg cm

class NonidealEOS(FastEOS):
 def __init__(self,charge_moment,strength=1.,**kwargs):
  super().__init__(**kwargs);self.charge_moment=float(charge_moment);self.strength=float(strength)
  self.coefficient=.9*E2_CGS*NA*(4*np.pi*NA/3)**(1/3)*self.charge_moment*self.strength
  self.calls=0;self.minimum_pressure=float('inf');self.minimum_cs2=float('inf');self.minimum_cv=float('inf');self.minimum_energy=float('inf')
 def correction(self,rho,ye):
  u=-self.coefficient*np.asarray(rho)**(1/3)*np.asarray(ye)**(4/3)
  p=np.asarray(rho)*u/3
  return u,p
 def state(self,rho,T,ye,yi):
  s=super().state(rho,T,ye,yi);u,p=self.correction(rho,ye)
  s.e=s.e+u;s.p=s.p+p;s.de_drho=s.de_drho+u/(3*rho);s.dp_drho=s.dp_drho+4*p/(3*rho)
  # p_T, u_T and c_v are unchanged because F_C has no T dependence.
  cs2=s.dp_drho+s.dp_dT*(s.p/rho**2-s.de_drho)/s.c_v
  s.gamma1=rho*cs2/s.p
  self.calls+=1;self.minimum_pressure=min(self.minimum_pressure,float(np.min(s.p)));self.minimum_cs2=min(self.minimum_cs2,float(np.min(cs2)));self.minimum_cv=min(self.minimum_cv,float(np.min(s.c_v)));self.minimum_energy=min(self.minimum_energy,float(np.min(s.e)))
  if np.any(s.p<=0) or np.any(cs2<=0) or np.any(s.c_v<=0):
   raise ValueError('Coulomb stress left positive-pressure/compressibility/heat-capacity domain')
  return s
 def diagnostics(self):
  return {'state_calls':self.calls,'minimum_pressure_dyn_cm2':self.minimum_pressure,'minimum_adiabatic_sound_speed_squared_cm2_s2':self.minimum_cs2,'minimum_cv_erg_g_K':self.minimum_cv,'minimum_specific_energy_erg_g':self.minimum_energy,'charge_weighted_Z_two_thirds':self.charge_moment,'strength':self.strength,'coefficient_cgs':self.coefficient}

def validate():
 eos=NonidealEOS(3.);rho=np.geomspace(.1,1e5,31);ye=np.linspace(.45,.7,31);h=1e-5
 u,p=eos.correction(rho,ye)
 ur=(eos.correction(rho*(1+h),ye)[0]-eos.correction(rho*(1-h),ye)[0])/(2*h*rho)
 pr=(eos.correction(rho*(1+h),ye)[1]-eos.correction(rho*(1-h),ye)[1])/(2*h*rho)
 uy=(eos.correction(rho,ye*(1+h))[0]-eos.correction(rho,ye*(1-h))[0])/(2*h*ye)
 checks={'density_energy_derivative_relative':float(np.max(abs(ur/(u/(3*rho))-1))),'density_pressure_derivative_relative':float(np.max(abs(pr/(4*p/(3*rho))-1))),'composition_derivative_relative':float(np.max(abs(uy/((4/3)*u/ye)-1))),'thermodynamic_identity_relative':float(np.max(abs((p-rho*rho*u/(3*rho))/p))),'archived_composition_identity_relative':float(np.max(abs((u+rho*u/(3*rho))/ye/((4/3)*u/ye)-1)))}
 assert max(checks.values())<1e-8,checks
 # Compare SI ion-sphere energy from independently stated constants.
 rhosi=rho*1000;mu=1.66053906660e-27;ae=(3/(4*np.pi*rhosi*ye/mu))**(1/3);usi=-.9*2.307077552e-28/ae*ye*eos.charge_moment/mu
 checks['independent_SI_energy_relative']=float(np.max(abs(u*1e-4/usi-1)));assert checks['independent_SI_energy_relative']<1e-9
 return {'passed':True,'checks':checks,'scope':'Finite differences of added free-energy term, exact Maxwell/composition identities, and independent SI normalization. Does not certify the source electron spline EOS or accuracy of strong-coupling approximation.'}
