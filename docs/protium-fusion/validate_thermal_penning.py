"""Independent integration and limiting-case checks, not experimental validation."""
from pathlib import Path
import json
from math import gamma,pi,sqrt
import numpy as np
from scipy.special import gammaincc
from scipy.integrate import quad,simpson
from thermal_neutron_screen import (DATA,ALPHA,thermal_integrals,moments,case,
                                  MASS,MEV,MU_PC,neutron_kinetic,Q)
from penning_equilibrium import equilibrium,at_wall,MRC,KAPPA
ROOT=Path(__file__).parent


def exact_linear(T,za,p):
 d=DATA['channels'][f'{za}_0'];e=ALPHA*np.array(d['energy_MeV']);s=np.array(d['sigma_barn'])
 lo=e[:-1];hi=e[1:];a=np.diff(s)/np.diff(e);b=s[:-1]-a*lo
 def integ(n):return T**(n+1)*gamma(n+1)*(gammaincc(n+1,lo/T)-gammaincc(n+1,hi/T))
 return float(sum(a*integ(p+1)+b*integ(p)))


def main():
 errs=[]
 for T in [.1,.3,.6,1.,3.]:
  for za in [1,1002,5010,6012,6013,7013,7014]:
   for power in [0,1]:
    n=thermal_integrals(T,za,power);a=exact_linear(T,za,power+1)
    if a:errs.append(abs(n/a-1))
 maxerr=max(errs);assert maxerr<2e-8,maxerr
 # Maxwellian normalization in SI and reaction-biased constant-cross-section E.
 T=.6;sigma=1e-28;kT=T*MEV
 norm=quad(lambda e:e*np.exp(-e/T),0,np.inf,epsabs=1e-12)[0]
 sv=sqrt(8*MEV/(pi*MU_PC))*T**-1.5*norm*sigma
 analytic=sigma*sqrt(8*kT/(pi*MU_PC));assert abs(sv/analytic-1)<1e-11
 # A rare-source Monte Carlo with an importance distribution for relative E.
 # Common-T COM motion has independent mean energy 3T/2.
 rng=np.random.default_rng(192301);n=1000000;proposal=1.7
 er=rng.gamma(1.5,proposal,n)
 d=DATA['channels']['7013_0'];sig=np.interp(er/ALPHA,d['energy_MeV'],d['sigma_barn'],left=0,right=0)
 w=sig*np.sqrt(er)*(proposal/T)**1.5*np.exp(-er*(1/T-1/proposal))
 kcom=rng.gamma(1.5,T,n);kn=neutron_kinetic(T,er,Q['pn'])
 mc_mean=float(sum(w*er)/sum(w));mc_pair=float(sum(w*(er+kcom))/sum(w))
 exact=moments(T,7013)['mean_reacting_Ecm_MeV'];assert abs(mc_mean/exact-1)<.003
 assert abs(mc_pair/(exact+1.5*T)-1)<.003
 # Weak charge Gaussian limit independently determined by the external well.
 dilute=equilibrium(1e-10,1e-10,.4,20,nout=3001)
 x=dilute['x'];h_analytic=dilute['h0']*np.exp(-dilute['ah']*x*x)
 c_analytic=dilute['c0']*np.exp(-dilute['ac']*x*x)
 dilute_error=max(np.max(abs(dilute['h']/h_analytic-1)),np.max(abs(dilute['c']/c_analytic-1)))
 assert dilute_error<1e-6
 # Integrate the solved densities independently to verify Gauss's law and
 # compare tighter ODE tolerances before checking field-scaling identities.
 h0=3/MRC/2;c0=3/MRC/12
 eq=equilibrium(h0,c0,.5,40,rtol=1e-9,nout=8001)
 eq2=equilibrium(h0,c0,.5,40,rtol=1e-11,nout=8001)
 x=eq['x'];q_integral=simpson((eq['h']+6*eq['c'])*x,x=x)
 gauss=abs((-x[-1]*eq['dpsi'][-1])/q_integral-1);assert gauss<1e-8
 tight=abs(eq['ihc'][-1]/eq2['ihc'][-1]-1);assert tight<1e-7
 useful=(eq['h']>1e-150)&(eq['c']>1e-150)
 observed=np.log(eq['c'][useful]/c0)-6*np.log(eq['h'][useful]/h0)
 predicted=.5*KAPPA**2*.5**2*(MRC-6)*x[useful]**2
 separation=float(np.max(abs(observed-predicted)));assert separation<1e-10
 a=at_wall(eq,24,B=100);b=at_wall(eq,24,B=200)
 scaling={'radius':abs(b['radius_m']/a['radius_m']*2-1),
          'rate':abs(b['reaction_rate_per_m_s']/a['reaction_rate_per_m_s']/4-1),
          'mass':abs(b['field_mass_per_m_kg']/a['field_mass_per_m_kg']-1),
          'wall_circulation':abs(b['wall_ion_energy_circulated_MeV_per_pn']/a['wall_ion_energy_circulated_MeV_per_pn']*2-1),
          'vacuum_radiation_per_event':abs(b['vacuum_ion_Larmor_MeV_per_pn']/a['vacuum_ion_Larmor_MeV_per_pn']-1)}
 assert max(scaling.values())<1e-12
 closure=max(abs(case(T)['energy_closure_MeV']) for T in [.1,.3,.5,.6,.8,1.,1.5,2.])
 assert closure<1e-8
 out={'all_checks_passed':True,'scope':'Model mathematics and specified energy ledger only; nuclear-channel consistency and reactor engineering remain unvalidated.',
 'max_gauss_legendre_vs_gamma_integral_relative_error':maxerr,
 'constant_sigma_Maxwellian_relative_error':abs(sv/analytic-1),
 'reaction_energy_Monte_Carlo_MeV':mc_mean,'reaction_energy_quadrature_MeV':exact,
 'reaction_pair_energy_Monte_Carlo_MeV':mc_pair,'reaction_pair_energy_analytic_MeV':exact+1.5*T,
 'dilute_Gaussian_max_relative_error':float(dilute_error),'Gauss_law_relative_error':float(gauss),
 'overlap_tolerance_relative_error':float(tight),'charge_mass_separation_identity_error':separation,
 'field_scaling_relative_errors':scaling,'maximum_nuclear_energy_closure_MeV':closure}
 (ROOT/'thermal-penning-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(out,indent=2))

if __name__=='__main__':main()
