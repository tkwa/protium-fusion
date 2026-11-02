"""Analytic loss screens for the conditional heavy-ion nonneutral reactor.

Repulsive Coulomb dipole bremsstrahlung follows a classical two-body orbit;
its thermal average agrees with Hardy & Lasenby (2017), appendix B, eq.B.6,
with the emitted vector set to the ordinary photon. Nuclear-size, quantum,
quadrupole and collective corrections remain outside this estimate.
"""
from math import pi,sqrt,exp
import json
import numpy as np
from scipy.integrate import quad
from catalyst_thermal import nucmass,case,MEV,C
from nonneutral_screen import EPS0,E_CHARGE
from pathlib import Path
ROOT=Path(__file__).parent
HBAR=1.054571817e-34
ME=9.1093837139e-31


def ion_brems(Z=22,A=49,T=1.7):
 m1=nucmass(1,1)*MEV/C**2;m2=nucmass(Z,A)*MEV/C**2
 mu=m1*m2/(m1+m2)
 qeff=E_CHARGE*(m2-Z*m1)/(m1+m2)
 K=Z*E_CHARGE**2/(4*pi*EPS0)
 # <v integral(E_gamma d sigma)> in J m^3/s at a common ion temperature.
 coefficient=2/3*qeff*qeff*K*T*MEV/(EPS0*mu*mu*C**3)
 thermal=case(f'{A}'+{6:'C',22:'Ti',23:'V',24:'Cr',25:'Mn'}[Z],T)
 sv=thermal['pn_sigma_v_m3_s']
 vthermal=sqrt(2*T*MEV/mu)
 return {'Z':Z,'A':A,'T_MeV':T,'dipole_effective_charge_e':qeff/E_CHARGE,
 'thermal_power_coefficient_J_m3_s':coefficient,'radiated_MeV_per_pn':coefficient/sv/MEV,
 'Coulomb_Sommerfeld_parameter_at_v_sqrt_2T_mu':K/(HBAR*vthermal),
 'zero_impact_closest_approach_at_E_T_m':K/(T*MEV),
 'relative_deBroglie_reduced_wavelength_m':HBAR/(mu*vthermal)}


def orbit_integral(u):
 if u<.03:
  # ((u^2+3)atan(u)-3u)/u^4, series avoids cancellation.
  return sum(((-1.)**k)*(3/(2*k+1)-1/(2*k-1))*u**(2*k-3) for k in range(2,10))
 return ((u*u+3)*np.arctan(u)-3*u)/u**4


def ends(src,length=20000,barrier_in_kT=25.,electrons_per_pn=1.,electron_clear_s=1e-4):
 T=src['T_MeV'];tj=T*MEV;m1=nucmass(1,1)*MEV/C**2
 nh=src['proton_inventory_per_m'];nt=src['catalyst_inventory_per_m'];rp=src['reaction_rate_per_m_s']
 Hperpn=2*nh*sqrt(tj/(2*pi*m1))*exp(-barrier_in_kT)/(rp*length)
 # This only estimates the thermal part and axial potential of escaping H.
 # Its rotation term and actual end transport require a finite reactor model.
 Ecir=Hperpn*(barrier_in_kT+2)*T
 ef=src['wall_E_V_m'];B=src['B_T'];beta=ef/(C*B)
 xmax=2*ME*ef/(E_CHARGE*(B*B-ef*ef/C**2)) if beta<1 else None
 return {'active_length_m':length,'barrier_kT':barrier_in_kT,'barrier_V':barrier_in_kT*T*1e6,
 'axis_plus_barrier_V':src['axis_to_wall_V']+barrier_in_kT*T*1e6,
 'axial_H_end_losses_per_pn':Hperpn,'end_H_energy_circulation_excluding_rotation_MeV_per_pn':Ecir,
 'electron_birth_per_pn_assumed':electrons_per_pn,'electron_clear_time_s_assumed':electron_clear_s,
 'electron_to_positive_charge_line_density':electrons_per_pn*rp*electron_clear_s/(nh+22*nt),
 'wall_E_over_cB':beta,'planar_cold_secondary_max_penetration_m':xmax,
 'planar_cold_secondary_max_energy_eV':ef*xmax if xmax is not None else None,
 'scope':'Ideal Maxwellian axial flux, prescribed electron residence, and locally planar crossed-field orbit; no finite electrode or kinetic clearing solution.'}

if __name__=='__main__':
 from bundle_solenoid import source
 numerical,error=quad(orbit_integral,0,np.inf,epsabs=1e-11,epsrel=1e-11)
 cases=[ion_brems(T=t) for t in [1.,1.3,1.6,1.7,1.9,2.2]]
 endcases=[ends(source(T=1.6,B=b,heat=.6),length=L) for b,L in [(100,2000),(100,20000),(180,30000)]]
 out={'classical_orbit_integral_numeric':numerical,'classical_orbit_integral_exact':2/3,
 'classical_orbit_integral_relative_error':numerical/(2/3)-1,'cases':cases,'end_and_electron_screens':endcases}
 (ROOT/'ion-radiation-and-ends.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(out,indent=2))
