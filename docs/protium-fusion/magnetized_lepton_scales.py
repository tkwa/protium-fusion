"""Map encounter scales; logarithmic intervals are not collision-rate bounds."""
from pathlib import Path
from math import pi,sqrt,log
import argparse,json,gzip
import numpy as np
from electron_wall_orbits import Cylinder,ME,MP,E_CHARGE as E,EPS0
from lepton_collision_screen import local_state
HBAR=1.054571817e-34
ROOT=Path(__file__).parent

def evaluate(root):
 p=next(p for p in json.loads((root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
 c=Cylinder(p['selected_source_profiles']['22049']);eq=c.eq
 prior=json.loads(gzip.decompress((root/'lepton-collision-screen.json.gz').read_bytes()))
 rows=[]
 for old in prior['cases']:
  if old['log_lambda']!=15 or not old['include_relative_flow']:continue
  f=old['reaction_fraction_inside'];r=old['local_state']['radius_m'];local=local_state(c,r)
  x=r/c.scale;psip=float(c.spline(x,1));psipp=float(c.spline(x,2))
  slopes=[(-psip-2*eq['ah']*x)/c.scale,(-eq['Z']*psip-2*eq['ac']*x)/c.scale]
  omegae=E*c.B/ME;energy=old['mean_energy_eV'];vesq=2*energy*E/ME
  rhoe=sqrt(2*vesq/3)/omegae
  positive=local['proton_density_m-3']+eq['Z']*local['catalyst_density_m-3']
  result={'reaction_fraction_inside':f,'radius_m':r,'local_state':{k:v for k,v in local.items() if k!='species'},
   'reference_electron_mean_energy_eV':energy,'electron_omega_rad_s':omegae,
   'electron_rms_perpendicular_gyroradius_m':rhoe,'electron_Landau_level_eV':HBAR*omegae/E,
   'electron_beta_if_ne_equals_positive_charge_density':omegae/sqrt(positive*E*E/(EPS0*ME)),
   'radial_E_gradient_length_m':abs(c.scale*psip/psipp),
   'proton_density_gradient_length_m':1/abs(slopes[0]),'catalyst_density_gradient_length_m':1/abs(slopes[1]),'species':[]}
  for (mass,n,z,temp),name in zip(local['species'],['H','Ti49']):
   mui=ME*mass/(ME+mass);vt=sqrt(2*temp*E/mass);omega=z*E*c.B/mass
   vrel=sqrt(vesq+local['relative_azimuthal_m_s']**2+3*temp*E/mass)
   vcold=sqrt(local['relative_azimuthal_m_s']**2+3*temp*E/mass)
   b90=z*E*E/(4*pi*EPS0*mui*vrel*vrel);bq=HBAR/(mui*vrel)
   bmin=max(b90,bq);bgyro=vrel/omegae;lam=local['ion_static_screening_m']
   spacing=(3/(4*pi*n))**(1/3);wp=sqrt(n*(z*E)**2/(EPS0*mass))
   result['species'].append({'name':name,'mass_kg':mass,'density_m-3':n,'charge':z,
    'ion_omega_rad_s':omega,'ion_individual_plasma_omega_rad_s':wp,'ion_beta':omega/wp,
    'ion_thermal_speed_m_s':vt,'ion_rms_perpendicular_gyroradius_m':vt/omega,
    'ion_gyroradius_over_static_screening':vt/omega/lam,
    'ion_Wigner_Seitz_radius_m':spacing,'ion_coupling_Gamma':(z*E)**2/(4*pi*EPS0*spacing*temp*E),
    'reference_relative_rms_speed_m_s':vrel,'classical_b90_m':b90,'reduced_deBroglie_m':bq,
    'cold_electron_relative_rms_speed_m_s':vcold,
    'ion_turn_during_cold_screening_encounter_rad':omega*lam/vcold,
    'reference_short_cutoff_m':bmin,'electron_one_radian_encounter_scale_m':bgyro,
    'screening_encounter_time_s':lam/vrel,'electron_turn_during_screening_encounter_rad':omegae*lam/vrel,
    'ion_turn_during_screening_encounter_rad':omega*lam/vrel,
    'static_log_interval_total':log(lam/bmin),
    'static_log_interval_below_0p1_electron_radian':max(0.,log(.1*bgyro/bmin)),
    'static_log_interval_above_one_electron_radian':max(0.,log(lam/bgyro))})
  rows.append(result)
 return {'scope':__doc__,'conventions':['Ion screening is the local static susceptibility length already used by the equilibrium; dynamic response may differ.','Relative speed is an RMS diagnostic using the earlier isotropic electron energy, ion thermal variance and transverse mean drift.','Quantum comparison uses reduced de Broglie wavelength hbar/(reduced mass * relative speed); this is a scale convention, not a fitted quantum collision operator.','No scalar logarithmic interval is assigned to a magnetized tensor or claimed as a bound on heating/transport.'],'cases':rows}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);a=ap.parse_args()
 out=evaluate(a.root);(a.outdir/'magnetized-lepton-scales.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 for row in out['cases']:
  print('fraction',row['reaction_fraction_inside'],'r',row['radius_m'],'rho_e',row['electron_rms_perpendicular_gyroradius_m'],'L_nTi',row['catalyst_density_gradient_length_m'],'L_E',row['radial_E_gradient_length_m'])
  for s in row['species']:print(s['name'],{k:v for k,v in s.items() if k in ['ion_beta','ion_coupling_Gamma','ion_gyroradius_over_static_screening','classical_b90_m','reduced_deBroglie_m','electron_one_radian_encounter_scale_m','ion_turn_during_screening_encounter_rad','static_log_interval_total','static_log_interval_below_0p1_electron_radian','static_log_interval_above_one_electron_radian']})
