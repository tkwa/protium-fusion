"""Low-energy N14(p,gamma) relevance screen; not an evaluated updated rate.

Adds a deliberately specified constant delta S=0.25 keV barn, either only
below the newly measured maximum Ecm, or over all energies as a sensitivity.
The latter is not the Chen R-matrix fit. Resonant/nonresonant REACLIB set
fractions are fit diagnostics, not a decomposition into measured channels.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
from scipy.integrate import quad
from scipy.constants import alpha,atomic_mass,c,e,Avogadro,k

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 data=json.loads(a.manifest.read_text());rate=next(r for r in data['rates'] if r['name']=='N14_p_to_O15_reaclib')
 mu_u=14/15;mu=mu_u*atomic_mass;MeV=1e6*e;S_J_m2=.25*1e3*e*1e-28;Eg=2*mu*c*c*(math.pi*alpha*7)**2/MeV;Emax=.260*14/15
 def delta_rate(T9,upper):
  kT=k*T9*1e9/MeV
  integral,error=quad(lambda E:math.exp(-E/kT-math.sqrt(Eg/E)) if E>0 else 0,0,upper,epsabs=1e-35,epsrel=1e-10,limit=300)
  return Avogadro*math.sqrt(8/(math.pi*mu))*(kT*MeV)**-1.5*S_J_m2*integral*MeV*1e6
 rows=[]
 for t in [.015,.1,.2,.3,.5,.684,1.,1.058,1.147,1.5,2.]:
  tf=np.array([1,1/t,t**(-1/3),t**(1/3),t,t**(5/3),math.log(t)])
  terms=[math.exp(np.dot(s['coefficients_a0_through_a6'],tf)) for s in rate['sets']];total=sum(terms);nr=sum(v for v,s in zip(terms,rate['sets']) if not s['resonant']);kt=k*t*1e9/MeV
  sub=delta_rate(t,Emax);allenergy=delta_rate(t,10.)
  rows.append({'T_GK':t,'saved_total_cm3_mol_s':total,'reaclib_set_values':terms,'nonresonant_labeled_fit_fraction':nr/total,'constant_deltaS_below_Emax_rate':sub,'constant_deltaS_below_Emax_relative_rate':sub/total,'constant_deltaS_all_energy_rate':allenergy,'constant_deltaS_all_energy_relative_rate':allenergy/total,'constant_S_kernel_fraction_below_Emax':sub/allenergy,'Gamow_peak_nonresonant_MeV':(Eg*kt*kt/4)**(1/3)})
 out={'scope':__doc__,'source_manifest_sha256':hashlib.sha256(a.manifest.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reaction':'N14(p,gamma)O15','mu_amu_mass_number_approximation':mu_u,'deltaS_keV_b':.25,'new_measurement_Ep_lab_MeV':[.110,.260],'new_measurement_Ecm_MeV':[.110*14/15,Emax],'rows':rows,'sources':{'published':'https://doi.org/10.1103/q756-hzmt','public_preprint':'https://arxiv.org/html/2410.16086v1'},'version_note':'Published abstract S(0)=1.93 +/-0.10 keV b; public preprint S(0)=1.92 +/-0.08. These are not interchangeable. Full published supplement not accessed; public preprint full text/supplement examined.','interpretation':'No updated hot rate is claimed. The finite-window screen shows why changing S(0) by15% does not justify uniformly scaling GK rates. High-energy R-matrix interference and resonance revisions require the full evaluated cross-section/rate, not this additive diagnostic.'}
 a.out.write_text(json.dumps(out,indent=2)+'\n')
 for r in rows:print(r['T_GK'],r['saved_total_cm3_mol_s'],r['nonresonant_labeled_fit_fraction'],r['constant_deltaS_below_Emax_relative_rate'],r['constant_deltaS_all_energy_relative_rate'])
if __name__=='__main__':main()
