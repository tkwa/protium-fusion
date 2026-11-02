"""Bounded photonuclear importance screen, not coupled shower/neutron transport."""
from pathlib import Path
import argparse,json,math,hashlib
NA=6.02214076e23;A_G_MOL=12.011;C13=.0107
NIST=[(5.,.02708),(6.,.02469),(8.,.02154),(20.,.01575)]
def mu(E):
 for (a,x),(b,y) in zip(NIST,NIST[1:]):
  if a<=E<=b:return math.exp(math.log(x)+(math.log(E/a)/math.log(b/a))*math.log(y/x))
 raise ValueError(E)
def first(E,sigma_mb,fraction):
 em_b=mu(E)*A_G_MOL/NA/1e-24;nuc_b=sigma_mb*.001*fraction
 return {'E_MeV':E,'cross_section_mb_per_target_isotope':sigma_mb,'target_number_fraction':fraction,'EM_mu_cm2_g':mu(E),'EM_cross_section_barn_per_carbon':em_b,'nuclear_branch_first_interaction':nuc_b/(em_b+nuc_b)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--prompt_fraction',type=float,default=.01406);a=ap.parse_args();f=a.prompt_fraction
 measured=[dict(first(6.4,.0941,C13),source='Edge1960, effective mixed gamma energy; uncertainty0.0100mb'),dict(first(6.61,.32,C13),source='Green/Donahue1964, weighted mixed gamma energy; uncertainty0.04mb')]
 low=[first(5.657,s,C13) for s in [.0941,.32,1.,10.]];high=[first(20.,s,1.) for s in [5.,10.,20.]]
 for row in high:row['prompt_fraction_times_first_branch']=f*row['nuclear_branch_first_interaction']
 d={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'natural_C13_number_fraction_assumed':C13,'nist_table_E_MeV_mu_cm2_g':NIST,'measured_nearby_examples':measured,'unmeasured_5p657_sensitivity':low,'unmeasured_20MeV_total_nuclear_sensitivity':high,'prompt_fraction_of_H_plus_U_scenario':f,'full_prompt_loss_thermal_reserve':f,'capture_if_all_other_energy_retained':1-f,'capture_if_rest_has_99percent_capture':.99*(1-f),'qualification':'First interaction in a homogeneous effectively thick slab only. Values5/10/20mb and assigning nearby6.4/6.61MeV cross sections to5.657MeV are sensitivities, not evaluated cross sections or confidence limits. Repeated Compton scattering, bremsstrahlung, secondary photons, nuclear branching and neutron capture/escape require coupled transport. The prompt fraction is a caller-supplied scenario, not recomputed for each hydro family. C13(gamma,n)C12 remains usable carbon; terminal isotope conversion is a separate ledger.','sources':{'Edge1960':'https://doi.org/10.1103/PhysRev.119.1643','GreenDonahue1964':'https://doi.org/10.1103/PhysRev.135.B701','Fultz1966':'https://doi.org/10.1103/PhysRev.143.790','Carchon1976':'https://doi.org/10.1103/PhysRevC.14.456','NIST':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html','original_figures_facsimile':'https://www.govinfo.gov/content/pkg/GOVPUB-C13-518f4d604594649e67f826fc87f6550e/pdf/GOVPUB-C13-518f4d604594649e67f826fc87f6550e.pdf'}}
 assert all(0<r['nuclear_branch_first_interaction']<1 for r in measured+low+high)
 a.out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'nearby':measured,'5p657':low,'20MeV':high,'reserve':f},indent=2))
if __name__=='__main__':main()
