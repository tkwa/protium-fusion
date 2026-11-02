#!/usr/bin/env python3
"""Rigid-electron-background ion-sphere diagnostic on saved initial profiles.
U_C=-0.9 sum(n_i Z_i^(5/3)) e²/(4pi eps0 a_e), P_C=U_C/3.
This is a strong-coupling scale, not a complete nonideal EOS or uncertainty bound.
"""
import argparse,hashlib,json,math,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path.insert(0,str(ROOT))
from fast_eos import FastEOS
MU=1.66053906660e-27;KE2=2.307077552e-28;KB=1.380649e-23
EOS=FastEOS()
def case(rho,T,mass,K,A,Z,x):
 ye=(Z+x)/(A+x);yi=(1+x)/(A+x);z53=(Z**(5/3)+x)/(A+x);ne=rho*ye/MU;ae=(3/(4*math.pi*ne))**(1/3)
 uc=-.9*KE2/ae*z53/MU;pc=rho*uc/3
 eos=EOS.state(rho/1000,T,np.full(len(T),ye),np.full(len(T),yi));P=eos.p*.1;U=float(mass@(eos.e*1e-4));Ec=float(mass@uc)
 ideal=rho*KB*T*(ye+yi)/MU;pi=rho*KB*T*yi/MU;gamma=Z**(5/3)*KE2/ae/(KB*T)
 return {'H_per_catalyst':x,'initial_internal_J_ideal_ion_Fermi_electron_LTE_radiation':U,'unchanged_initial_kinetic_J':K,'prepared_J':U+K,'ion_sphere_correction_J':Ec,'ion_sphere_per_total_prepared':Ec/(U+K),'ion_sphere_per_internal':Ec/U,'mass_weighted_abs_Pc_over_EOS_pressure':float(mass@(abs(pc)/P)/sum(mass)),'max_abs_Pc_over_EOS_pressure':float(max(abs(pc)/P)),'volume_weighted_absPc_over_P_ratio':float(sum(abs(pc)*mass/rho)/sum(P*mass/rho)),'mass_weighted_abs_Pc_over_nondegenerate_thermal_pressure':float(mass@(abs(pc)/ideal)/sum(mass)),'mass_weighted_abs_Pc_over_ion_thermal_pressure':float(mass@(abs(pc)/pi)/sum(mass)),'mass_weighted_catalyst_Gamma':float(mass@gamma/sum(mass)),'catalyst_Gamma_min_max':[float(min(gamma)),float(max(gamma))],'scope':'Same initial rho,T,radial velocity and totalmass. Newmix internal/preparedenergy recalculated; source pressure is EOS total material plusLTE radiation, including degeneracy. Pressure correction is rigid-background ion-sphere density derivative, not complete electron-ion exchange/correlation.'}
def main():
 paths={'A':('/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-A-compression_cap',14,7,1),'B1':('/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-B1-compression_cap',15,7,1),'C13':('/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-C13-compression_cap',13,6,1),'B2':('/var/tmp/tkwa-chatgpt/protium-b2-completion/B2-64-compression-cap',12,6,2)}
 rows={};checks=[]
 for stage,(stem,A,Z,x0) in paths.items():
  d=json.load(open(stem+'.json'));f=np.load(stem+'.npz');nn=len(f['species']);nz=int(f['zones']);nv=int(f['variables_per_zone']);z=f['state'][:,0].reshape(nz,nv);m=f['cell_mass_kg'];r=np.r_[0.,z[:,nn+1]*d['R0_m']];rho=m/(4*math.pi/3*np.diff(r**3));T=z[:,nn]*d['T0_K'];v=z[:,nn+2]*float(f['velocity_unit_cm_s'])*.01;nodem=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K=float(.5*nodem@(v*v))
  cases=[case(rho,T,m,K,A,Z,x) for x in sorted(set([x0,2,4]))];original=next(c for c in cases if c['H_per_catalyst']==x0)
  checks.extend([abs(K/d['initial_kinetic_energy_J']-1)<1e-12,abs(original['initial_internal_J_ideal_ion_Fermi_electron_LTE_radiation']/d['initial_internal_energy_J']-1)<1e-10])
  for c in cases:c['ion_sphere_per_reference_fixed_prepared']=c['ion_sphere_correction_J']/(original['initial_internal_J_ideal_ion_Fermi_electron_LTE_radiation']+K)
  rows[stage]={'initial_ratio':x0,'cases':cases,'input_sha256':{s:hashlib.sha256(Path(stem+s).read_bytes()).hexdigest() for s in ['.json','.npz']}}
 return {'scope':__doc__,'sources':{'ion_liquid_EOS':'https://doi.org/10.1103/PhysRevE.62.8554','strong_coupling_context':'https://academic.oup.com/mnras/article/490/4/5839/5610243'},'rows':rows,'checks':{'count':len(checks),'passed':sum(checks),'all_passed':all(checks)},'source_EOS_sha256':hashlib.sha256((ROOT/'electron-eos.npz').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':['The -0.9Gamma ion-sphere term overstates attraction versus the full liquidOCP at moderatecoupling; it is a scale not rigorous upperbound.','No quantumion, electronscreening, exchange, boundstates or composition-dependent ionization calculation.','Local pressurechanges can matter even when volumeintegrated correction is small; no yieldresponse inferred.']}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default=str(Path(__file__).with_name('coulomb-initial.json')));a=p.parse_args();o=main();Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
