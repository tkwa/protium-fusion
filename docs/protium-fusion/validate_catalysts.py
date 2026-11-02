"""Independent quadrature/Poisson checks and conserved main-cycle identities."""
from pathlib import Path
from math import gamma,pi,sqrt
import json,gzip
import numpy as np
from scipy.special import gammaincc
from scipy.integrate import simpson
from catalyst_catalog import thermal,nucmass
from catalyst_thermal import CAT,CURVES,case,MEV,C
from catalyst_penning import equilibrium,at_wall
from penning_equilibrium import equilibrium as carbon_equilibrium
ROOT=Path(__file__).parent


def exact_rate(tab,T,Z,A,power=0):
 alpha=nucmass(Z,A)/(nucmass(Z,A)+nucmass(1,1));mu=nucmass(Z,A)*nucmass(1,1)/(nucmass(Z,A)+nucmass(1,1))*MEV/C**2
 e=alpha*np.array(tab['energy_MeV']);s=np.array(tab['sigma_barn'])
 lo=e[:-1];hi=e[1:];slope=np.diff(s)/np.diff(e);offset=s[:-1]-slope*lo
 def integ(p):return T**(p+1)*gamma(p+1)*(gammaincc(p+1,lo/T)-gammaincc(p+1,hi/T))
 val=float(np.sum(slope*integ(2+power)+offset*integ(1+power)))
 return sqrt(8*MEV/(pi*mu))*T**-1.5*val*1e-28


def main():
 errors=[];count=0
 for r in CAT.values():
  z,a=r['Z'],r['A'];tab=CURVES[r['source_file']]['channels'][f'{(z+1)*1000+a}_0']
  for T in [.2,.6,1.,2.5]:
   for p in [0,1]:
    val=thermal(tab,T,z,a,p);ind=exact_rate(tab,T,z,a,p)
    if ind:errors.append(abs(val/ind-1));count+=1
 assert max(errors)<1e-7,max(errors)
 # Reproduce the separately implemented C13 equilibrium at identical inputs.
 old=carbon_equilibrium(.16,.004,.5,40,rtol=1e-10,nout=3001)
 new=equilibrium(6,13,.16,.004,.5,40,rtol=1e-10,nout=3001)
 carbon_diff=max(float(np.max(abs(old[k]-new[k]))) for k in ['psi','dpsi','ih','ic','ihc','ihcx2'])
 assert carbon_diff<1e-11,carbon_diff
 checks=[]
 for name in ['49Ti','51V','53Cr','55Mn']:
  b=json.loads((ROOT/f'catalyst-penning-{name}.json').read_text())['best'][0]
  z,a=CAT[name]['Z'],CAT[name]['A']
  eq=equilibrium(z,a,b['h0'],b['c0'],b['xi'],40,rtol=1e-9,nout=8001)
  tight=equilibrium(z,a,b['h0'],b['c0'],b['xi'],40,rtol=1e-11,nout=8001)
  x=eq['x'];q=simpson((eq['h']+z*eq['c'])*x,x=x)
  gauss=abs((-x[-1]*eq['dpsi'][-1])/q-1)
  tol=abs(eq['ihc'][-1]/tight['ihc'][-1]-1)
  sel=(eq['h']>1e-150)&(eq['c']>1e-150)
  observed=np.log(eq['c'][sel]/eq['c0'])-z*np.log(eq['h'][sel]/eq['h0'])
  predicted=.5*eq['kappa']**2*eq['xi']**2*(eq['M']-z)*x[sel]**2
  separation=float(np.max(abs(observed-predicted)))
  d=at_wall(eq,22,name);d2=at_wall(eq,22,name,B=200);ds=at_wall(eq,22,name,strength=5e6)
  scalings=[abs(d2['reaction_rate_per_m_s']/d['reaction_rate_per_m_s']/4-1),
    abs(d2['radius_m']/d['radius_m']*2-1),abs(d2['field_mass_per_m_kg']/d['field_mass_per_m_kg']-1),
    abs(d2['vacuum_ion_Larmor_MeV_per_pn']/d['vacuum_ion_Larmor_MeV_per_pn']-1),
    abs(ds['field_mass_only_net_W_kg']/d['field_mass_only_net_W_kg']/5-1)]
  closures=[abs(case(name,T,f,g)['energy_closure_MeV']) for T in [.4,1.,2.] for f in [0.,.5,1.] for g in [0.,.5,1.]]
  assert gauss<1e-8 and tol<1e-7 and separation<1e-10 and max(scalings)<1e-12 and max(closures)<1e-8
  checks.append({'target':name,'Gauss_law_relative_error':gauss,'overlap_tolerance_relative_error':tol,
    'charge_mass_identity_error':separation,'max_scaling_error':max(scalings),'maximum_cycle_closure_MeV':max(closures)})
 out={'all_checks_passed':True,'scope':'Mathematics and conditional main-cycle ledger; no experimental trap validation or exclusive side-cycle closure.',
 'independent_integrals_checked':count,'maximum_gauss_vs_incomplete_gamma_error':max(errors),
 'generalized_vs_original_C13_equilibrium_absolute_error':carbon_diff,'equilibrium_checks':checks}
 (ROOT/'catalyst-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
