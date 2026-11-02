"""Compare interpolation of the cold surrogate with its source free-electron EOS.
This validates interpolation only, not a cold material or its ionization state.
"""
from pathlib import Path
import json
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1];sys.path.insert(0,str(ROOT))
from fast_eos import FastEOS,KB_CGS,NA,ARAD_CGS
import pynucastro as p
from pynucastro.eos.electron_eos import ElectronEOS

def main():
    original=np.load(ROOT/'electron-eos.npz');new=np.load(HERE/'electron-eos-cold-surrogate.npz')
    n_old=len(original['log_temperature']);n_new=len(new['log_temperature']);offset=n_new-n_old
    checks=[]
    def check(name,value,tol):checks.append({'name':name,'error':float(value),'maximum':tol,'passed':bool(value<=tol)})
    for key in ('log_pressure','log_energy_density'):
        check('hot nodes unchanged '+key,float(np.max(abs(new[key][:,offset:]-original[key]))),0.)
    model=FastEOS(HERE/'electron-eos-cold-surrogate.npz')
    source=ElectronEOS();composition=p.Composition(['p']);composition.set_all(1.)
    ye,yi=8/15,2/15
    def exact(rho,T):
        ee,pp=source.pe_state(rho*ye,T,composition,compute_derivs=False,eta_guess_max=1e10)
        pressure=ee.p+pp.p+rho*KB_CGS*NA*T*yi+ARAD_CGS*T**4/3
        energy=(ee.e+pp.e)*ye+1.5*KB_CGS*NA*T*yi+ARAD_CGS*T**4/rho
        return pressure,energy
    points=[]
    rng=np.random.default_rng(82419)
    for logrho,logT in zip(rng.uniform(-2,6,21),rng.uniform(2.5,5.2,21)):
        rho,T=10**logrho,10**logT
        st=model.state(rho,T,ye,yi);pressure,energy=exact(rho,T)
        h=.005;pt,et=exact(rho,T*(1+h));pm,em=exact(rho,T*(1-h))
        cv=(et-em)/(2*h*T)
        pr,er=exact(rho*(1+h),T);pl,el=exact(rho*(1-h),T)
        dPdr=(pr-pl)/(2*h*rho)
        errors={'pressure':abs(float(st.p)/pressure-1),'energy':abs(float(st.e)/energy-1),
                'heat_capacity':abs(float(st.c_v)/cv-1),'density_derivative':abs(float(st.dp_drho)/dPdr-1)}
        points.append({'rho_g_cm3':rho,'temperature_K':T,'errors':errors})
    for key,tol in [('pressure',1e-4),('energy',1e-4),('heat_capacity',.005),('density_derivative',.001)]:
        check(key+' maximum of 21 independent source comparisons',max(q['errors'][key] for q in points),tol)
    out={'checks':checks,'points':points,'passed':sum(c['passed'] for c in checks),'total':len(checks),
         'qualification':'Only the interpolation of a fully ionized free-electron/ideal-ion/radiation surrogate is checked; cold-material validity is absent.'}
    (HERE/'cold-eos-validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':out['passed'],'total':out['total'],'checks':checks},indent=2))
    if out['passed']!=out['total']:raise SystemExit(1)

if __name__=='__main__':main()
