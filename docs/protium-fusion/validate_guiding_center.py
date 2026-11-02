"""Independent weak-encounter, Maxwellian, derivative and balance checks."""
from pathlib import Path
import argparse,json,time
from math import pi,sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from guiding_center_collisions import force_tensor_derivatives,species_table,coefficients,state_moments
from guiding_center_encounters import encounter_test
from lepton_collision_screen import ion_coefficients
from cold_magnetized_response import speed_density,conditional_mu2
from electron_wall_orbits import ME,MP,E_CHARGE as EV,EPS0

ROOT=Path(__file__).parent


def independent_maxwell_tensor(w,vth=1.,n=128):
    # Direct integration of (I-u_hat u_hat)/|u| over shifted Gaussian ions.
    # Polar axis is w; relative speed coordinates remove the u=0 singularity.
    x,wx=leggauss(n);mu,wm=leggauss(n)
    high=w/vth+9.;r=(x+1)*high/2;wr=wx*high/2
    weight=np.exp(-r[:,None]**2-(w/vth)**2+2*r[:,None]*w/vth*mu[None,:])
    common=2/(sqrt(pi)*vth)*(wr*r)
    dp=common@(weight@((1-mu*mu)*wm))
    dt=common@(weight@((1+mu*mu)*wm/2))
    return np.array([dp,dt])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    start=time.monotonic();checks=[]
    def check(name,ok,metric=None):checks.append({'name':name,'passed':bool(ok),'metric':metric})
    for w in [0.,1e-4,.01,.3,1.,3.,8.]:
        p,t,pp,tp=force_tensor_derivatives(w,1.,1.)
        direct=independent_maxwell_tensor(w)
        err=float(np.linalg.norm(direct-[p,t])/np.linalg.norm([p,t]))
        check(f'Direct Maxwell impulse covariance w/vth={w}',err<2e-10,err)
        if w>0:
            h=1e-4*max(w,.01)
            vals=[np.array(force_tensor_derivatives(w+i*h,1.,1.)[:2]) for i in [-2,-1,1,2]]
            fd=(vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*h)
            err=float(np.linalg.norm(fd-[pp,tp])/max(np.linalg.norm([pp,tp]),1e-12))
            check(f'Analytic covariance derivative w/vth={w}',err<2e-6,err)
    for mass,z in [(MP,1),(49*MP,22)]:
        species=[(mass,1e19,z,1.6e6)];table=species_table(species,4.)
        G,ratio,vth2=table[0]
        for a in [.01,.3,1.,3.,10.]:
            w=a*sqrt(vth2);p,t,pp,tp=force_tensor_derivatives(w,G,vth2)
            nu,*_=ion_coefficients(w,mass,1e19,z,1.6e6,4.)
            A=.5*(pp+2*(p-t)/w)-ratio/vth2*p*w
            error=float(abs(A/(-nu*w)-1))
            check(f'Full-space divergence plus polarization gives NRL drag Z={z} w/vth={a}',error<2e-8,error)
        for U in [0.,-6e6]:
            for vz in [0.,1e6,1e7,3e7]:
                omega=-2.4e13;c=coefficients(vz,U,table)
                A,D,div=state_moments(vz,U,table,omega)
                # Units are rescaled before evaluating positive semidefiniteness.
                S=np.diag([omega,omega,1.]);scaled=S@D@S
                ev=np.linalg.eigvalsh(scaled)
                check(f'Positive diffusion Z={z}, U={U}, vz={vz}',ev.min()>-1e-12*ev.max(),float(ev.min()/ev.max()))
                # Local formal Gibbs gradient; this is not a global source profile.
                grad=np.array([ME*omega*U/(1.6e6*EV),0.,-ME*vz/(1.6e6*EV)])
                expected=.5*D@grad;actual=A-.5*div
                scale=max(np.linalg.norm(S@expected),1.)
                err=float(np.linalg.norm(S@(actual-expected))/scale)
                check(f'Local fluctuation/polarization balance Z={z}, U={U}, vz={vz}',err<2e-9,err)
            # Independent cold force-pulse quadrature at the same constant log.
            vth=sqrt(2*1.6e6*EV/mass);g=4*pi*1e19*(z*EV**2/(4*pi*EPS0))**2/ME/EV*4
            power=g*quad(lambda u:speed_density(u,vth,abs(U))/u*.5*(1-float(conditional_mu2(np.array([u]),vth,abs(U))[0])),0,abs(U)+10*vth,epsabs=1e-22)[0]
            theory=.5*ME/EV*coefficients(0.,U,table)[4]
            err=float(abs(power/theory-1))
            check(f'Cold axial heating from direct impulse flux Z={z}, U={U}',err<2e-9,err)
    direct=[]
    for eta,omega,nb in [(0.,-1000.,48),(.1,-1000.,48),(.1,1000.,64)]:
        result=encounter_test(eta,omega,nb);direct.append(result)
        for key,tolerance in [('diffusion_relative_error',1e-7),('drift_relative_error',1e-7),('second_order_energy_error',2e-8),('first_order_energy_error',2e-8)]:
            check(f'Two-body variational {key}: mass ratio={eta}, omega={omega}, nodes={nb}',result['solver_success'] and result[key]<tolerance,result[key])
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'independent_encounters':direct,'seconds':time.monotonic()-start}
    (args.outdir/'guiding-center-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass in',out['seconds'],'s')
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
