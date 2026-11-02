"""Independent second-order, softened screened two-body GC encounter test.

Electron m=1; ion mass=1/eta. The ion is unmagnetized on encounter times.
Force is the gradient of -epsilon*exp(-sqrt(r*r+a*a))/sqrt(r*r+a*a).
The O(epsilon) relative displacement drives O(epsilon**2) mean impulse.
This avoids subtracting two nearly equal finite-strength scattering outcomes.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import k1
from numpy.polynomial.legendre import leggauss


def encounter_test(eta=0.,omega=-1000.,nb=48,nphi=8,rtol=3e-9):
    w=np.array([0.,.8,.6]);u=np.linalg.norm(w);a=.1
    nodes,weights=leggauss(nb)
    lo,hi=np.log(1e-4*a),np.log(20.)
    b=np.exp(lo+(nodes+1)*(hi-lo)/2)
    wb=weights*(hi-lo)/2*b*b
    phi=(np.arange(nphi)+.5)*2*np.pi/nphi
    e1=np.array([1.,0.,0.]);e2=np.cross(w/u,e1)
    bhat=e1[:,None]*np.cos(phi)+e2[:,None]*np.sin(phi)
    bvec=(b[None,:,None]*bhat[:,None,:]).reshape(3,-1)
    weights=np.repeat(wb,nphi)*2*np.pi/nphi*u
    n=len(weights)
    def rhs(t,y):
        y=y.reshape(9,n);dr=y[:3];impulse=y[3:6]
        r=bvec+w[:,None]*t;R=np.sqrt(np.sum(r*r,axis=0)+a*a)
        h=np.exp(-R)*(R**-3+R**-2)
        hprime=-np.exp(-R)*(3*R**-5+3*R**-4+R**-3)
        force=-h*r
        second=-h*dr-hprime*r*np.sum(r*dr,axis=0)
        relative=eta*impulse
        relative=relative.copy()
        relative[0]+=force[1]/omega
        relative[1]-=force[0]/omega
        relative[2]+=impulse[2]
        return np.vstack([relative,force,second]).ravel()
    start=time.monotonic()
    sol=solve_ivp(rhs,[-24/u,24/u],np.zeros(9*n),method='DOP853',rtol=rtol,atol=rtol*.01)
    final=sol.y[:,-1].reshape(9,n);I=final[3:6];J=final[6:9]
    D=(I*weights)@I.T;A=J@weights
    Leff=quad(lambda b:b**3/(b*b+a*a)*k1(np.sqrt(b*b+a*a))**2,0,np.inf,epsabs=1e-12)[0]
    G=4*np.pi*Leff
    expectedD=G/u*(np.eye(3)-np.outer(w,w)/u**2)
    derivative=np.array([-G*w[0]*(u*u-3*w[2]**2)/u**5,
                         -G*w[1]*(u*u-3*w[2]**2)/u**5,
                         -3*G*(w[0]**2+w[1]**2)*w[2]/u**5])
    expectedA=.5*derivative-eta*G*w/u**3
    # Energy through second order, including ion recoil, for each impact.
    closure=w@J+.5*(I[2]**2+eta*np.sum(I*I,axis=0))
    scale=np.max(np.abs(w@J))+.5*np.max(I[2]**2+eta*np.sum(I*I,axis=0))
    return {'eta':eta,'omega':omega,'nb':nb,'nphi':nphi,'rtol':rtol,
        'diffusion':D.tolist(),'expected_diffusion':expectedD.tolist(),
        'drift':A.tolist(),'expected_drift':expectedA.tolist(),
        'diffusion_relative_error':float(np.linalg.norm(D-expectedD)/np.linalg.norm(expectedD)),
        'drift_relative_error':float(np.linalg.norm(A-expectedA)/np.linalg.norm(expectedA)),
        'second_order_energy_error':float(np.max(np.abs(closure))/scale),
        'first_order_energy_error':float(np.max(np.abs(w@I))/np.max(np.abs(I))),
        'solver_success':bool(sol.success),'nfev':sol.nfev,'seconds':time.monotonic()-start}


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent)
    args=ap.parse_args()
    results=[]
    for eta,omega,nb in [(0.,-1000.,48),(.1,-1000.,48),(.1,1000.,64)]:
        row=encounter_test(eta,omega,nb)
        results.append(row);print(json.dumps(row),flush=True)
    (args.outdir/'gc-encounter-tests.json').write_text(json.dumps(results,indent=2)+'\n')
