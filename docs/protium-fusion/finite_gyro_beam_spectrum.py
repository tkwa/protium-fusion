"""Spectral operator for monoenergetic straight background ions.

Integrating the resonance delta over polar wavevector angle permits direct
comparison with a binary-encounter ensemble without Maxwellian sampling noise.
"""
from math import pi
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv


def diffusion(vperp,w,Omega,kappa=1.,kcut=4.,G0=4*pi,nk=112,nphi=128):
    w=np.array(w,dtype=float);speed=np.linalg.norm(w);what=w/speed;Omega=abs(Omega)
    e1=np.cross(what,np.array([0.,0.,1.]))
    if np.linalg.norm(e1)<.1:e1=np.cross(what,np.array([1.,0.,0.]))
    e1/=np.linalg.norm(e1);e2=np.cross(what,e1)
    x,wx=leggauss(nk);phi=(np.arange(nphi)+.5)*2*pi/nphi
    transverse=e1[:,None]*np.cos(phi)+e2[:,None]*np.sin(phi)
    maxk=6*kcut;nmax=int(np.floor(maxk*speed/Omega));D=np.zeros((4,4))
    for n in range(-nmax,nmax+1):
        mink=max(kappa*np.exp(-12),abs(n)*Omega/speed)
        if mink>=maxk:continue
        lo=np.log(mink);hi=np.log(maxk);k=np.exp(lo+(x+1)*(hi-lo)/2)
        weight=wx*(hi-lo)/2*k**4/(k*k+kappa*kappa)**2*np.exp(-(k/kcut)**2)
        mu=-n*Omega/(k*speed)
        direction=what[:,None,None]*mu[None,:,None]+np.sqrt(1-mu*mu)[None,:,None]*transverse[:,None,:]
        argument=k[:,None]*vperp/Omega*np.sqrt(np.maximum(0.,1-direction[2]**2))
        W=weight[:,None]*jv(n,argument)**2*(2*pi/nphi)
        g=np.concatenate([direction,np.broadcast_to((n*Omega/k)[None,:,None],(1,nk,nphi))])
        for i in range(4):
            for j in range(i,4):D[i,j]+=np.sum(W*g[i]*g[j])
    return G0/(pi*speed)*(D+np.triu(D,1).T)


def coefficients(vperp,w,Omega,eta=0.,**kw):
    w=np.array(w,dtype=float);h=1e-4*max(1.,np.linalg.norm(w));D=diffusion(vperp,w,Omega,**kw)
    derivative=[]
    for j in range(3):
        shift=np.zeros(3);shift[j]=h
        derivative.append((diffusion(vperp,w+shift,Omega,**kw)-diffusion(vperp,w-shift,Omega,**kw))/(2*h))
    if vperp==0:raise ValueError('Use a finite perpendicular speed for this prototype derivative test.')
    de=1e-4*vperp*vperp/2
    epsilon_derivative=(diffusion(np.sqrt(vperp*vperp+2*de),w,Omega,**kw)-diffusion(np.sqrt(vperp*vperp-2*de),w,Omega,**kw))/(2*de)
    A=.5*(derivative[2][:,2]+epsilon_derivative[:,3])+eta*.5*sum(derivative[j][:,j] for j in range(3))
    return {'diffusion':D.tolist(),'drift':A.tolist(),'w':w.tolist(),'vperp':vperp,'Omega':Omega,'eta':eta}


if __name__=='__main__':
    import argparse,json,time
    from pathlib import Path
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);ap.add_argument('--eta',type=float,default=0.)
    args=ap.parse_args()
    start=time.monotonic();result=coefficients(1.,[0.,.8,.6],1.,eta=args.eta)
    result['seconds']=time.monotonic()-start
    (args.outdir/f'finite-gyro-beam-spectrum-eta{args.eta:g}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
