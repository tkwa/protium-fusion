"""Independent finite-gyro weak binary encounters for a smooth Yukawa field.

The potential Fourier amplitude is -4*pi/(k²+kappa²)*exp(-k²/(2*kcut²)).
The ion is unmagnetized on encounter times, with first-order recoil retained.
Electron trajectories retain full Lorentz rotation at both perturbative orders.
"""
from math import pi,sqrt
import numpy as np
from scipy.special import erfc,erfcx
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss


def radial_field(r,kappa=1.,kcut=4.):
    """Return V, h, h'/r for force F=-h*r and gradient -h I-(h'/r) rr."""
    r=np.asarray(r,dtype=float);shape=r.shape;rr=r.ravel()
    a=kappa/(sqrt(2)*kcut);b=kcut/sqrt(2)
    V=np.empty_like(rr);h=V.copy();hp=V.copy()
    small=rr*kcut<.08
    if np.any(small):
        # H(r)=sum c_n r^(2n+1), V=-H/r. Recurrence from its radial ODE.
        c=[-kappa*erfcx(a)+2*b/sqrt(pi)]
        term=1.
        for n in range(7):
            c.append((kappa*kappa*c[-1]-4*b**3/sqrt(pi)*term)/((2*n+3)*(2*n+2)))
            term*=-b*b/(n+1)
        x=rr[small]**2
        V[small]=-sum(c[n]*x**n for n in range(len(c)))
        h[small]=-sum(2*n*c[n]*x**(n-1) for n in range(1,len(c)))
        hp[small]=-sum(4*n*(n-1)*c[n]*x**(n-2) for n in range(2,len(c)))
    if np.any(~small):
        t=rr[~small];arg=a-b*t;gauss=np.exp(-(b*t)**2)
        plus=gauss*erfcx(a+b*t)
        minus=np.empty_like(t);positive=arg>=0
        minus[positive]=gauss[positive]*erfcx(arg[positive])
        minus[~positive]=np.exp(a*a-kappa*t[~positive])*erfc(arg[~positive])
        H=.5*(minus-plus)
        Hd=-kappa*.5*(minus+plus)+2*b/sqrt(pi)*gauss
        Hdd=kappa*kappa*H-4*b**3*t/sqrt(pi)*gauss
        N=H-t*Hd
        V[~small]=-H/t;h[~small]=N/t**3;hp[~small]=-Hdd/t**3-3*N/t**5
    return V.reshape(shape),h.reshape(shape),hp.reshape(shape)


def ensemble(vperp=1.,w=np.array([0.,.8,.6]),omega=-1.,eta=0.,nb=40,nphi=8,ngyro=8,rtol=3e-8,kappa=1.,kcut=4.):
    import time
    start=time.monotonic();w=np.array(w,dtype=float);u=np.linalg.norm(w);what=w/u
    e1=np.cross(what,np.array([0.,0.,1.]))
    if np.linalg.norm(e1)<.1:e1=np.cross(what,np.array([1.,0.,0.]))
    e1/=np.linalg.norm(e1);e2=np.cross(what,e1)
    nodes,weights=leggauss(nb);rho=vperp/abs(omega)
    lo=np.log(1e-5/kcut);hi=np.log(24/kappa+rho)
    b=np.exp(lo+(nodes+1)*(hi-lo)/2);wb=weights*(hi-lo)/2*b*b
    phi=(np.arange(nphi)+.5)*2*pi/nphi
    gyro=(np.arange(ngyro)+.5)*2*pi/ngyro
    bhat=e1[:,None]*np.cos(phi)+e2[:,None]*np.sin(phi)
    bvec=np.broadcast_to(b[None,:,None,None]*bhat[:,None,:,None],(3,nb,nphi,ngyro)).reshape(3,-1)
    phase=np.broadcast_to(gyro[None,None,:],(nb,nphi,ngyro)).ravel()
    weight=np.repeat(wb,nphi*ngyro)*u*2*pi/nphi/ngyro
    n=len(weight)
    T=(28/kappa+rho)/u
    def rhs(t,y):
        y=y.reshape(15,n);dr=y[:3];v1=y[3:6];I1=y[6:9];v2=y[9:12]
        theta=phase-omega*t
        gyr=np.array([-vperp/omega*np.sin(theta),vperp/omega*np.cos(theta),np.zeros(n)])
        r=bvec+w[:,None]*t+gyr
        _,h,hp=radial_field(np.sqrt(np.sum(r*r,axis=0)),kappa,kcut)
        F=-h*r;F2=-h*dr-hp*r*np.sum(r*dr,axis=0)
        a1=F.copy();a2=F2.copy()
        a1[0]+=omega*v1[1];a1[1]-=omega*v1[0]
        a2[0]+=omega*v2[1];a2[1]-=omega*v2[0]
        return np.vstack([v1+eta*I1,a1,F,a2,F2]).ravel()
    sol=solve_ivp(rhs,[-T,T],np.zeros(15*n),method='DOP853',rtol=rtol,atol=rtol*.003,t_eval=[T])
    final=sol.y[:,-1].reshape(15,n);v1=final[3:6];I1=final[6:9];v2=final[9:12];I2=final[12:15]
    theta=phase-omega*T;v0=np.array([vperp*np.cos(theta),vperp*np.sin(theta)])
    P1=np.sum(v0*v1[:2],axis=0)
    P2=np.sum(v0*v2[:2],axis=0)+.5*np.sum(v1[:2]**2,axis=0)
    q1=np.vstack([I1,P1]);q2=np.vstack([I2,P2])
    D=(q1*weight)@q1.T;A=q2@weight
    closure1=w@I1+P1
    closure2=w@I2+P2+.5*(I1[2]**2+eta*np.sum(I1*I1,axis=0))
    norm1=np.max(np.abs(w@I1))+np.max(np.abs(P1))
    norm2=np.max(np.abs(w@I2))+np.max(np.abs(P2))+.5*np.max(I1[2]**2+eta*np.sum(I1*I1,axis=0))
    return {'vperp':vperp,'w':w.tolist(),'omega':omega,'eta':eta,'nb':nb,'nphi':nphi,'ngyro':ngyro,'rtol':rtol,
        'diffusion':D.tolist(),'drift':A.tolist(),'first_order_mean':(q1@weight).tolist(),
        'first_order_energy_error':float(np.max(np.abs(closure1))/max(norm1,1e-100)),
        'second_order_energy_error':float(np.max(np.abs(closure2))/max(norm2,1e-100)),
        'solver_success':bool(sol.success),'nfev':sol.nfev,'seconds':time.monotonic()-start}


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent)
    ap.add_argument('--nb',type=int,default=40);ap.add_argument('--nphi',type=int,default=8);ap.add_argument('--ngyro',type=int,default=8)
    ap.add_argument('--eta',type=float,default=0.);ap.add_argument('--omega',type=float,default=-1.)
    args=ap.parse_args()
    result=ensemble(nb=args.nb,nphi=args.nphi,ngyro=args.ngyro,eta=args.eta,omega=args.omega)
    (args.outdir/f'finite-gyro-encounters-{args.nb}-{args.nphi}-{args.ngyro}-eta{args.eta:g}-omega{args.omega:g}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
