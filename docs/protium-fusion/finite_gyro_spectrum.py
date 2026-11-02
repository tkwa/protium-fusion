"""Exploratory static-screened weak helical-orbit fluctuation operator.

Dimensionless or consistently dimensional inputs. Straight Maxwellian ions,
uniform B, and fixed spectral cutoffs. The ion response is retained through
Maxwellian polarization at mass ratio eta. Dynamic dielectric, ion curvature,
electron relativistic effects and nonlinear recollisions are outside this model.

Coordinates for the returned tensor are (integrated Fy/me, vz, vperp**2/2).
The first is an impulse ledger, mapped to GC X by division by signed omega.
The invariant phase-space density is dX dY dvz d(vperp**2/2).
"""
from math import pi,sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv,jvp
from scipy.integrate import quad


def spectral_log(kappa,kcut):
    return quad(lambda x:np.exp(4*x)/(np.exp(2*x)+kappa*kappa)**2*np.exp(-np.exp(2*x)/kcut**2),np.log(kappa)-16,np.log(kcut)+np.log(8),epsabs=1e-11)[0]


def helical_coefficients(vperp,vz,U,vth,Omega,eta,kappa,kcut,G_per_log=1.,nk=80,nmu=48,nphi=48):
    Omega=abs(Omega)
    if Omega==0:raise ValueError('Use the analytic straight-orbit operator at zero field.')
    x,wx=leggauss(nk);mu,wm=leggauss(nmu)
    lo=np.log(kappa)-12;hi=np.log(kcut)+np.log(6.)
    k=np.exp(lo+(x+1)*(hi-lo)/2)
    wk=wx*(hi-lo)/2*k**4/(k*k+kappa*kappa)**2*np.exp(-(k/kcut)**2)
    phi=(np.arange(nphi)+.5)*2*pi/nphi
    ny=np.sqrt(1-mu*mu)[:,None]*np.sin(phi)
    # k, mu, phi axes; each harmonic is accumulated independently.
    yy=np.broadcast_to(ny[None,:,:],(nk,nmu,nphi))
    zz=np.broadcast_to(mu[None,:,None],yy.shape)
    base=zz*vz-yy*U
    weights=wk[:,None,None]*wm[None,:,None]*(2*pi/nphi)
    argument=k[:,None]*np.sqrt(1-mu*mu)[None,:]*vperp/Omega
    nmax=1 if vperp==0 else int(np.ceil(argument.max()+12+4*argument.max()**(1/3)))
    if nmax>2000:raise ValueError('Harmonic count exceeds prototype limit; use close-encounter matching.')
    D=np.zeros((3,3));div=np.zeros(3);pol=np.zeros(3)
    for n in range(-nmax,nmax+1):
        if vperp==0:
            J2=np.ones_like(argument) if n==0 else np.zeros_like(argument)
            deriv=(-1. if n==0 else .5)*k[:,None]**2*(1-mu*mu)[None,:]/Omega**2
        else:
            J=jv(n,argument);J2=J*J
            deriv=J*jvp(n,argument)*argument/(.5*vperp*vperp)
        harmonic=np.broadcast_to((n*Omega/k)[:,None,None],yy.shape)
        resonance=base+harmonic
        gaussian=np.exp(-(resonance/vth)**2)
        common=weights*gaussian
        g=[yy,zz,harmonic]
        q=common*J2[:,:,None]
        dv=q*(-2*zz*resonance/vth**2)
        de=common*deriv[:,:,None]
        for i in range(3):
            for j in range(i,3):
                D[i,j]+=np.sum(q*g[i]*g[j])
            div[i]+=np.sum(dv*g[i]*zz+de*g[i]*harmonic)
            pol[i]+=-eta/vth**2*np.sum(q*g[i]*resonance)
    D=D+np.triu(D,1).T
    prefactor=G_per_log/(pi**1.5*vth)
    return {'drift':prefactor*(.5*div+pol),'diffusion':prefactor*D,
            'divergence':prefactor*div,'polarization':prefactor*pol,
            'harmonic_limit':nmax,'spectral_log':spectral_log(kappa,kcut)}


def unmagnetized_gyro_average(vperp,vz,U,vth,eta,logarithm,G_per_log=1.,nphi=512):
    from guiding_center_collisions import force_tensor_derivatives
    D=np.zeros((3,3));A=np.zeros(3)
    for phi in (np.arange(nphi)+.5)*2*pi/nphi:
        v=np.array([vperp*np.cos(phi),vperp*np.sin(phi),vz]);w=v-np.array([0.,U,0.]);speed=np.linalg.norm(w)
        dp,dt,dpp,dtp=force_tensor_derivatives(speed,G_per_log*logarithm,vth*vth)
        forceD=dt*np.eye(3)+(dp-dt)*np.outer(w,w)/speed**2 if speed else dt*np.eye(3)
        forceA=(.5*(dpp+2*(dp-dt)/speed)-eta/vth**2*dp*speed)*w/speed if speed else np.zeros(3)
        T=np.array([[0.,1.,0.],[0.,0.,1.],[v[0],v[1],0.]])
        D+=T@forceD@T.T/nphi
        A+=(T@forceA+np.array([0.,0.,.5*(forceD[0,0]+forceD[1,1])]))/nphi
    return {'drift':A,'diffusion':D}
