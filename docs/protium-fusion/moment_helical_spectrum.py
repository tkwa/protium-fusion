"""Azimuthal-moment reduction of the weak helical collision operator.

Retains the original magnetic-axis polar quadrature. The azimuth integral is
reduced to three Gaussian moments, and Bessel values are shared between ion
species at the same electron state and wavevector grid. No physical term is
removed. Positive/negative cyclotron harmonics are paired by k -> -k symmetry.
"""
from functools import lru_cache
from math import pi
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv
from finite_gyro_spectrum import spectral_log


@lru_cache(maxsize=2)
def wave_grid(vperp,Omega,kappa,kcut,nk,nmu):
    x,wx=leggauss(nk);mu,wm=leggauss(nmu)
    lo=np.log(kappa)-12;hi=np.log(kcut)+np.log(6.)
    k=np.exp(lo+(x+1)*(hi-lo)/2)
    wk=wx*(hi-lo)/2*k**4/(k*k+kappa*kappa)**2*np.exp(-(k/kcut)**2)
    argument=k[:,None]*np.sqrt(1-mu*mu)[None,:]*vperp/Omega
    nmax=1 if vperp==0 else int(np.ceil(argument.max()+12+4*argument.max()**(1/3)))
    if nmax>2000:raise ValueError('Harmonic count exceeds prototype limit.')
    B=None if vperp==0 else jv(np.arange(nmax+2)[:,None,None],argument[None,:,:])
    return k,mu,wk[:,None]*wm[None,:],argument,B,nmax


def helical_coefficients(vperp,vz,U,vth,Omega,eta,kappa,kcut,G_per_log=1.,nk=96,nmu=96,nphi=24):
    Omega=abs(Omega)
    if Omega==0:raise ValueError('Use the analytic straight-orbit operator at zero field.')
    k,mu,weight,argument,B,nmax=wave_grid(vperp,Omega,kappa,kcut,nk,nmu)
    mu=mu[None,:];sine=np.sqrt(1-mu*mu)
    az=np.cos((np.arange(nphi)+.5)*pi/nphi)
    D=np.zeros((3,3));div=np.zeros(3)
    for n in range(nmax+1):
        harmonic=(n*Omega/k)[:,None];a=mu*vz+harmonic;b=U*sine
        gaussian=np.exp(-((a[:,:,None]-b[:,:,None]*az[None,None,:])/vth)**2)
        I0=np.sum(gaussian,axis=2)*(2*pi/nphi)
        I1=np.sum(gaussian*az[None,None,:],axis=2)*(2*pi/nphi)
        I2=np.sum(gaussian*az[None,None,:]**2,axis=2)*(2*pi/nphi)
        if vperp==0:
            J2=np.ones_like(argument) if n==0 else np.zeros_like(argument)
            derivative=(-1. if n==0 else .5)*k[:,None]**2*sine*sine/Omega**2
        else:
            J=B[n];J2=J*J
            derivative=(n*J2-argument*J*B[n+1])/(.5*vperp*vperp)
        W=weight*(1 if n==0 else 2);Q=W*J2;E=W*derivative
        D[0,0]+=np.sum(Q*sine*sine*I2)
        D[0,1]+=np.sum(Q*sine*mu*I1)
        D[0,2]+=np.sum(Q*sine*harmonic*I1)
        D[1,1]+=np.sum(Q*mu*mu*I0)
        D[1,2]+=np.sum(Q*mu*harmonic*I0)
        D[2,2]+=np.sum(Q*harmonic*harmonic*I0)
        R0=a*I0-b*I1;R1=a*I1-b*I2
        div[0]+=np.sum(-2*Q*mu*mu*sine/vth**2*R1+E*sine*harmonic*I1)
        div[1]+=np.sum(-2*Q*mu**3/vth**2*R0+E*mu*harmonic*I0)
        div[2]+=np.sum(-2*Q*mu*mu*harmonic/vth**2*R0+E*harmonic*harmonic*I0)
    prefactor=G_per_log/(pi**1.5*vth)
    D=prefactor*(D+np.triu(D,1).T);div=prefactor*div
    pol=-eta/vth**2*(D@np.array([-U,vz,1.]))
    return {'drift':.5*div+pol,'diffusion':D,'divergence':div,'polarization':pol,
            'harmonic_limit':nmax,'spectral_log':spectral_log(kappa,kcut),
            'quadrature_backend':'magnetic-axis polar nodes; three Chebyshev azimuth moments; shared Bessel grid'}
