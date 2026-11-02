"""Resonance-aligned quadrature for the weak helical collision operator.

The physical operator is unchanged from finite_gyro_spectrum.py. Polar nodes
follow the Maxwellian resonance, and azimuth uses Chebyshev integration because
all retained projections depend only on sin(phi). Positive/negative harmonics
are paired by inversion of the wavevector. nmu counts Gaussian-variable nodes;
nphi counts distinct Chebyshev nodes, equivalent to twice as many periodic nodes.
"""
from math import pi,sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import jv
from finite_gyro_spectrum import spectral_log,helical_coefficients as original


def helical_coefficients(vperp,vz,U,vth,Omega,eta,kappa,kcut,G_per_log=1.,nk=96,nmu=48,nphi=48,tail=7.):
    Omega=abs(Omega);w=sqrt(U*U+vz*vz)
    if Omega==0:raise ValueError('Use the analytic straight-orbit operator at zero field.')
    if w/vth<1e-6:
        result=original(vperp,vz,U,vth,Omega,eta,kappa,kcut,G_per_log,nk,nmu,2*nphi)
        result['quadrature_backend']='original at nearly zero relative guiding-center speed'
        return result
    x,wx=leggauss(nk);t0,wt0=leggauss(nmu)
    az=np.cos((np.arange(nphi)+.5)*pi/nphi)
    wy=-U/w;wz=vz/w;kmax=6*kcut;rho=vperp/Omega
    maxarg=kmax*rho
    nmax=1 if vperp==0 else min(int(np.ceil(maxarg+12+4*maxarg**(1/3))),int(np.floor(kmax*(w+tail*vth)/Omega)))
    if nmax>2000:raise ValueError('Harmonic count exceeds prototype limit.')
    D=np.zeros((3,3));div=np.zeros(3)
    for n in range(nmax+1):
        kmin=max(kappa*np.exp(-7),n*Omega/(w+tail*vth))
        if kmin>=kmax:continue
        lo=np.log(kmin);hi=np.log(kmax)
        k=np.exp(lo+(x+1)*(hi-lo)/2)
        wk=wx*(hi-lo)/2*k**4/(k*k+kappa*kappa)**2*np.exp(-(k/kcut)**2)
        harmonic=n*Omega/k
        low=np.maximum(-tail,(-w+harmonic)/vth)
        high=np.minimum(tail,(w+harmonic)/vth)
        t=(low+high)[:,None]/2+(high-low)[:,None]/2*t0[None,:]
        wt=wt0[None,:]*(high-low)[:,None]/2*vth/w*np.exp(-t*t)
        mu=(vth*t-harmonic[:,None])/w
        transverse=np.sqrt(np.maximum(0.,1-mu*mu))
        yy=mu[:,:,None]*wy+transverse[:,:,None]*wz*az[None,None,:]
        zz=mu[:,:,None]*wz-transverse[:,:,None]*wy*az[None,None,:]
        kp2=k[:,None,None]**2*np.maximum(0.,1-zz*zz)
        argument=np.sqrt(kp2)*rho
        if vperp==0:
            J2=np.ones_like(argument) if n==0 else np.zeros_like(argument)
            deriv=(-1. if n==0 else .5)*kp2/Omega**2
        else:
            J=jv(n,argument);J2=J*J
            deriv=(n*J2-argument*J*jv(n+1,argument))/(.5*vperp*vperp)
        multiplicity=1 if n==0 else 2
        weights=multiplicity*wk[:,None,None]*wt[:,:,None]*(2*pi/nphi)
        q=weights*J2
        harmonic=np.broadcast_to(harmonic[:,None,None],yy.shape)
        resonance=np.broadcast_to((vth*t)[:,:,None],yy.shape)
        drift_integrand=weights*(-2*zz*zz*resonance/vth**2*J2+deriv*harmonic)
        g=[yy,zz,harmonic]
        for i in range(3):
            for j in range(i,3):D[i,j]+=np.sum(q*g[i]*g[j])
            div[i]+=np.sum(drift_integrand*g[i])
    prefactor=G_per_log/(pi**1.5*vth)
    D=prefactor*(D+np.triu(D,1).T);div=prefactor*div
    pol=-eta/vth**2*(D@np.array([-U,vz,1.]))
    return {'drift':.5*div+pol,'diffusion':D,'divergence':div,'polarization':pol,
        'harmonic_limit':nmax,'spectral_log':spectral_log(kappa,kcut),
        'quadrature_backend':'relative-velocity-aligned Gaussian variable and Chebyshev azimuth'}
