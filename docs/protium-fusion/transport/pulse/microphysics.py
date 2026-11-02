"""Diagnostic microphysics; this module does not evolve a pulse.

All public quantities use SI. Scattering, true absorption, Coulomb corrections,
and conductivity are deliberately separate. Models outside their domain are
flagged by the caller; no extrapolated coefficient is a certified bound.
"""
import numpy as np
from scipy.special import expit
from numpy.polynomial.legendre import leggauss

KB=1.380649e-23; H=6.62607015e-34; HB=H/(2*np.pi)
C=299792458.; ME=9.1093837139e-31; MU=1.66053906660e-27
E=1.602176634e-19; EPS=8.8541878128e-12
EC=E*E/(4*np.pi*EPS); SIGT=6.6524587051e-29
ARAD=7.565733250033928e-16

def fermions(T, ne, order=128, iterations=55):
    """Ideal relativistic Fermi gas, charge-neutrality with thermal pairs.

    Returns eta=(mu-mc2)/kT, electron density, positron density. Fixed Gaussian
    quadrature in momentum; verify separately against refined quadrature and
    analytic limits. This is a diagnostic, independent of the hydro EOS table.
    """
    T,ne=np.broadcast_arrays(T,ne); shape=T.shape
    theta=(KB*T/(ME*C*C)).ravel(); net=ne.ravel()
    xf=HB*(3*np.pi*np.pi*net)**(1/3)/(ME*C)
    tf=(np.sqrt(1+xf*xf)-1)/theta
    low=-1/theta; high=tf+20
    # The upper limit remains fixed during each neutrality root solve.
    umax=np.sqrt((1+theta*(high+45))**2-1)
    nodes,weights=leggauss(order)
    pref=(ME*C/HB)**3/np.pi**2
    eta_out=np.empty_like(theta); nm_out=np.empty_like(theta); np_out=np.empty_like(theta)
    for start in range(0,len(theta),2048):
        sl=slice(start,start+2048)
        u=.5*umax[sl,None]*(nodes+1)
        wt=.5*umax[sl,None]*weights*u*u
        kinetic=(np.sqrt(1+u*u)-1)/theta[sl,None]
        lo=low[sl].copy();hi=high[sl].copy()
        for _ in range(iterations):
            mid=.5*(lo+hi)
            nm=pref*np.sum(wt*expit(mid[:,None]-kinetic),axis=1)
            np_=pref*np.sum(wt*expit((-2/theta[sl]-mid)[:,None]-kinetic),axis=1)
            less=nm-np_<net[sl]
            lo=np.where(less,mid,lo);hi=np.where(less,hi,mid)
        eta_out[sl]=mid;nm_out[sl]=nm;np_out[sl]=np_
    return tuple(x.reshape(shape) for x in [eta_out,nm_out,np_out])

def poutanen(T,eta,leptons,rho):
    """Poutanen 2017 eq36/Table1, 2--300 keV, -10<=eta<=7.

    Pair scattering is added by using the total lepton density. The fit was
    reported for electron scattering; positron fraction is returned elsewhere
    so pair-rich application is visible. Values outside fit domain are NaN.
    """
    kev=KB*T/E/1000
    x=np.exp(.682*eta-.0454*eta*eta)
    f1=1+.24*x+.0043*x*x
    f2=1+.050*x-.0067*x*x
    f3=1-.037*x+.0031*x*x
    valid=(kev>=2)&(kev<=300)&(eta>=-10)&(eta<=7)&(f2>0)
    lam=f1*(1+(kev/np.maximum(43.3*f2,1e-100))**(.885*f3))
    return np.where(valid,SIGT*leptons/rho/lam,np.nan),valid

def freefree(T,ne,z2ni,x=5.,gaunt=1.2):
    """Classical nonrelativistic free--free emissivity and absorption.

    alpha_nu=3.692e8 T^-1/2 ne sum(Z²ni) nu^-3(1-exp(-x))g in cgs.
    Returned absorption [m^-1] includes stimulated emission; emissivity [W/m3]
    is frequency-integrated. Degeneracy, relativity, collective dielectric and
    bound-state corrections are not hidden in a fitted gray parameter.
    """
    freq=x*KB*T/H
    alpha=3.692e-2*T**(-.5)*ne*z2ni*freq**(-3)*(-np.expm1(-x))*gaunt
    emission=1.426e-40*np.sqrt(T)*ne*z2ni*gaunt
    alpha_planck=emission/(C*ARAD*T**4)
    return alpha,emission,alpha_planck

def ocp_internal(gamma):
    """Gamma*f'(Gamma), Chugunov & DeWitt 2009 eq24 / PC OCP fit.

    Analytic derivative; linear mixing caller includes only ion-ion rigid-
    background free energy. Electron polarization and mixing corrections absent.
    """
    a1=-.907;a2=.62954;a3=-np.sqrt(3)/2-a1/np.sqrt(a2)
    b1=.00456;b2=211.6;b3=-1e-4;b4=.00462
    return (a1*gamma**1.5/np.sqrt(a2+gamma)
            +a3*gamma**1.5/(1+gamma)
            +b1*gamma**2/(b2+gamma)+b3*gamma**2/(b4+gamma**2))

def coulomb(T,rho,Y,Z):
    """Return ion-ion energy J/kg, pressure Pa and dominant coupling."""
    ne=rho/MU*np.sum(Y*Z,axis=-1)
    ge=EC*(4*np.pi*ne/3)**(1/3)/(KB*T)
    gamma=ge[...,None]*Z**(5/3)
    u=KB*T/MU*np.sum(Y*ocp_internal(gamma),axis=-1)
    dominant=np.argmax(Y*Z,axis=-1)
    gd=np.take_along_axis(gamma,dominant[...,None],axis=-1)[...,0]
    return u,rho*u/3,gd,ge

def conductivity_envelope(T,ne,z2ni,lnlambda=1.):
    """Max of classical Lorentz-gas and degenerate e-ion conductivity.

    lnLambda=1 is a high-conductivity sensitivity choice, not a proven lower
    bound for the effective collision integral. No e-e resistance is included.
    There is no claim that max is a physically correct interpolation.
    """
    zeff=z2ni/ne
    # Spitzer-Harm Lorentz gas, zero-current epsilon=0.4, delta_T=1.
    classical=20*(2/np.pi)**1.5*.4*KB*(KB*T)**2.5/(EC**2*zeff*np.sqrt(ME)*lnlambda)
    pf=HB*(3*np.pi**2*ne)**(1/3)
    mf=np.sqrt(ME*ME+(pf/C)**2);vf=pf/mf
    nuei=4*np.pi*EC**2*z2ni*lnlambda/(pf*pf*vf)
    degenerate=np.pi**2*KB*KB*T*ne/(3*mf*nuei)
    return np.maximum(classical,degenerate),classical,degenerate

def weighted_quantile(values,weights,probs=(.01,.5,.99)):
    v=np.asarray(values).ravel();w=np.asarray(weights).ravel()
    keep=np.isfinite(v)&np.isfinite(w)&(w>0)
    v=v[keep];w=w[keep]
    if not len(v):return [None]*len(probs)
    ix=np.argsort(v);v=v[ix];w=w[ix]
    q=(np.cumsum(w)-.5*w)/np.sum(w)
    return np.interp(probs,q,v).tolist()
