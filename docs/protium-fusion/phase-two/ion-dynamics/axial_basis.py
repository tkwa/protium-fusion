"""Cylindrical finite-element utilities for the passive axial response model.

The original single-time-origin correlation prototype was rejected because
sampling noise produced nonpassive weak-band estimates. The retained solver
uses orbit-harmonic positive spectral weights in axial_spectral.py.
"""
import numpy as np
EPS0=8.8541878128e-12
C_LIGHT=299792458.

def matrices(nodes,k):
    """Exact FEM weak Poisson stiffness/mass on cylindrical radial elements."""
    n=len(nodes);K=np.zeros((n,n));M=np.zeros_like(K)
    xs,ws=np.polynomial.legendre.leggauss(3)
    for i in range(n-1):
        lo,hi=nodes[i:i+2];d=hi-lo;r=(xs+1)*d/2+lo;w=ws*d/2*2*np.pi*r
        b=np.stack([(hi-r)/d,(r-lo)/d]);der=np.array([-1/d,1/d])
        K[i:i+2,i:i+2]+=np.sum(w)*np.outer(der,der)
        M[i:i+2,i:i+2]+=np.einsum('aq,bq,q->ab',b,b,w)
    return EPS0*(K+k*k*M),K,M

def grid(R,ncore=24,nouter=8):
    return np.r_[np.linspace(0,.008,ncore+1),np.linspace(.008,R,nouter+1)[1:]]

