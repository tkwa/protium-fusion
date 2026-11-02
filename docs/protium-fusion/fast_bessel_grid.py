"""Integer J_n(x) table with stable forward recurrence and a tiny-tail bound.

Forward recurrence is used only for n<x; the transition/evanescent region
uses scipy.special.jv. This follows the stability distinction in DLMF10.74(iv).
For n>x, the generating-function contour gives the absolute bound
exp(sqrt(n*n-x*x)-n*acosh(n/x)); values with bound below exp(-30) are omitted.
"""
import numpy as np
from scipy.special import j0,j1,jv


def _next_order(n,x,previous,current,log_tolerance):
    result=np.zeros_like(x)
    forward=x>n
    result[forward]=2*(n-1)/x[forward]*current[forward]-previous[forward]
    other=(x>0)&(~forward)
    xs=x[other]
    logarithm=np.sqrt(np.maximum(0.,n*n-xs*xs))-n*np.arccosh(n/xs)
    keep=np.zeros_like(other);keep[other]=logarithm>=log_tolerance
    result[keep]=jv(n,x[keep])
    return result


def integer_bessel_grid(argument,nmax,log_tolerance=-30.):
    x=np.asarray(argument,dtype=float)
    if np.any(x<0):raise ValueError('This table expects nonnegative real arguments.')
    out=np.zeros((nmax+1,*x.shape));out[0]=j0(x)
    if nmax==0:return out
    out[1]=j1(x)
    for n in range(2,nmax+1):
        out[n]=_next_order(n,x,out[n-2],out[n-1],log_tolerance)
    return out


def integer_bessel_pairs(argument,nmax,log_tolerance=-30.):
    """Yield (n,J_n,J_(n+1)) while retaining only two previous orders."""
    x=np.asarray(argument,dtype=float)
    if np.any(x<0):raise ValueError('This table expects nonnegative real arguments.')
    current=j0(x);following=j1(x)
    for n in range(nmax+1):
        yield n,current,following
        if n<nmax:
            newest=_next_order(n+2,x,current,following,log_tolerance)
            current,following=following,newest
