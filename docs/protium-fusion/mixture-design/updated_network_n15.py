"""Paired NACRE-II N15(p,alpha) evaluation sensitivity; carbon revision fixed.

Published adopted/low/high curves are evaluation limits without an assigned
probability level here. Runtime arrays and uncached kernels avoid global mutation.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from numba import njit
from scipy.interpolate import PchipInterpolator
import pulse_model as pm
import updated_network as carbon
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

@njit(cache=False)
def apply_n15(rates,t,pair,knots,poly,edges,active):
    if active:
        x=np.log(t)
        if x<=knots[0]:lf=edges[0]
        elif x>=knots[-1]:lf=edges[1]
        else:
            j=np.searchsorted(knots,x)-1;dx=x-knots[j]
            lf=((poly[0,j]*dx+poly[1,j])*dx+poly[2,j])*dx+poly[3,j]
        f=np.exp(lf);rates[pair[0]]*=f;rates[pair[1]]*=f
    return rates

@njit(cache=False)
def kernel(temp,rho,y,coeff,rid,pref,dexp,react,stoich,q,qnu,zdel,unit,pairs,knots,poly,edges,active,n15pair,n15knots,n15poly,n15edges):
    nz,nn=y.shape;nr=len(pref)
    dy=np.zeros((nz,nn));qt=np.zeros(nz);qn=np.zeros(nz);qc=np.zeros(nz)
    for i in range(nz):
        t=temp[i]/1e9
        rates=carbon.evaluate_rates(t,coeff,rid,nr,pairs,knots,poly,edges,active)
        rates=apply_n15(rates,t,n15pair,n15knots,n15poly,n15edges,active)
        for j in range(nr):
            f=rates[j]*pref[j]*rho[i]**dexp[j]
            for k in range(3):
                n=react[j,k]
                if n>=0:f*=y[i,n]
            for n in range(nn):dy[i,n]+=stoich[n,j]*f
            qt[i]+=q[j]*f*unit;qn[i]+=qnu[j]*f*unit;qc[i]+=zdel[j]*f
    return dy,qt,qn,qc

def setup(table_path=None,quantile='median',mode='updated',interpolation='pchip',n15_curve='adopted'):
    assert n15_curve in ['adopted','low','high']
    args,meta=carbon.setup(table_path,quantile,mode,interpolation)
    path=HERE/'nacreii-n15.json';d=json.loads(path.read_text());names=[r.fname for r in pm.RATES]
    pair=np.array([names.index('N15_p_to_He4_C12_reaclib'),names.index('C12_He4_to_p_N15_reaclib')],dtype=np.int64)
    rows=d['rows'];x=np.log([r['T9'] for r in rows]);old=np.array([pm.RATES[pair[0]].eval(r['T9']*1e9) for r in rows]);lf=np.log(np.array([r[n15_curve] for r in rows])/old);p=PchipInterpolator(x,lf)
    if interpolation=='linear':p.c[:2]=0.;p.c[2]=np.diff(lf)/np.diff(x);p.c[3]=lf[:-1]
    meta['carbon_kernel_sha256']=meta['kernel_sha256'];meta['kernel_sha256']=sha(__file__)
    meta['nitrogen15_evaluation']={'data_sha256':sha(path),'doi':d['doi'],'table':35,'curve':n15_curve,'temperature_range_GK':[float(np.exp(x[0])),float(np.exp(x[-1]))],'interpretation':d['confidence'],'reverse_policy':'Same forward/reverse factor preserves archived equilibrium ratio; no new Q or species.'}
    meta['revision']+='-plus-nacreii-n15-'+n15_curve
    meta['scope']='Two measured carbon capture updates and one paired N15(p,alpha) evaluation sensitivity; all other physics unchanged.'
    return args+(pair,x,p.c,np.array([lf[0],lf[-1]])),meta

def install(radial,table_path=None,quantile='median',mode='updated',interpolation='pchip',n15_curve='adopted'):
    args,meta=setup(table_path,quantile,mode,interpolation,n15_curve)
    radial.network=lambda temp,rho,y:kernel(temp,rho,y,*args)
    return meta
