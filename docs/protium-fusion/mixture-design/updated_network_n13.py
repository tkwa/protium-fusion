"""Distinct Li2006 Eq17 N13 capture sensitivity, combined with carbon revision.

The source fit differs strongly from JINA's later lg06 version5 above1 GK.
This module does not assert the library is wrong. It isolates the consequence
of using the cited authors' explicit 0.01–10 GK fit. Reverse rates receive the
same new/old factor. Source fit residual is not experimental uncertainty.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from numba import njit
import pulse_model as pm
import updated_network as carbon
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

@njit(cache=False)
def kernel(temp,rho,y,coeff,rid,pref,dexp,react,stoich,q,qnu,zdel,unit,pairs,knots,poly,edges,active,n13pair,n13coeff):
    nz,nn=y.shape;nr=len(pref)
    dy=np.zeros((nz,nn));qt=np.zeros(nz);qn=np.zeros(nz);qc=np.zeros(nz)
    for i in range(nz):
        t=temp[i]/1e9
        rates=carbon.evaluate_rates(t,coeff,rid,nr,pairs,knots,poly,edges,active)
        if active:
            # Clip evaluation to the stated fit interval, retaining the original
            # rate's endpoint-scaled temperature behavior outside it.
            tn=min(max(t,.01),10.)
            tf=np.array([1.,1/tn,tn**(-1/3),tn**(1/3),tn,tn**(5/3),np.log(tn)])
            newrate=0.;oldrate=0.
            for row in n13coeff:newrate+=np.exp(row@tf)
            for s in range(len(coeff)):
                if rid[s]==n13pair[0]:oldrate+=np.exp(coeff[s]@tf)
            factor=newrate/oldrate
            rates[n13pair[0]]*=factor;rates[n13pair[1]]*=factor
        for j in range(nr):
            f=rates[j]*pref[j]*rho[i]**dexp[j]
            for k in range(3):
                n=react[j,k]
                if n>=0:f*=y[i,n]
            for n in range(nn):dy[i,n]+=stoich[n,j]*f
            qt[i]+=q[j]*f*unit;qn[i]+=qnu[j]*f*unit;qc[i]+=zdel[j]*f
    return dy,qt,qn,qc

def setup(table_path=None,quantile='median',mode='updated',interpolation='pchip'):
    args,meta=carbon.setup(table_path,quantile,mode,interpolation)
    data=json.loads((Path(table_path) if table_path else carbon.HERE/'updated-carbon-rates.json').read_text())
    grids=[np.array(data['tables'][a]['rows'])[:,0] for a,b in carbon.PAIRS]
    assert all(np.array_equal(grids[0],g) for g in grids[1:])
    d=json.loads((HERE/'n13-source-fit.json').read_text());names=[r.fname for r in pm.RATES]
    pair=np.array([names.index(d['forward']),names.index(d['reverse'])],dtype=np.int64)
    co=np.array(d['coefficients'])
    meta['carbon_kernel_sha256']=meta['kernel_sha256'];meta['kernel_sha256']=sha(__file__)
    meta['nitrogen_source_fit']={'data_sha256':sha(HERE/'n13-source-fit.json'),'doi':d['doi'],'source_equation':17,
      'temperature_range_GK':[.01,10.],'interpretation':'Separate source-fit sensitivity; discrepancy with JINA lg06 version5 is unresolved. Published<5% fit residual is not a rate uncertainty.',
      'uncertainty_selector_scope':'Quantile selector applies only to the two carbon captures; N13 source fit held fixed.'}
    meta['revision']+='-plus-li2006-eq17-sensitivity'
    meta['scope']='Three capture pairs changed; N13 source-fit discrepancy separate from new measured carbon normalization.'
    return args+(pair,co),meta

def install(radial,table_path=None,quantile='median',mode='updated',interpolation='pchip'):
    args,meta=setup(table_path,quantile,mode,interpolation)
    radial.network=lambda temp,rho,y:kernel(temp,rho,y,*args)
    return meta
