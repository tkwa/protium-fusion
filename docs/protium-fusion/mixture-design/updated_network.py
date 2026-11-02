"""Runtime carbon-rate revision; frozen 27-species stoichiometry and hydro unchanged.

install(radial, table_path=None, quantile='median', mode='updated') returns
provenance and replaces only this process's radial.network callable. Every
numerical array is a runtime argument to an uncached Numba kernel. No cached
module constants are mutated. Both directions of each capture share a factor.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from numba import njit
from scipy.interpolate import PchipInterpolator
import pulse_model as pm
HERE=Path(__file__).resolve().parent
PAIRS=[('C12_p_to_N13_reaclib','N13_to_p_C12_reaclib'),('C13_p_to_N14_reaclib','N14_to_p_C13_reaclib')]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

@njit(cache=False)
def evaluate_rates(t,coeff,rid,nr,pairs,knots,poly,edges,active):
    tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
    rates=np.zeros(nr)
    for s in range(len(coeff)):
        exponent=0.
        for k in range(7):exponent+=coeff[s,k]*tf[k]
        rates[rid[s]]+=np.exp(exponent)
    if active:
        x=np.log(t)
        for p in range(len(pairs)):
            if x<=knots[0]:lf=edges[p,0]
            elif x>=knots[-1]:lf=edges[p,1]
            else:
                j=np.searchsorted(knots,x)-1;dx=x-knots[j]
                lf=((poly[p,0,j]*dx+poly[p,1,j])*dx+poly[p,2,j])*dx+poly[p,3,j]
            fac=np.exp(lf)
            rates[pairs[p,0]]*=fac;rates[pairs[p,1]]*=fac
    return rates

@njit(cache=False)
def kernel(temp,rho,y,coeff,rid,pref,dexp,react,stoich,q,qnu,zdel,unit,pairs,knots,poly,edges,active):
    nz,nn=y.shape;nr=len(pref)
    dy=np.zeros((nz,nn));qt=np.zeros(nz);qn=np.zeros(nz);qc=np.zeros(nz)
    for i in range(nz):
        rates=evaluate_rates(temp[i]/1e9,coeff,rid,nr,pairs,knots,poly,edges,active)
        for j in range(nr):
            f=rates[j]*pref[j]*rho[i]**dexp[j]
            for k in range(3):
                n=react[j,k]
                if n>=0:f*=y[i,n]
            for n in range(nn):dy[i,n]+=stoich[n,j]*f
            qt[i]+=q[j]*f*unit;qn[i]+=qnu[j]*f*unit;qc[i]+=zdel[j]*f
    return dy,qt,qn,qc

def setup(table_path=None,quantile='median',mode='updated',interpolation='pchip'):
    path=Path(table_path) if table_path else HERE/'updated-carbon-rates.json'
    d=json.loads(path.read_text());column=d['columns'].index(quantile)
    names=[r.fname for r in pm.RATES];pairs=np.array([[names.index(a),names.index(b)] for a,b in PAIRS],dtype=np.int64)
    polynomials=[];edges=[];ratios=[];first_grid=None
    for k,(a,b) in enumerate(PAIRS):
        rows=np.array(d['tables'][a]['rows']);x=np.log(rows[:,0]);old=np.array([pm.RATES[pairs[k,0]].eval(t*1e9) for t in rows[:,0]])
        if first_grid is None:first_grid=rows[:,0].copy()
        else:assert np.array_equal(first_grid,rows[:,0]), "Carbon table temperature grids differ"
        lf=np.log(rows[:,column]/old);p=PchipInterpolator(x,lf)
        if interpolation=='linear':
            p.c[:2]=0.;p.c[2]=np.diff(lf)/np.diff(x);p.c[3]=lf[:-1]
        polynomials.append(p.c);edges.append([lf[0],lf[-1]]);ratios.append(np.exp(lf))
    args=(pm.COEFF,pm.RATE_INDEX,pm.PREF,pm.DEXP,pm.REACT,pm.STOICH,pm.Q,pm.QNU,pm.ZZ@pm.STOICH,pm.MEV/pm.MU*1e4,pairs,x,np.array(polynomials),np.array(edges),mode!='identity')
    meta={'revision':d['revision'],'mode':mode,'quantile':quantile,'table_sha256':sha(path),'kernel_sha256':sha(__file__),
          'source_doi':d['paper_doi'],'temperature_range_GK':[float(np.exp(x[0])),float(np.exp(x[-1]))],
          'interpolation':interpolation+' in log(new/old) versus log(T9); frozen REACLIB shape supplies within-cell temperature dependence',
          'outside_table':'Constant endpoint new/old ratio; relevant trajectories remain within table during burn.',
          'reverse_policy':'Identical multiplicative factor for each forward/reverse pair, preserving the original library equilibrium ratio and Q.',
          'pairs':[{'forward':a,'reverse':b,'indices':pairs[k].tolist()} for k,(a,b) in enumerate(PAIRS)],
          'scope':'Two carbon proton-capture pairs only; 27 species, other 190 rates, EOS, opacity, hydro and weak-energy convention unchanged.'}
    return args,meta

def install(radial,table_path=None,quantile='median',mode='updated',interpolation='pchip'):
    args,meta=setup(table_path,quantile,mode,interpolation)
    radial.network=lambda temp,rho,y:kernel(temp,rho,y,*args)
    return meta
