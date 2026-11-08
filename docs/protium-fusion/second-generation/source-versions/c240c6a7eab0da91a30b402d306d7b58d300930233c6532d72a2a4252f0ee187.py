#!/usr/bin/env python3
"""Extended-network postprocessing on prescribed trajectories, with no hydro feedback.

Tracer calculations follow one heavy nucleus in the saved n/p/alpha bath.
Bulk-mixture tests instead evolve all 66 abundances at prescribed T,rho.
Neither calculation establishes a new self-consistent pulse.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

P=Path(__file__).resolve().parent
D=json.loads((P/'rates.json').read_text())
NAMES=D['species']; A=np.array(D['A']); Z=np.array(D['Z']); MASS=np.array(D['mass_MeV'])
NN=len(NAMES); NR=len(D['rates'])
S=np.zeros((NN,NR)); REACT=np.full((NR,max(len(r['react']) for r in D['rates'])),-1,dtype=np.int64)
CO=[]; RI=[]
for j,r in enumerate(D['rates']):
    for k,i in enumerate(r['react']): REACT[j,k]=i; S[i,j]-=1
    for i in r['product']: S[i,j]+=1
    for c in r['coefficients']: CO.append(c); RI.append(j)
CO=np.array(CO); RI=np.array(RI,dtype=np.int64)
PREF=np.array([r['prefactor'] for r in D['rates']]); DE=np.array([r['density_exponent'] for r in D['rates']])
YE=np.array([r['ye_weighting'] for r in D['rates']]); Q=-MASS@S
HEAVY=np.flatnonzero(A>=10); HI={i:k for k,i in enumerate(HEAVY)}; NH=len(HEAVY)
WATCH=['N15_p_to_O16_reaclib','O16_p_to_F17_reaclib','O17_p_to_He4_N14_reaclib',
       'O17_p_to_F18_reaclib','F17_p_to_Ne18_reaclib','F18_p_to_He4_O15_reaclib',
       'O18_p_to_He4_N15_reaclib','O18_p_to_F19_reaclib','F19_p_to_He4_O16_reaclib',
       'F19_p_to_Ne20_reaclib','C12_He4_to_O16_reaclib','O15_He4_to_Ne19_reaclib',
       'O14_He4_to_p_F17_reaclib','Ne18_He4_to_p_Na21_reaclib','O16_He4_to_Ne20_reaclib',
       'F18_p_to_Ne19_reaclib','N13_He4_to_p_O16_reaclib']
WI=np.array([next(j for j,r in enumerate(D['rates']) if r['name']==n) for n in WATCH])
CARBON=json.loads((P.parent/'mixture-design/updated-carbon-rates.json').read_text())
TSRC=np.full(NR,-1,dtype=np.int64); TDST=np.full(NR,-1,dtype=np.int64)
for j,r in enumerate(D['rates']):
    hr=[i for i in r['react'] if i in HI]; hp=[i for i in r['product'] if i in HI]
    if len(hr)==1 and len(hp)==1 and all(A[i]<=4 or i==hr[0] for i in r['react']):
        TSRC[j]=HI[hr[0]]; TDST[j]=HI[hp[0]]

@njit(cache=True)
def rate_values(T):
    t=max(T/1e9,1e-5)
    tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
    out=np.zeros(NR)
    for k in range(len(CO)):
        x=0.
        for a in range(7): x+=CO[k,a]*tf[a]
        out[RI[k]]+=np.exp(min(690.,max(-745.,x)))
    return out

@njit(cache=True)
def flows(y,rho,rr):
    out=rr*PREF*rho**DE
    for j in range(NR):
        for k in range(REACT.shape[1]):
            i=REACT[j,k]
            if i>=0: out[j]*=y[i]
        if YE[j]: out[j]*=np.dot(Z,y)
    return out

@njit(cache=True)
def interpolate(x, grid, vals):
    k=max(0,min(len(grid)-2,np.searchsorted(grid,x)-1))
    w=max(0.,min(1.,(x-grid[k])/(grid[k+1]-grid[k])))
    return vals[k]*(1-w)+vals[k+1]*w

@njit(cache=True)
def tracer_generator(rr,rho,bath):
    G=np.zeros((NH,NH)); C=np.zeros((len(WI),NH))
    for j in range(NR):
        a=TSRC[j]; b=TDST[j]
        if a<0: continue
        lam=rr[j]*PREF[j]*rho**DE[j]
        for k in range(REACT.shape[1]):
            i=REACT[j,k]
            if i>=0 and A[i]<10: lam*=max(0.,bath[i])
        if YE[j]:lam*=bath[NN]
        if a!=b: G[b,a]+=lam; G[a,a]-=lam
        for k in range(len(WI)):
            if WI[k]==j:C[k,a]=lam
    return G,C

def cold_map():
    G,_=tracer_generator(rate_values(1e6),0.,np.zeros(NN+1))
    unstable=np.flatnonzero(-np.diag(G)>1e-18); stable=np.array([k for k in range(NH) if k not in unstable])
    B=np.zeros((NH,NH)); B[stable,stable]=1
    # Integrated transient residence times give the exact absorbing branch probabilities.
    B[np.ix_(stable,unstable)]=G[np.ix_(stable,unstable)]@np.linalg.inv(-G[np.ix_(unstable,unstable)])
    return B,{NAMES[HEAVY[k]]:float(np.log(2)/(-G[k,k])) for k in unstable}

def trajectory_arrays(stage,index,stretch=1.):
    t=np.array(stage['time_s']); T=np.array(stage['temperature_K'][index]);rho=np.array(stage['density_g_cm3'][index])
    bath=np.zeros((len(t),NN+1));bath[:,:3]=np.array(stage['light_Y'][index]).T
    bath[:,NN]=stage['electron_Y'][index]
    rr=np.array([rate_values(v) for v in T]); x=t/t[-1]
    # Match the selected source's current carbon medians and reverse policy.
    for name,reverse in [('C12_p_to_N13_reaclib','N13_to_p_C12_reaclib'),
                         ('C13_p_to_N14_reaclib','N14_to_p_C13_reaclib')]:
        j=next(i for i,r in enumerate(D['rates']) if r['name']==name)
        k=next(i for i,r in enumerate(D['rates']) if r['name']==reverse)
        table=np.array(CARBON['tables'][name]['rows']); knots=np.log(table[:,0])
        old=np.array([rate_values(t*1e9)[j] for t in table[:,0]])
        poly=PchipInterpolator(knots,np.log(table[:,4]/old))
        factor=np.exp(poly(np.clip(np.log(np.maximum(T/1e9,1e-5)),knots[0],knots[-1])))
        rr[:,j]*=factor;rr[:,k]*=factor
    return x,rho,bath,rr,t[-1]*stretch

@njit(cache=True)
def tracer_rhs(t,x,grid,rhos,baths,rr,duration):
    G,C=tracer_generator(interpolate(t,grid,rr),interpolate(t,grid,rhos),interpolate(t,grid,baths))
    return np.concatenate((G@x[:NH],C@x[:NH]))*duration

@njit(cache=True)
def tracer_jac(t,x,grid,rhos,baths,rr,duration):
    G,C=tracer_generator(interpolate(t,grid,rr),interpolate(t,grid,rhos),interpolate(t,grid,baths))
    J=np.zeros((NH+len(WI),NH+len(WI)));J[:NH,:NH]=G;J[NH:,:NH]=C
    return J*duration

def tracer(stage,index,isotope,stretch=1.,rtol=2e-7):
    grid,rhos,baths,rr,duration=trajectory_arrays(stage,index,stretch)
    x0=np.zeros(NH+len(WI));x0[HI[NAMES.index(isotope)]]=1.
    args=(grid,rhos,baths,rr,duration)
    sol=solve_ivp(tracer_rhs,(0.,1.),x0,method='BDF',rtol=rtol,atol=1e-13,jac=tracer_jac,args=args)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:NH,-1], sol.y[NH:,-1], abs(sol.y[:NH,-1].sum()-1)

@njit(cache=True)
def bulk_rhs(t,x,grid,rhos,rr,duration):
    f=flows(x[:NN],interpolate(t,grid,rhos),interpolate(t,grid,rr))
    return np.concatenate((S@f,f[WI]))*duration

@njit(cache=True)
def bulk_jac(t,x,grid,rhos,rr,duration):
    rv=interpolate(t,grid,rr);rho=interpolate(t,grid,rhos);y=x[:NN]
    df=np.zeros((NR,NN));ye=np.dot(Z,y)
    for j in range(NR):
        pref=rv[j]*PREF[j]*rho**DE[j]
        for k in range(REACT.shape[1]):
            i=REACT[j,k]
            if i<0:continue
            val=pref
            for l in range(REACT.shape[1]):
                n=REACT[j,l]
                if n>=0 and l!=k:val*=y[n]
            df[j,i]+=val*(ye if YE[j] else 1.)
        if YE[j]:
            val=pref
            for k in range(REACT.shape[1]):
                i=REACT[j,k]
                if i>=0:val*=y[i]
            df[j]+=val*Z
    J=np.zeros((NN+len(WI),NN+len(WI)));J[:NN,:NN]=S@df;J[NN:,:NN]=df[WI]
    return J*duration

def bulk(stage,index,ratios,stretch=1.,rtol=2e-7):
    grid,rhos,baths,rr,duration=trajectory_arrays(stage,index,stretch)
    y0=np.zeros(NN)
    for name,v in ratios.items():y0[NAMES.index(name)]=v
    y0/=A@y0; norm=sum(y0[HEAVY]); x0=np.r_[y0,np.zeros(len(WI))]
    sol=solve_ivp(bulk_rhs,(0.,1.),x0,method='BDF',rtol=rtol,atol=1e-14,args=(grid,rhos,rr,duration),jac=bulk_jac)
    if not sol.success:raise RuntimeError(sol.message)
    y=sol.y[:NN,-1]
    return y/norm,sol.y[NN:,-1]/norm,abs(A@y-1),float((MASS@(y0-y))/norm)

def run_bulk(bins):
    with gzip.open(P/'trajectories.json.gz','rt') as f:stages=json.load(f)
    B,_=cold_map(); rows=[]
    for stage,s in stages.items():
        pure={n:v/sum(v for n,v in s['initial_Y'].items() if n!='p') for n,v in s['initial_Y'].items()}
        oxy={n:v*(.99 if n!='p' else 1.) for n,v in pure.items()};oxy['o16']=.01
        cases=[('control',pure,1.),('oxygen1pct',oxy,1.),
               ('mixed',dict(p=2.,c12=.25,c13=.25,n14=.25,n15=.24,o16=.01),1.)]
        if stage=='n15':cases +=[('large_time_surrogate',pure,100**(1/3)),('slow_time_surrogate',pure,2.)]
        for label,ratios,stretch in cases:
            y=np.zeros(NN);counts=np.zeros(len(WI));err=0.;energy=0.
            for i in sample(s,bins):
                a,c,e,q=bulk(s,i,ratios,stretch);y+=a/bins;counts+=c/bins;err=max(err,e);energy+=q/bins
            cold=B@y[HEAVY]
            row=dict(stage=stage,case=label,initial_number_ratios=ratios,time_stretch=stretch,
                hot_per_initial_heavy={NAMES[i]:float(v) for i,v in enumerate(y) if abs(v)>1e-13},
                cold_heavy_per_initial_heavy={NAMES[HEAVY[i]]:float(v) for i,v in enumerate(cold) if abs(v)>1e-13},
                watched_events_per_initial_heavy=dict(zip(WATCH,counts.tolist())),
                baryon_residual=err,hot_total_mass_Q_MeV_per_initial_heavy=energy,
                beyond_F_nuclei_per_initial_heavy=float(sum(cold[i] for i,n in enumerate(HEAVY) if Z[n]>=10)))
            rows.append(row);print(stage,label,'baryon',err,'Q',energy,flush=True)
    out=dict(scope='66-species nonlinear abundance integration on fixed saved T/rho; current carbon medians. Hydrodynamic, opacity, source energy and EOS feedback of mixtures and longer durations remain uncomputed.',bins=bins,cases=rows,
             code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (P/f'bulk-{bins}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

def sample(stage,bins):
    w=np.array(stage['mass_fraction']);c=np.cumsum(w)
    # Equal-mass stratified midpoint quadrature; repeat at a finer bin count.
    return [min(len(w)-1,int(np.searchsorted(c,(i+.5)/bins))) for i in range(bins)]

def run(bins):
    start=time.monotonic()
    with gzip.open(P/'trajectories.json.gz','rt') as f:stages=json.load(f)
    B,halflives=cold_map();rows=[]
    for stage,s in stages.items():
        for isotope in ['o16','o17','o18','f19']:
            hot=np.zeros(NH); counts=np.zeros(len(WI)); residual=0.
            for i in sample(s,bins):
                h,c,r=tracer(s,i,isotope);hot+=h/bins;counts+=c/bins;residual=max(residual,r)
            cold=B@hot
            row=dict(stage=stage,isotope=isotope,hot={NAMES[HEAVY[i]]:float(v) for i,v in enumerate(hot) if v>1e-15},
                cold={NAMES[HEAVY[i]]:float(v) for i,v in enumerate(cold) if v>1e-15},
                watched_events=dict(zip(WATCH,counts.tolist())), probability_residual=residual)
            row['quartet_return']=sum(cold[HI[NAMES.index(n)]] for n in ['c12','c13','n14','n15'])
            row['extended_pool']=sum(cold[i] for i,n in enumerate(HEAVY) if Z[n] in [6,7,8,9])
            row['beyond_F_fraction']=sum(cold[i] for i,n in enumerate(HEAVY) if Z[n]>=10)
            rows.append(row);print(stage,isotope,'return',row['quartet_return'],'beyondF',row['beyond_F_fraction'],flush=True)
    out=dict(scope=__doc__,bins=bins,half_lives_s=halflives,tracer=rows,
             code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             input_sha256={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['rates.json','trajectories.json.gz']},
             wall_s=time.monotonic()-start)
    (P/f'network-{bins}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bins',type=int,default=8);ap.add_argument('--bulk',action='store_true');a=ap.parse_args()
    run_bulk(a.bins) if a.bulk else run(a.bins)
