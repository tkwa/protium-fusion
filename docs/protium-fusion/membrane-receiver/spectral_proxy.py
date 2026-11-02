#!/usr/bin/env python3
"""Planck-shape proxy for saved grey escape; not a spectral radiation calculation.
Run with python3 spectral_proxy.py. Numpy only. No network or hydro operation.
NIST mass coefficients transcribed 2026-10-30 via web tool, 1--1000 keV.
Below 1 keV and above 1 MeV, report zero/full-retention bounds rather than extrapolate.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
KB_KEV=8.617333262e-8
# Energy keV; mu/rho and mu_en/rho in cm^2/g. Si has its actual K-edge jump.
E=[1,1.5,2,3,4,5,6,8,10,15,20,30,40,50,60,80,100,150,200,300,400,500,600,800,1000]
DATA={
'C':[[2211,2209],[700.2,699],[302.6,301.6],[90.33,89.63],[37.78,37.23],[19.12,18.66],[10.95,10.54],[4.576,4.242],[2.373,2.078],[.8071,.5627],[.442,.2238],[.2562,.06614],[.2076,.03343],[.1871,.02397],[.1753,.02098],[.161,.02037],[.1514,.02147],[.1347,.02449],[.1229,.02655],[.1066,.02870],[.09546,.02950],[.08715,.02969],[.08058,.02956],[.07076,.02885],[.06361,.02792]],
'O':[[4590,4576],[1549,1545],[694.9,692.6],[217.1,215.8],[93.15,92.21],[47.9,47.15],[27.7,27.08],[11.63,11.16],[5.952,5.565],[1.836,1.545],[.8651,.6179],[.3779,.1729],[.2585,.07530],[.2132,.04414],[.1907,.03207],[.1678,.02468],[.1551,.02355],[.1361,.02506],[.1237,.02679],[.107,.02877],[.09566,.02953],[.08729,.02971],[.0807,.02957],[.07087,.02887],[.06372,.02794]],
'B':[[1229,1228],[376.6,375.9],[159.7,159.1],[46.67,46.17],[19.27,18.86],[9.683,9.332],[5.538,5.223],[2.346,2.072],[1.255,1.006],[.4827,.2698],[.3014,.1084],[.2063,.03506],[.1793,.02084],[.1665,.01737],[.1583,.01680],[.1472,.01785],[.1391,.01940],[.1243,.02255],[.1136,.02453],[.09862,.02654],[.08834,.02731],[.08065,.02749],[.07460,.02737],[.06549,.02671],[.05890,.02586]],
'N':[[3311,3306],[1083,1080],[476.9,475.5],[145.6,144.7],[61.66,60.94],[31.44,30.86],[18.09,17.59],[7.562,7.170],[3.879,3.545],[1.236,.9715],[.6178,.3867],[.3066,.1099],[.2288,.05051],[.198,.03217],[.1817,.02548],[.1639,.02211],[.1529,.02231],[.1353,.02472],[.1233,.02665],[.1068,.02873],[.09557,.02952],[.08719,.02969],[.08063,.02956],[.07081,.02886],[.06364,.02792]],
'Si':[[1570,1567],[535.5,533.1],[2777,2669],[978.4,951.6],[452.9,442.7],[245,240],[147,143.9],[64.68,63.13],[33.89,32.89],[10.34,9.794],[4.464,4.076],[1.436,1.164],[.7012,.4782],[.4385,.243],[.3207,.1434],[.2228,.06896],[.1835,.04513],[.1448,.03086],[.1275,.02905],[.1082,.02932],[.09614,.02968],[.08748,.02971],[.08077,.02951],[.07082,.02875],[.06361,.02778]]}
AM={'C':12.011,'O':15.999,'Si':28.0855,'B':10.81,'N':14.007}
MATS={'silica':{'rho':2200,'cp':1000,'parts':{'Si':1,'O':2}},'carbon':{'rho':2000,'cp':1000,'parts':{'C':1}},'BN':{'rho':2100,'cp':1000,'parts':{'B':1,'N':1}},'SiC':{'rho':3210,'cp':1000,'parts':{'Si':1,'C':1}}}
SOURCES={x:f'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z{z:02}.html' for x,z in [('C',6),('O',8),('Si',14),('B',5),('N',7)]}
SOURCES['definition']='https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def coeff(x,grid):
    ee=np.array(E,float);dd=np.array(DATA[x],float)
    if x=='Si':
        # Explicit edge on both sides avoids interpolating through the physical jump.
        ee=np.insert(ee,2,[1.83890*(1-1e-12),1.83890]);dd=np.insert(dd,2,[[309.2,307.0],[3192,3059]],axis=0)
    return np.stack([np.exp(np.interp(np.log(grid),np.log(ee),np.log(dd[:,j]))) for j in [0,1]])*.1

def spectrum(grid,kt):
    x=grid[None,:]/kt[:,None]
    # High x is negligible; clipping prevents overflow without affecting support.
    return grid[None,:]**3/np.expm1(np.minimum(x,700))/(math.pi**4/15*kt[:,None]**4)

def planck_parts(kt,df,grid):
    # Weighted incident energy density, normalized by full analytic Planck integral.
    sp=np.sum(spectrum(grid,kt)*df[:,None],axis=0)
    lowg=np.geomspace(1e-7,1,1800);hig=np.geomspace(1000,100000,1000)
    low=float(np.trapz(np.sum(spectrum(lowg,kt)*df[:,None],axis=0),lowg))
    high=float(np.trapz(np.sum(spectrum(hig,kt)*df[:,None],axis=0),hig))
    return sp,low,high

def calculate(points=6000):
    pulse=HERE/'pulse-loading.json';pl=json.loads(pulse.read_text());grid=np.unique(np.r_[np.geomspace(1,1000,points),1.8389*(1-1e-10),1.8389, E])
    mc={}
    for m,p in MATS.items():
        mass=sum(AM[x]*n for x,n in p['parts'].items());mc[m]=sum(AM[x]*n/mass*coeff(x,grid) for x,n in p['parts'].items())
    stage_rows={};checks=[]
    for name,s in pl['stages'].items():
        pr=s['prompt'];f=np.array(pr['cumulative_fluence_J_m2']);df=np.diff(f);temp=np.array(pr['outer_temperature_K']);dtime=np.diff(np.array(pr['emission_time_s']))
        assert np.all(df>=-1e-9);df=np.maximum(df,0);kt=.5*(temp[:-1]+temp[1:])*KB_KEV
        assert np.min(kt)>0
        sets={}
        for label,ktcase in [('midpoint',kt),('bin_start',temp[:-1]*KB_KEV),('bin_end',temp[1:]*KB_KEV),('half_temperature_sensitivity',kt*.5),('double_temperature_sensitivity',kt*2)]:
            sp,low,high=planck_parts(ktcase,df,grid);resolved=float(np.trapz(sp,grid));norm=(resolved+low+high)/sum(df)-1
            checks.append(abs(norm)<5e-6)
            rows=[]
            for m,p in MATS.items():
                mu,men=mc[m]
                for nm in [10,30,100,300,1000,3000]:
                    am=p['rho']*nm*1e-9;lost=-np.expm1(-mu*am);dep=men/mu*lost
                    heat=float(np.trapz(sp*dep,grid));inter=float(np.trapz(sp*lost,grid));delta=heat/(am*p['cp'])
                    # mu_en deposition, and a looser all-interaction transfer upper screen.
                    rows.append({'material':m,'thickness_nm':nm,'areal_mass_kg_m2':am,'energy_transfer_proxy_resolved_J_m2':heat,'all_interaction_energy_upper_resolved_J_m2':inter,'sensible_delta_K_local_retention1_resolved':delta,'resolved_surface_sensible_delta_K_upper_no_conduction':float(np.trapz(sp*men,grid)/p['cp']),'sensible_delta_K_unresolved_tail_full_absorption_addition':(low+high)/(am*p['cp']),'transfer_fraction_resolved':heat/sum(df),'uncollided_interaction_fraction_resolved':inter/sum(df)})
                    checks.extend([0<=heat<=inter+1e-12,inter<=resolved+1e-8])
            sets[label]={'escape_weighted_kT_keV':float(np.dot(df,ktcase)/sum(df)),'fluence_below_1keV_J_m2':low,'fluence_above_1MeV_J_m2':high,'fraction_below_1keV':low/sum(df),'fraction_above_1MeV':high/sum(df),'normalization_relative_error':norm,'rows':rows}
        stage_rows[name]={'prompt_fluence_J_m2':float(sum(df)),'duration_s':float(sum(dtime)),'peak_kT_keV':float(max(temp)*KB_KEV),'quadrature':sets}
    return {'scope':'Planck spectral-shape proxy assigned to each saved grey escape increment at outer-cell temperature; not a predicted escaping spectrum. No new hydro, radiation transport, or lifetime inference.','temperature_scope':'A/B2/C13 weighted kT13--16keV; B1 is softer. Bin-endpoint variants are sampling checks; temperature factors are model sensitivities, not confidence intervals.','energy_scope':'mu_en/rho bulk energy-transfer proxy with mixture mass additivity. Local retention factor is an independent unknown; this does not simulate secondary electrons, incident plasma electrons, sheath return, or photon multiple scattering. All-interaction energy is a looser bound for each single crossing.','input_sha256':{'pulse-loading.json':sha(pulse)},'code_sha256':sha(__file__),'table_source_urls':SOURCES,'tabulated_data':{'energy_keV':E,'coefficient_units':'cm2/g','rows_mu_mu_en':DATA,'silicon_extra_edge_rows':[[1.83890,309.2,307],[1.83890,3192,3059]]},'radius_m':pl['radius_m'],'material_assumptions':MATS,'energy_grid_points':len(grid),'stages':stage_rows,'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':all(checks)}}
if __name__=='__main__':
    out=calculate();(HERE/'spectral-proxy.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'summary':{s:{'weighted_kT':d['quadrature']['midpoint']['escape_weighted_kT_keV'],'soft_fraction':d['quadrature']['midpoint']['fraction_below_1keV'],'silica10nm':next(x for x in d['quadrature']['midpoint']['rows'] if x['material']=='silica' and x['thickness_nm']==10),'silica1um':next(x for x in d['quadrature']['midpoint']['rows'] if x['material']=='silica' and x['thickness_nm']==1000)} for s,d in out['stages'].items()}},indent=2))
