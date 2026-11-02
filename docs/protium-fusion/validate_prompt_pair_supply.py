"""Independent spectrum, source-ledger and finite-decay checks."""
from pathlib import Path
import json,gzip,argparse,math
import numpy as np
from scipy.integrate import quad
from prompt_pair_supply import Coefficients,fold_spectrum,MULT,THRESHOLD,KEYS,functions
from residual_lepton_births import half_life

def validate(data,root):
    checks=[]
    def check(name,ok,detail=None):
        checks.append({'name':name,'pass':bool(ok),'detail':detail})
        if not ok:raise AssertionError((name,detail))
    co=json.loads(gzip.decompress((data/'bricc-pair-coefficients.json.gz').read_bytes()))
    fine=json.loads(gzip.decompress((data/'bricc-pair-coefficients-refined.json.gz').read_bytes()))
    c=Coefficients(co);cf=Coefficients(fine)
    for z,values in c.data.items():
        check(f'positive monotone coefficient Z{z}',np.all(values>0) and np.all(np.diff(values,axis=0)>0))
    # Known direct BrIcc XML query, independently read before extraction.
    check('Z22 2MeV IPF multipoles',np.array_equal(c.values(22,2.),np.array([6.400e-4,2.505e-4,3.136e-4,1.347e-4,1.710e-4,7.738e-5])))
    # Spectrum with a discrete line plus a separately normalized flat continuum
    # crossing both thresholds and the upper coefficient boundary.
    lineE=2.;lineP=.2;lo=.5;hi=9.;density=.8/(hi-lo)
    sub={'LAW':1,'LEP':1,'INT':np.array([2]),'LIP':0,'E':np.array([1e6]),'ND':np.array([1]),
         'Ep':[np.array([lineE,lo,hi])*1e6], 'b':[np.array([[lineP],[density/1e6],[0.]])],
         'yields':{'INT':np.array([2]),'Eint':np.array([1e6]),'yi':np.array([1.])}}
    h=fold_spectrum(sub,22,c)
    check('mixed line continuum normalization',abs(h['photon_count'][0]-1)<2e-14)
    check('mixed line continuum energy',abs(h['photon_energy_MeV'][0]-(lineP*lineE+.8*(lo+hi)/2))<2e-14)
    for key,a,b in [('below_pair_count',lo,THRESHOLD),('near_threshold_count',THRESHOLD,1.1),('tabulated_count',1.1,8),('above_table_count',8,hi)]:
        target=density*(b-a)+(lineP if a<=lineE<b else 0)
        check('analytic region '+key,abs(h[key][0]-target)<2e-14)
    # Adaptive integration provides a second algorithm for the nonlinear
    # alpha/(1+alpha) integral. Split at known interpolation knots.
    for i,m in enumerate(MULT):
        g=lambda e: float(c.values(22,e)[i]/(1+c.values(22,e)[i]))
        q=sum(quad(g,a,b,epsabs=1e-15)[0] for a,b in zip(c.grid[:-1],c.grid[1:]))
        target=lineP*g(lineE)+density*q
        check('adaptive branch quadrature '+m,abs(h['pair_branch_'+m][0]-target)<2e-14)
    # Linear continuum density, with exact polynomial count/energy checks.
    sub['LEP']=2;sub['ND']=np.array([0]);sub['Ep']=[np.array([0.,10.])*1e6];sub['b']=[np.array([[0.],[.2/1e6]])]
    h=fold_spectrum(sub,23,c)
    check('linear spectrum norm',abs(h['photon_count'][0]-1)<2e-14)
    check('linear spectrum energy',abs(h['photon_energy_MeV'][0]-20/3)<2e-14)
    a=json.loads((data/'prompt-pair-supply.json').read_text());b=json.loads((data/'prompt-pair-supply-refined.json').read_text())
    original=json.loads((root/'network-ion-budget.json').read_text());refinement={}
    for d in [a,b]:
        tot=d['per_primary']
        check('reproduce previous photon energy order'+str(d['quadrature_order']),abs(tot['photon_energy_MeV']-original['prompt_photon_MeV_per_primary'])<1e-12)
        check('partition photon count',abs(sum(tot[k] for k in ['below_pair_count','near_threshold_count','tabulated_count','above_table_count'])-tot['photon_count'])<1e-13)
        check('partition photon energy',abs(sum(tot[k] for k in ['below_pair_energy_MeV','near_threshold_energy_MeV','tabulated_energy_MeV','above_table_energy_MeV'])-tot['photon_energy_MeV'])<1e-13)
        for m in MULT:
            check('pair threshold energy identity '+m,abs(tot['pair_transition_energy_MeV_'+m]-tot['pair_kinetic_MeV_'+m]-THRESHOLD*tot['pair_branch_'+m])<2e-16)
            check('ratio exceeds branch '+m,0<tot['pair_branch_'+m]<tot['pair_ratio_'+m]<tot['tabulated_count'])
        tail=d['discarded_MT5_photons_per_primary']
        check('MT5 photon partition',abs(sum(tail[k] for k in ['below_pair_count','near_threshold_count','tabulated_count','above_table_count'])-tail['photon_count'])<1e-18)
        check('MT5 energy partition',abs(sum(tail[k] for k in ['below_pair_energy_MeV','near_threshold_energy_MeV','tabulated_energy_MeV','above_table_energy_MeV'])-tail['photon_energy_MeV'])<1e-18)
        check('MT5 ordinary pair contribution is small',tail['pair_branch_E1']/tot['pair_branch_E1']<2e-6)
        for st in d['stages']:
            old=next(x for x in original['stages'] if x['ZA']==st['ZA'])
            check('stage energy '+str(st['ZA']),abs(st['per_stage_pn']['photon_energy_MeV']-old['photon_MeV_per_pn'])<1e-11)
    for m in MULT:
        key='pair_branch_'+m;err=abs(a['per_primary'][key]/b['per_primary'][key]-1);refinement[m]=err
        check('coefficient refinement '+m,err<1e-5,err)
    d=json.loads((data/'residual-lepton-births.json').read_text())
    check('direct residual total matches extraction',abs(d['direct_residuals_per_primary']-original['extraction_events_per_primary'])<1e-12)
    check('all encountered weak half-lives parsed',len(d['missing_half_lives'])==0)
    for raw,expect in [('2.5 ms',.0025),('1.5 d',129600.),('81 ky',81000*365.25*86400),('271 Py',271e15*365.25*86400)]:
        check('half-life units '+raw,math.isclose(half_life({'half_life_raw':raw})['seconds'],expect,rel_tol=1e-14))
    for population in ['ground_state','fastest_known_weak_state']:
        ss=[s for s in d['scenarios'] if s['population']==population]
        check('residence monotonic '+population,np.all(np.diff([s['total_neutral_rate_proxy_per_primary'] for s in ss])>0))
        range_ok=[];linear_ok=[]
        for s in ss:
            for item in s['items']:
                row=next(r for r in d['weak_candidate_rows'] if r['ZA']==item['ZA'])
                x=math.log(2)*s['residence_s']/item['half_life_s']
                val=item['neutral_rate_decay_proxy_per_primary']
                range_ok.append(0<=val<=row['births_per_primary'])
                # Small-time rate provides an independent first-order bound.
                linear_ok.append(val<=row['births_per_primary']*x*(1+1e-14))
        check('decay probability range '+population,all(range_ok))
        check('linear survival upper bound '+population,all(linear_ok))
    return {'checks':checks,'passed':len(checks),'coefficient_grid_relative_changes':refinement,
            'scope':'Numerical integration and data accounting; not qualification of missing nuclear channels or collectors.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path(__file__).parent);ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--out',type=Path);args=ap.parse_args()
    d=validate(args.data,args.root);out=args.out or args.data/'prompt-pair-validation.json';out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(d['passed'],'checks passed',d['coefficient_grid_relative_changes'])
