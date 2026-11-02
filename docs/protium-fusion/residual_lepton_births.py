"""Finite-residence weak-decay sensitivity for directly produced heavy residuals.

Uses NUBASE neutral-atom half-lives, treats every EC/B+ entry as if it emitted
one positron, and does not count cold capture/decay chains as in-plasma births.
An all-fastest-known-weak-isomer assignment is an intentionally severe population
scenario. Neither is a bound on charge-state-dependent weak rates.
"""
from pathlib import Path
import argparse,json,math,re,hashlib
from collections import defaultdict
from catalyst_catalog import nubase_records,NUBASE
from catalyst_network import load_data,stage,HOT,label
ROOT=Path(__file__).parent
UNITS={'ys':1e-24,'zs':1e-21,'as':1e-18,'fs':1e-15,'ps':1e-12,'ns':1e-9,'us':1e-6,'ms':1e-3,'s':1.,'m':60.,'h':3600.,'d':86400.,'y':365.25*86400.,'ky':1e3*365.25*86400.,'My':1e6*365.25*86400.,'Gy':1e9*365.25*86400.,'Ty':1e12*365.25*86400.,'Py':1e15*365.25*86400.}

def half_life(r):
    raw=r['half_life_raw']
    if raw=='stbl':return None
    m=re.fullmatch(r'([<>~]?)([0-9]+(?:\.[0-9]*)?(?:[eE][+\-]?[0-9]+)?)[#]?\s*([A-Za-z]+)',raw)
    if not m or m[3] not in UNITS:return None
    return {'seconds':float(m[2])*UNITS[m[3]],'qualifier':m[1],'raw':raw}

def run(root):
    base=json.loads((root/'network-ion-budget.json').read_text());data=load_data(root/'beam-data/tendl2025-proton-recycling.zip');nb=nubase_records();rates=defaultdict(float)
    for za in HOT:
        w=base['network']['processor_pn_rates'][label(za)]
        for e in stage(za,data[za],base['T_MeV'])['events']:
            if e['mt']!=103:rates[tuple(e['residual'])]+=w*e['rate']
    rows=[];missing=[]
    for za,rate in sorted(rates.items(),key=lambda x:-x[1]):
        if rate==0:continue
        states=[]
        for key,r in nb.items():
            if key[:2]!=za or key[2]>=8:continue
            mode=r['decay_modes'];has_plus=('B+' in mode or 'EC=' in mode or 'e+=' in mode);has_minus='B-' in mode
            if not (has_plus or has_minus) or r['half_life_raw']=='stbl':continue
            half=half_life(r)
            if half is None:missing.append({'ZA':list(za),'state':key[2],'rate_per_primary':rate,'raw':r['half_life_raw'],'modes':mode});continue
            states.append({'state':key[2],'half_life':half,'decay_modes':mode,
                'plus_or_EC':has_plus,'minus':has_minus})
        if states:rows.append({'ZA':list(za),'births_per_primary':rate,'states':states})
    scenarios=[]
    for residence in [1e-6,1e-3,.1,1.,10.,100.]:
        for population in ['ground_state','fastest_known_weak_state']:
            items=[]
            for row in rows:
                cand=[s for s in row['states'] if population!='ground_state' or s['state']==0]
                if not cand:continue
                s=min(cand,key=lambda x:x['half_life']['seconds'])
                decay=row['births_per_primary']*(-math.expm1(-math.log(2)*residence/s['half_life']['seconds']))
                items.append({'ZA':row['ZA'],'state':s['state'],'half_life_s':s['half_life']['seconds'],
                              'neutral_rate_decay_proxy_per_primary':decay,
                              'positron_or_EC_proxy':decay if s['plus_or_EC'] else 0.,
                              'beta_minus_proxy':decay if s['minus'] else 0.})
            total=sum(i['neutral_rate_decay_proxy_per_primary'] for i in items)
            scenarios.append({'residence_s':residence,'population':population,'total_neutral_rate_proxy_per_primary':total,
                'positron_or_EC_proxy_per_primary':sum(i['positron_or_EC_proxy'] for i in items),
                'beta_minus_proxy_per_primary':sum(i['beta_minus_proxy'] for i in items),
                'items':sorted(items,key=lambda i:-i['neutral_rate_decay_proxy_per_primary'])})
    return {'scope':__doc__,'inputs':{'NUBASE_sha256':hashlib.sha256(NUBASE.read_bytes()).hexdigest()},
        'direct_residuals_per_primary':sum(rates.values()),'weak_candidate_rows':rows,'missing_half_lives':missing,
        'scenarios':scenarios,'limitations':['EC entries do not create positrons in general. The plus/EC column is an intentionally inclusive neutral-rate proxy.',
         'Fully stripped ions suppress bound-electron capture; continuum capture and bound-state beta decay require actual charge-state and lepton-density models.',
         'Every weak decay in the fastest-state scenario uses the total half-life, including competing IT. Its unknown production fraction is set to one.',
         'Hot retained catalyst decays, delayed chains following a first in-plasma decay, and missing states are outside this direct-residual calculation.',
         'The residence is assigned; selective extraction on that timescale is not established.']}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=ROOT/'residual-lepton-births.json');args=ap.parse_args()
    d=run(args.root);args.out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print('missing',d['missing_half_lives'])
    for s in d['scenarios']:print(s['residence_s'],s['population'],s['total_neutral_rate_proxy_per_primary'],s['items'][:2])
