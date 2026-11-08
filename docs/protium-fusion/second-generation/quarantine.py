#!/usr/bin/env python3
"""Test repeated Ne/Na/Mg/Al return in actual C13 host baths, with Si28 terminal.
This exposes competing onward captures instead of inferring recycle as 1/p from
one isolated favorable reaction. The A<=28 truncation is an explicit boundary.
"""
import gzip
import json
from pathlib import Path
import numpy as np
import network as n
P=Path(__file__).resolve().parent

def main():
    with gzip.open(P/'trajectories.json.gz','rt') as f:st=json.load(f)['c13']
    B,_=n.cold_map()
    states=['ne20','ne21','ne22','na23','mg24','mg25','mg26','al26','al27']
    rows=[]
    for isotope in states:
        hot=np.zeros(n.NH);counts=np.zeros(n.NT)
        for i in n.sample(st,8):
            h,c,e=n.tracer(st,i,isotope);hot+=h/8;counts+=c/8
        cold=B@hot
        rows.append(dict(isotope=isotope,cold={n.NAMES[n.HEAVY[j]]:float(v) for j,v in enumerate(cold) if v>1e-16},
            hot_ledger=dict(zip(n.TRACK_NAMES,counts[len(n.WI):].tolist()))))
        print(isotope,'one pass return',sum(cold[n.HI[n.NAMES.index(v)]] for v in ['c12','c13','n14','n15','o16','o17','o18','f19']),flush=True)
    M=np.array([[r['cold'].get(v,0.) for r in rows] for v in states])
    ret=np.array([sum(r['cold'].get(v,0.) for v in ['c12','c13','n14','n15','o16','o17','o18','f19']) for r in rows])
    results=[]
    for isotope in ['ne20','ne21','ne22']:
        initial=np.array([float(v==isotope) for v in states]);visits=np.linalg.solve(np.eye(len(states))-M,initial)
        results.append(dict(isotope=isotope,visits=float(sum(visits)),return_to_extended_pool_probability=float(ret@visits),
            terminal_other_probability=float(1-ret@visits),expected_C13_allocated_preparation_MeV=float(sum(visits))*.4575048140317259,
            conservation_matrix_residual=float(max(abs((np.eye(len(states))-M)@visits-initial)))))
    (P/'quarantine-results.json').write_text(json.dumps(dict(scope=__doc__,bins=8,stage='c13',states=states,
        one_pass=rows,results=results,condition_number=float(np.linalg.cond(np.eye(len(states))-M)),
        qualification='Si28 and products outside the declared transient states are quarantined. The unresolved isotope boundary, rate updates, neutron bath and finite sorting/processing conditions limit any irreversible-loss inference.'),indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
