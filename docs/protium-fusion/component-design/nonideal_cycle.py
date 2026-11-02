"""Saved-ledger comparison of all16 combinations from the eight Coulomb cases.
Perfect eventual CNO/H return, same nuclear ledger and heat/driver efficiencies.
No further hydro, receiver or plant-mass calculation.
"""
from pathlib import Path
import itertools,json,hashlib
from material_cycle import load_stage,solve_cycle
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 comparison=json.loads((HERE/'nonideal_comparison.json').read_text());pairs={}
 for r in comparison['cases']:
  pairs[(r['stage'],r['H_per_catalyst'])]={'control':load_stage(r['control_path']),'stress':load_stage(r['archived_record'])}
 rows=[]
 for hs in itertools.product([1,3],[1,2],[2,2.5],[1,4]):
  cases={}
  for mode in ['control','stress']:
   stages={s:pairs[(s,h)][mode] for s,h in zip(['A','B1','B2','C13'],hs)}
   try:cases[mode]=solve_cycle(stages,carbon_return=1.,hydrogen_return=1.,eta_driver=.6,eta_heat=.6,capture=.99)
   except Exception as exc:cases[mode]={'feasible':False,'error':repr(exc)}
  rows.append({'mixture_H_per_catalyst':dict(zip(['A','B1','B2','C13'],hs)),'cases':cases})
 ranks={mode:sorted([{'mixture':r['mixture_H_per_catalyst'],'preaux_net_MeV_per_He':r['cases'][mode]['preaux_net_MeV_per_He']} for r in rows if r['cases'][mode]['feasible']],key=lambda r:r['preaux_net_MeV_per_He'],reverse=True) for mode in ['control','stress']}
 out={'scope':__doc__,'code_sha256':sha(__file__),'ledger_sha256':sha(HERE/'material_cycle.py'),'comparison_sha256':sha(HERE/'nonideal_comparison.json'),'cases':rows,'ranked_preaux':ranks,'interpretation':'Pre-auxiliary electric energy per helium is only one design metric. Source/cadence/receiver/stock/structure mass and processing are unchanged inputs rather than reoptimized here. Perfect eventual return does not imply instantaneous separation or absence of gas inventory.'}
 (HERE/'nonideal_cycle.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'ranks':ranks,'feasible_counts':{k:len(v) for k,v in ranks.items()}},indent=2))
if __name__=='__main__':main()
