"""Compare supplied completed pulse recipes at separate H/CNO return fractions.

Stage-map JSON maps four recipe keys to a path or {record,initial_Y,
initial_counts,trace_path}. Paths are relative to --root unless absolute.
This tool never infers input stoichiometry from a label or filename.
"""
import argparse,json,hashlib,math
from pathlib import Path
import material_cycle as mc

def compare(stage_map,root,*,targets=(.0056,.0063,.007),hydrogen_returns=(1.,.9999,.999,.99,.9),
            carbon_returns=(.99,.999,.9995,.9999,1.),eta_driver=.6,eta_heat=.6,capture=.99):
 root=Path(root).resolve();stages={}
 def resolve(p):
  p=Path(p);return p if p.is_absolute() else root/p
 for key,spec in stage_map.items():
  spec={'record':spec} if isinstance(spec,str) else dict(spec)
  path=resolve(spec.pop('record'));spec['trace_path']=resolve(spec['trace_path']) if spec.get('trace_path') else None
  stages[key]=mc.load_stage(path,**spec)
 kwargs={'eta_driver':eta_driver,'eta_heat':eta_heat,'capture':capture}
 half={'n13':597.9978026337243,'o14':70.59978360562287,'o15':121.99951623253882}
 def decorate(case):
  if not case.get('feasible'):return case
  case['processing']={'cold_output_atoms_per_He':{},'gross_input_atoms_per_He':{},'hot_radioactive_birth_atoms_per_He':{}}
  for key,s in stages.items():
   a=case['stage_attempts_per_He'][key]
   mc.add(case['processing']['cold_output_atoms_per_He'],s['cold_output'],a)
   mc.add(case['processing']['gross_input_atoms_per_He'],s['input_counts'],a)
   mc.add(case['processing']['hot_radioactive_birth_atoms_per_He'],{n:s['hot_output'].get(n,0.) for n in half},a)
  case['processing']['ideal_steady_radioactive_inventory_atoms_per_He_s']= {n:v*half[n]/math.log(2) for n,v in case['processing']['hot_radioactive_birth_atoms_per_He'].items()}
  case['processing']['decay_inventory_scope']='Multiply by He/s for all-flight-plus-catcher parent inventory at steady production; no additional reactions or selective loss. Post-wall allocation needs each new arrival history.'
  return case
 targetcases=[]
 for target in targets:
  for rh in [None,*hydrogen_returns]:
   try:c=decorate(mc.solve_supply_target(stages,target=target,total_budget=.007,hydrogen_return=rh,**kwargs))
   except (ValueError,ArithmeticError,mc.np.linalg.LinAlgError) as e:c={'feasible':False,'reason':str(e)}
   targetcases.append({'fuel_C_H_target':target,'structural_C_H_reserve':.007-target,'hydrogen_return_mode':'coupled_to_CNO' if rh is None else 'independent','hydrogen_return_requested':rh,'result':c})
 scans=[]
 for rc in carbon_returns:
  for rh in hydrogen_returns:
   try:c=decorate(mc.solve_cycle(stages,carbon_return=rc,hydrogen_return=rh,**kwargs))
   except (ValueError,ArithmeticError,mc.np.linalg.LinAlgError) as e:c={'feasible':False,'reason':str(e)}
   scans.append({'CNO_return':rc,'H_return':rh,'result':c})
 return {'scope':__doc__,'normalization':'One NET helium4 nucleus, actual atomic input/output masses; physical preparation/capture efficiencies remain explicit requirements.',
  'loaded_stages':stages,'perfect_eventual_return_case':decorate(mc.solve_cycle(stages,carbon_return=1.,hydrogen_return=1.,**kwargs)),'target_cases':targetcases,'independent_return_scan':scans,'network_half_lives_s':half,
  'positive_He_and_H_required':True,'module_sha256':mc.hashfile(mc.__file__),'code_sha256':mc.hashfile(__file__)}

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--stage-map',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--eta-driver',type=float,default=.6);ap.add_argument('--eta-heat',type=float,default=.6);ap.add_argument('--capture',type=float,default=.99)
 a=ap.parse_args();out=compare(json.loads(a.stage_map.read_text()),a.root,eta_driver=a.eta_driver,eta_heat=a.eta_heat,capture=a.capture);out['stage_map_sha256']=mc.hashfile(a.stage_map);a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'output':str(a.output),'target_cases':len(out['target_cases']),'return_scan_cases':len(out['independent_return_scan'])}))
