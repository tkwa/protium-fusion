"""Integrate one explicit four-stage variant through material/receiver/plant ledgers.

No run discovery or implicit mixing of source, mesh, EOS, or nuclear-rate revisions.
Metadata adapters preserve the exact physical records and initial profiles while
supplying declared physical duration and trace paths to the existing reducer.
Archived results embed the reduced inputs for replay without raw hydro traces.
"""
from pathlib import Path
import argparse,copy,hashlib,importlib.util,json,math,sys

CAT={'n14':'A','n15':'B1','c12':'B2','c13':'C13'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def resolve(root,p):
 p=Path(p);return p.resolve() if p.is_absolute() else (root/p).resolve()
def close(a,b):return math.isclose(float(a),float(b),rel_tol=1e-8,abs_tol=1e-14)
def first(*values):return next((v for v in values if v is not None),None)

def extract_explicit(stage_map_path,research,workdir,receiver,material,*,provenance_path=None,thermal_points=65536):
 definition=read(stage_map_path);mapping=definition.get('stage_map',definition)
 if set(mapping)!=set(CAT):raise ValueError('Stage map must name exactly n14,n15,c12,c13')
 provenance=read(provenance_path) if provenance_path else {};detail=provenance.get('stages',{});selected=[];archive={};checks=[]
 sources={str(stage_map_path):sha(stage_map_path)}
 if provenance_path:sources[str(provenance_path)]=sha(provenance_path)
 for cat,tag in CAT.items():
  spec={'record':mapping[cat]} if isinstance(mapping[cat],str) else dict(mapping[cat]);extra=detail.get(cat,{})
  rec=resolve(research,spec['record']);d=read(rec)
  if d.get('stage')!=tag:raise ValueError(f'{cat}: stage label disagrees with map')
  initial=resolve(research,first(spec.get('initial'),extra.get('initial_profile'),str(rec.with_suffix('.initial.json'))))
  inputs=resolve(research,first(spec.get('inputs'),extra.get('input_record'),str(rec.with_suffix('.inputs.json'))))
  inp=read(inputs);profile=read(initial)
  trace=resolve(research,first(spec.get('trace'),extra.get('raw_trace_path'),d.get('trace_path'),str(rec.with_suffix('.npz'))))
  eos=spec.get('eos','FastEOS')
  if eos!='FastEOS':raise ValueError(f'{cat}: {eos} needs its own EOS reconstruction adapter')
  # Refuse explicit nonideal metadata rather than silently apply the ideal table.
  for owner in [d,inp]:
   for key in ['coulomb_scale','coulomb_multiplier','coulomb_strength','nonideal_eos']:
    if owner.get(key):raise ValueError(f'{cat}: nonideal marker {key} requires an EOS adapter')
   if owner.get('eos_model') not in [None,'FastEOS','ideal','ideal_FastEOS']:raise ValueError(f'{cat}: unknown eos_model')
  hashes={str(x):sha(x) for x in [rec,initial,inputs,trace]}
  for path,key in [(rec,'result_sha256'),(initial,'initial_profile_sha256'),(inputs,'input_record_sha256'),(trace,'raw_trace_sha256')]:
   for expected in [spec.get(key),extra.get(key)]:
    if expected and sha(path)!=expected:raise ValueError(f'{cat}: declared hash mismatch for {path.name}')
  if extra.get('result') and resolve(research,extra['result'])!=rec:raise ValueError(f'{cat}: provenance refers to another result')
  duration=first(spec.get('duration_s'),inp.get('duration_s'),inp.get('physical_duration_s'),extra.get('physical_duration_s'),d.get('mixture_inputs',{}).get('duration_s'))
  if duration is None or duration<=0:raise ValueError(f'{cat}: explicit requested physical duration required')
  if not d.get('success') or d.get('termination_reason')!='final_time' or not close(d.get('elapsed_s',-1),duration):raise ValueError(f'{cat}: requested physical endpoint was not completed')
  for candidate in [inp.get('duration_s'),inp.get('physical_duration_s'),extra.get('physical_duration_s'),d.get('mixture_inputs',{}).get('duration_s')]:
   if candidate is not None and not close(candidate,duration):raise ValueError(f'{cat}: contradictory physical durations')
  adapter=copy.deepcopy(d);adapter['trace_path']=str(trace)
  if not isinstance(adapter.get('mixture_inputs'),dict):adapter['mixture_inputs']={}
  adapter['mixture_inputs']['duration_s']=duration
  input_adapter=copy.deepcopy(inp);input_adapter['duration_s']=duration
  # Add no invented composition. The loader can extract the actual trace mean.
  path=workdir/(cat+'__'+rec.stem+'.json');write(path,adapter);write(path.with_suffix('.initial.json'),profile);write(path.with_suffix('.inputs.json'),input_adapter)
  s=material.load_stage(path,trace_path=trace)
  for declared in [inp.get('initial_Y'),extra.get('initial_Y'),spec.get('initial_Y')]:
   if declared is not None and max(abs(declared.get(n,0)-s['initial_Y'].get(n,0)) for n in set(declared)|set(s['initial_Y']))>1e-12:raise ValueError(f'{cat}: input composition disagrees with actual trace')
  r=receiver.reduce_case(path,research.parents[1],radius=1e9,thermal_points=thermal_points)
  if not r.get('completed_requested_duration') or not r.get('usable_ballistic'):raise ValueError(f'{cat}: receiver projection is not qualified for completed ballistic accounting')
  if r['radius_m']!=1e9:raise ValueError('Plant adapter requires receiver caches at explicit 1e9 m reference')
  # Keep physical fields bit-identical; adapter changes only path/duration metadata.
  physicalkeys=[k for k in d if k not in ['trace_path','mixture_inputs']]
  if any(adapter[k]!=d[k] for k in physicalkeys):raise ValueError('Adapter changed a physical field')
  checks.extend([{'name':cat+' completed physical duration','passed':True},{'name':cat+' source hashes and physical adapter','passed':True},{'name':cat+' receiver energy reconstruction','passed':True}])
  entry={'name':rec.stem,'path':str(rec),'record':adapter,'receiver':r,'stage_data':s}
  selected.append(entry);archive[cat]={'record_path':str(rec),'record_sha256':sha(rec),'initial_profile':profile,'input_record':inp,'input_sha256':hashes,'eos':'FastEOS','metadata_adapter':{'trace_path':str(trace),'requested_duration_s':duration,'scope':'Only file-location and physical-duration fields normalized; original source/profile energy and composition preserved.'},'selected':entry}
  write(workdir/(cat+'.receiver.json'),r)
  sources.update(hashes)
 return selected,archive,checks,sources

def result_checks(selected,cycle,plant):
 checks=[]
 for k in ['baryon_residual_per_He','mass_defect_Q_residual_MeV_per_He','energy_partition_residual_MeV_per_He']:
  value=cycle[k];checks.append({'name':k,'value':value,'passed':abs(value)<1e-7})
 if 'net_W' in plant:
  h=plant['nuclear_non_neutrino_W'];res=(h-plant['escaping_source_W']-plant['net_W']-plant['cold_rejection_W'])/h
  checks.append({'name':'independent global energy closure','value':res,'passed':abs(res)<1e-10})
  checks.append({'name':'independent mass closure','value':sum(plant['mass_breakdown_kg'].values())/plant['mass_kg']-1,'passed':math.isclose(sum(plant['mass_breakdown_kg'].values()),plant['mass_kg'],rel_tol=1e-10)})
  mean=sum(cycle['stage_pulse_data'][k]['pulses_per_He']*next(s['receiver']['timing_at_radius']['mass999']['occupied_s'] for s in selected if s['record']['stage']==CAT[k]) for k in CAT)*plant['helium_per_s']*plant['radius_m']/1e9
  checks.append({'name':'independent main-window occupancy','value':mean-plant['sequential_occupancy_mass999'],'passed':abs(mean-plant['sequential_occupancy_mass999'])<1e-12})
 return checks

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--research-root',type=Path,required=True);src=ap.add_mutually_exclusive_group(required=True);src.add_argument('--stage-map',type=Path);src.add_argument('--replay',type=Path);ap.add_argument('--provenance',type=Path);ap.add_argument('--family',default='explicit-variant');ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reference',type=Path);ap.add_argument('--radius',type=float,default=1.1e9);ap.add_argument('--thermal-points',type=int,default=65536);ap.add_argument('--carbon-return',type=float,default=1.);ap.add_argument('--hydrogen-return',type=float,default=1.);a=ap.parse_args();research=a.research_root.resolve();a.work_dir.mkdir(parents=True,exist_ok=True)
 sys.path.insert(0,str(research/'component-design'));import material_cycle as mc
 receiver=module('variant_receiver',research/'mixture-design/receiver.py');plantmod=module('variant_plant',research/'mixture-design/plant_compare.py')
 sourcecodes=[Path(__file__),research/'component-design/material_cycle.py',research/'mixture-design/receiver.py',research/'mixture-design/plant_compare.py',research/'cold-sink-redesign/station.py',research/'open-loop-redesign/operating-cases.json',research/'cold-sink-redesign/cases.json'];codehash={str(p):sha(p) for p in sourcecodes}
 receiver.self_check()
 if a.replay:
  old=read(a.replay);archive=old['stage_inputs'];selected=[archive[k]['selected'] for k in CAT];checks=[{'name':'explicit reduced-input replay','passed':True}];sources={str(a.replay):sha(a.replay)}
 else:selected,archive,checks,sources=extract_explicit(a.stage_map.resolve(),research,a.work_dir,receiver,mc,provenance_path=a.provenance.resolve() if a.provenance else None,thermal_points=a.thermal_points)
 cycle=plantmod.cycle_for(selected,rc=a.carbon_return,rh=a.hydrogen_return);plant=plantmod.best_rate(selected,cycle,a.radius)
 checks+=result_checks(selected,cycle,plant);reference={}
 if a.reference:
  ref=read(a.reference);ref=ref['top_detailed'][0] if 'top_detailed' in ref else ref
  rp=ref.get('plant',ref.get('capacity_plant'));rc=ref['cycle'];same=plantmod.station(selected,cycle,rp['helium_per_s'],a.radius)
  reference={'path':str(a.reference),'sha256':sha(a.reference),'same_He_rate_plant':same,
   'capacity_relative_changes':{k:plant[k]/rp[k]-1 for k in ['net_W','mass_kg','net_W_kg','helium_per_s','fresh_C_H_mass_ratio'] if k in plant and k in rp},
   'cycle_relative_changes':{k:cycle[k]/rc[k]-1 for k in ['U_prepared_MeV_per_He','H_nonnu_MeV_per_He','fresh_C_H_mass_ratio','gross_fuel_kg_per_He','gross_H_output_atoms_per_He']},
   'same_rate_relative_changes':{k:same[k]/rp[k]-1 for k in ['net_W','mass_kg','net_W_kg'] if k in same}}
 for p,h in codehash.items():
  if sha(p)!=h:raise RuntimeError('A consumed implementation changed during integration')
 for c in checks:c['passed']=bool(c['passed'])
 if not all(c['passed'] for c in checks):raise ArithmeticError('Variant integration check failed')
 out={'schema':1,'family':a.family,'scope':__doc__,'stage_inputs':archive,'cycle':cycle,'capacity_plant':plant,'reference_comparison':reference,'parameters':{'radius_m':a.radius,'receiver_reference_radius_m':1e9,'carbon_eventual_return':a.carbon_return,'hydrogen_eventual_return':a.hydrogen_return,'thermal_points':a.thermal_points},'checks':checks,'source_sha256':sources,'code_sha256':codehash,
 'qualification':'Conditional integrated accounting for these explicitly named pulse states. Finite grey-opacity spherical hydro, source assembly, receiver thermal/material hardware and global support assumptions remain separately qualified. Surface-temperature proxies are reported diagnostics, not imposed apparatus limits.'}
 write(a.output,out);print(json.dumps({'family':a.family,'checks':len(checks),'plant':{k:plant.get(k) for k in ['feasible','net_W_kg','net_W','mass_kg','helium_per_s','sequential_occupancy_mass999','radiator_area_fraction','fresh_C_H_mass_ratio']},'cycle_prepared_MeV_He':cycle['U_prepared_MeV_per_He'],'comparison':reference.get('capacity_relative_changes'),'output':str(a.output)},indent=2))
if __name__=='__main__':main()
