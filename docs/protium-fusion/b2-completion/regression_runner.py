"""Four-stage corrected-viscosity regression against frozen small-stage inputs."""
from pathlib import Path
import argparse, hashlib, json, sys, time, traceback

p=argparse.ArgumentParser()
p.add_argument('stage',choices=['A','B1','C13','He'])
p.add_argument('--root',type=Path,required=True)
p.add_argument('--model-directory',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--rtol',type=float,default=1e-7)
p.add_argument('--atol',type=float,default=1e-12)
p.add_argument('--mode',choices=['volume_gate','compression_cap','divergence','legacy'],default='compression_cap')
a=p.parse_args()
sys.path.insert(0,str(a.root));sys.path.insert(0,str(a.model_directory))
from radial_completion import simulate
s=a.stage;radius=3 if s=='He' else 30
old_dir=a.root/'inputs/pulse-profiles'
name=f'scaled-hot-{s}-64-R{radius}'
old_path=old_dir/(name+'.json');initial_path=old_dir/(name+'-initial.json');inputs_path=old_dir/(name+'-inputs.json')
old=json.loads(old_path.read_text());initial=json.loads(initial_path.read_text());inputs=json.loads(inputs_path.read_text())
assert inputs['radius_m']==radius and inputs['zones']==64
assert old['R0_m']==radius and old['zones']==64 and old['radiation_diffusion']
a.output.mkdir(parents=True,exist_ok=True)
stem=a.output/f'regression-{s}-{a.mode}'
sources=[Path(__file__),a.model_directory/'radial_completion.py',a.root/'radial_pulse.py',a.root/'pulse_model.py',a.root/'fast_eos.py',a.root/'electron-eos.npz',old_path,initial_path,inputs_path]
configuration=dict(stage=s,radius_m=radius,zones=64,rtol=a.rtol,atol=a.atol,mode=a.mode,final_time_units=4,radiation=True,trace_samples=120,original_profile_loaded_exactly=True,source_hashes={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in sources})
stem.with_suffix('.inputs.json').write_text(json.dumps(configuration,indent=2)+'\n')
t=time.monotonic()
try:
 r=simulate(s,old['T0_K'],old['rho0_kg_m3'],radius,zones=64,nuclear=True,reduced_a=False,radiation=True,kappa=old['opacity_m2_kg'],initial_profile=initial,final_time_units=4,rtol=a.rtol,atol=a.atol,carbon13_fraction=old.get('initial_carbon13_number_fraction',0.),viscosity_mode=a.mode,trace_samples=120,save_trace=stem.with_suffix('.npz'))
 r['diagnostic_scope']='Regression of the saved 64-zone small prescribed hot inflow with corrected compression-only scalar viscosity; full network and radiation retained. No new driver or multidimensional validation.'
 stem.with_suffix('.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
 keys=['stage','success','termination_reason','elapsed_sound_times','elapsed_s','initial_energy_J','nuclear_energy_J','integrated_nuclear_energy_J','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','max_cell_baryon_residual','minimum_mass_fraction','wall_seconds','nfev','njev']
 print(json.dumps({k:r.get(k) for k in keys},indent=2),flush=True)
except Exception as error:
 stem.with_suffix('.failure.json').write_text(json.dumps(dict(error=repr(error),wall_seconds=time.monotonic()-t,traceback=traceback.format_exc()),indent=2)+'\n');raise
