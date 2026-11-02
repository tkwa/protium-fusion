"""Completion experiment; original model modules remain frozen."""
import argparse,hashlib,json,sys,time,traceback
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--zones',type=int,default=64);ap.add_argument('--rtol',type=float,default=1e-7);ap.add_argument('--atol',type=float,default=1e-12);ap.add_argument('--mode',choices=['legacy','volume_gate','divergence'],default='volume_gate');ap.add_argument('--name',default='B2-64-volume-gate');ap.add_argument('--time',type=float,default=4);ap.add_argument('--no-radiation',action='store_true');a=ap.parse_args();sys.path.insert(0,str(a.root));sys.path.insert(0,str(Path(__file__).parent))
from radial_completion import simulate
initial=json.loads((a.root/f'compression-B2-{a.zones}-initial.json').read_text());initial['outer_radius_m']=[r*.1 for r in initial['outer_radius_m']];initial['description']='Original B2 prescribed hot inflow, radii scaled by0.1 to30m reference; new compression-only scalar viscosity, no driver model.'
a.output.mkdir(parents=True,exist_ok=True);stem=a.output/a.name
sources=[Path(__file__),Path(__file__).with_name('radial_completion.py'),a.root/'radial_pulse.py',a.root/'fast_eos.py',a.root/'electron-eos.npz',a.root/'pulse_model.py',a.root/f'compression-B2-{a.zones}-initial.json']
config={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()};config['source_hashes']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};stem.with_suffix('.inputs.json').write_text(json.dumps(config,indent=2)+'\n')
try:
 r=simulate('B2',3e8,1e8,30,zones=a.zones,nuclear=True,reduced_a=False,radiation=not a.no_radiation,initial_profile=initial,final_time_units=a.time,rtol=a.rtol,atol=a.atol,viscosity_mode=a.mode,save_trace=stem.with_suffix('.npz'))
 r['diagnostic_scope']='Completion of the reduced spherical specified-hot-inflow calculation. No multidimensional driver or plant validation.'
 stem.with_suffix('.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:r[k] for k in ['success','termination_reason','elapsed_sound_times','elapsed_s','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','energy_residual_fraction_initial','minimum_mass_fraction','wall_seconds','nfev','njev','integrated_nuclear_energy_J','nuclear_energy_J']},indent=2),flush=True)
except Exception as e:
 stem.with_suffix('.failure.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()},indent=2)+'\n');raise
