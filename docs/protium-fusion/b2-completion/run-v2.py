"""Completion experiment with conservative mass/thermal profile refinement."""
import argparse,hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--zones',type=int,default=64);ap.add_argument('--source-zones',type=int,default=64)
ap.add_argument('--rtol',type=float,default=1e-7);ap.add_argument('--atol',type=float,default=1e-12)
ap.add_argument('--mode',choices=['legacy','volume_gate','divergence'],default='divergence')
ap.add_argument('--name',default='B2-64-divergence');ap.add_argument('--time',type=float,default=4)
ap.add_argument('--no-radiation',action='store_true');ap.add_argument('--energy-match',action='store_true')
ap.add_argument('--trace-samples',type=int,default=240)
a=ap.parse_args();sys.path.insert(0,str(a.root));sys.path.insert(0,str(Path(__file__).parent))
from radial_completion import simulate
source=a.root/f'compression-B2-{a.source_zones}-initial.json'
old=json.loads(source.read_text());initial=old.copy();initial['outer_radius_m']=[r*.1 for r in old['outer_radius_m']]
oldr=np.r_[0.,initial['outer_radius_m']];oldv=np.r_[0.,old['velocity_m_s']]
oldq=np.linspace(0,1,a.source_zones+1)**3;oldm=np.diff(oldq)
oldnode=np.r_[.5*(oldm[:-1]+oldm[1:]),.5*oldm[-1]]
newq=np.linspace(0,1,a.zones+1)**3;newm=np.diff(newq)
newnode=np.r_[.5*(newm[:-1]+newm[1:]),.5*newm[-1]]
if a.zones != a.source_zones:
 if a.zones%a.source_zones:raise ValueError('Only nested integer refinements are supported')
 newr=np.interp(newq,oldq,oldr**3)**(1/3)
 parent=np.searchsorted(oldq,.5*(newq[:-1]+newq[1:]),side='right')-1
 initial['outer_radius_m']=newr[1:].tolist()
 initial['T_K']=np.asarray(old['T_K'])[parent].tolist()
 initial['velocity_m_s']=np.interp(newr,oldr,oldv)[1:].tolist()
 initial['mass_fraction_of_total']=newm.tolist()
 initial.pop('rho_kg_m3',None);initial.pop('desired_product_per_initial_catalyst',None)
oldK=.5*oldnode@oldv[1:]**2;newK=.5*newnode@np.asarray(initial['velocity_m_s'])**2
velocity_scale=float(np.sqrt(oldK/newK)) if a.energy_match else 1.
initial['velocity_m_s']=(np.asarray(initial['velocity_m_s'])*velocity_scale).tolist()
initial['description']='Original64 B2 hot inflow scaled0.1 in radius; nested mass-volume refinement with inherited parent temperature and radial-linear nodal velocity.'
regrid={'source_zones':a.source_zones,'zones':a.zones,'density_conservative':True,'parent_temperature_preserved':True,'velocity_interpolation':'piecewise linear in radius','specific_kinetic_before_J_kg':float(oldK),'specific_kinetic_after_unscaled_J_kg':float(newK),'relative_kinetic_quadrature_change':float(newK/oldK-1),'velocity_energy_match_scale':velocity_scale,'energy_match_requested':a.energy_match}
a.output.mkdir(parents=True,exist_ok=True);stem=a.output/a.name
initial_path=stem.with_suffix('.initial.json');initial_path.write_text(json.dumps(initial,indent=2)+'\n')
sources=[Path(__file__),Path(__file__).with_name('radial_completion.py'),a.root/'radial_pulse.py',a.root/'fast_eos.py',a.root/'electron-eos.npz',a.root/'pulse_model.py',source]
config={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()};config['source_hashes']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};config['regridding']=regrid;stem.with_suffix('.inputs.json').write_text(json.dumps(config,indent=2)+'\n')
try:
 r=simulate('B2',3e8,1e8,30,zones=a.zones,nuclear=True,reduced_a=False,radiation=not a.no_radiation,initial_profile=initial,final_time_units=a.time,rtol=a.rtol,atol=a.atol,viscosity_mode=a.mode,save_trace=stem.with_suffix('.npz'),trace_samples=a.trace_samples)
 r['diagnostic_scope']='Completion of reduced spherical specified-hot-inflow model. No multidimensional driver or plant validation.'
 r['initial_profile_path']=str(initial_path);r['initial_profile_sha256']=hashlib.sha256(initial_path.read_bytes()).hexdigest();r['regridding']=regrid
 stem.with_suffix('.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:r[k] for k in ['success','termination_reason','elapsed_sound_times','elapsed_s','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','energy_residual_fraction_initial','minimum_mass_fraction','wall_seconds','nfev','njev','integrated_nuclear_energy_J','nuclear_energy_J']},indent=2),flush=True)
except Exception as e:
 stem.with_suffix('.failure.json').write_text(json.dumps({'error':repr(e),'traceback':traceback.format_exc()},indent=2)+'\n');raise
