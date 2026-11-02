"""Bounded scale discriminator, using existing hot inflows and unmodified radial code."""
import argparse,hashlib,json,sys,time,traceback
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['A','B1','B2','C13','He']);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--rtol',type=float,default=1e-6);ap.add_argument('--tag',default='');ap.add_argument('--no-radiation',action='store_true');ap.add_argument('--radius',type=float,default=30.);ap.add_argument('--rate-floor',type=float)
a=ap.parse_args();sys.path.insert(0,str(a.root));from radial_pulse import simulate
import radial_pulse as rp
if a.rate_floor:
 original_network=rp.network
 rp.network=lambda temp,rho,y:original_network(np.maximum(temp,a.rate_floor),rho,y)
stage=a.stage;nz=64;radius=a.radius;a.output.mkdir(parents=True,exist_ok=True)
source=a.root/(f'compression-{stage}-{nz}'+('-full' if stage=='A' else '')+'-burning.json');initial_source=a.root/f'compression-{stage}-{nz}-initial.json'
old=json.loads(source.read_text());profile=json.loads(initial_source.read_text());old_r=old['R0_m']
profile['outer_radius_m']=(np.array(profile['outer_radius_m'])*radius/old_r).tolist()
profile['description']=f'Spatially scaled prescribed hot {stage} inflow, reference radius {old_r:g} to {radius:g} m; density, temperature, velocity unchanged. No preparation model.'
stem=f'scaled-hot-{stage}-{nz}-R{radius:g}'+a.tag;(a.output/(stem+'-initial.json')).write_text(json.dumps(profile,indent=2)+'\n')
inputs={'stage':stage,'source_radius_m':old_r,'radius_m':radius,'zones':nz,'run_options':{'rtol':a.rtol,'atol':1e-11,'radiation':not a.no_radiation,'rate_floor_K':a.rate_floor,'final_time_units':4},'source_hashes':{str(p.relative_to(a.root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,initial_source,a.root/'radial_pulse.py',a.root/'pulse_model.py',a.root/'fast_eos.py']}}
(a.output/(stem+'-inputs.json')).write_text(json.dumps(inputs,indent=2)+'\n')
start=time.monotonic()
try:
 r=simulate(stage,old['T0_K'],old['rho0_kg_m3'],radius,zones=nz,nuclear=True,reduced_a=False,radiation=not a.no_radiation,kappa=old['opacity_m2_kg'],initial_profile=profile,final_time_units=4,save_trace=a.output/(stem+'.npz'),rtol=a.rtol,atol=1e-11,carbon13_fraction=old.get('initial_carbon13_number_fraction',0.))
 r['network_rate_evaluation_temperature_floor_K']=a.rate_floor
 r['effective_termination_reason']=r['termination_reason'] if r['success'] else 'solver_failure'
 r['diagnostic_scope']='Scale discriminator with prescribed hot inflow; no driver, multidimensional validation or grid convergence at this scale.'
 (a.output/(stem+'.json')).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
 fields=['stage','mass_kg','initial_energy_J','nuclear_energy_J','weak_energy_escaped_upper_J','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','elapsed_s','energy_residual_fraction_initial','success','termination_reason','wall_seconds','minimum_mass_fraction','max_cell_baryon_residual','output_mass_fractions']
 print(json.dumps({k:r[k] for k in fields},indent=2),flush=True)
except Exception as e:
 failure={'stage':stage,'error':repr(e),'wall_seconds':time.monotonic()-start,'traceback':traceback.format_exc()}
 (a.output/(stem+'-failure.json')).write_text(json.dumps(failure,indent=2)+'\n');print(json.dumps(failure,indent=2),flush=True);raise
