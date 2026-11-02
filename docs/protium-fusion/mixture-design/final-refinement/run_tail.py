"""Authorized B2 H2 256-cell late-outflow control: same initial state, 2x duration."""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
P=Path(__file__).resolve().parent;R=P.parent.parent;sys.path[:0]=[str(R),str(R/'b2-completion'),str(R/'mixture-design')]
import radial_completion as radial
import pulse_model as pm
from fast_eos import FastEOS
import updated_network
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 src=Path('/home/tkwa/src/protium-fusion/output/mixture-rate-runs/root256/B2-H2-fixed_energy-z256-k0p001-t1-rates-updated-median.json');d=json.loads(src.read_text());inp=json.loads(src.with_suffix('.inputs.json').read_text());profile=json.loads(src.with_suffix('.initial.json').read_text());y=np.array([inp['initial_Y'].get(s,0.) for s in pm.SPECIES]);_,target,_=radial.initial('B2');radial.initial=lambda stage:(y.copy(),target,float(y[pm.SPECIES.index('c12')]))
 duration=2*d['elapsed_s'];inp['duration_s']=duration;inp['physical_duration_s']=duration;inp['duration_factor']=2.;inp['late_outflow_control']={'original_record':str(src),'original_record_sha256':sha(src),'original_initial_sha256':sha(src.with_suffix('.initial.json')),'original_end_time_s':d['elapsed_s'],'requested_end_time_s':duration,'scope':'Exact original initial state rerun with doubled physical duration, to resolve nonmonotone late outflow. A duration sensitivity, not a matched-time mesh or composition result.'}
 meta=updated_network.install(radial,quantile='median');meta['runner_sha256']=sha(__file__);inp['rate_update']=meta
 eos=FastEOS();rho=d['rho0_kg_m3']/1000;e=eos.state(rho,d['T0_K'],float(pm.ZZ@y),float(sum(y)));tu=d['R0_m']*100/math.sqrt(e.gamma1*e.p/rho)
 raw=Path('/home/tkwa/src/protium-fusion/output/final-refinement/B2-H2-z256-tail2');write(raw.with_suffix('.inputs.json'),inp);write(raw.with_suffix('.initial.json'),profile)
 out=radial.simulate('B2',d['T0_K'],d['rho0_kg_m3'],d['R0_m'],zones=256,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=profile,final_time_units=duration/tu,rtol=d['rtol'],atol=d['atol'],art_linear=d['artificial_viscosity_linear'],art_quad=d['artificial_viscosity_quadratic'],viscosity_mode='compression_cap',save_trace=raw.with_suffix('.npz'),trace_samples=400)
 assert abs(out['initial_energy_J']/d['initial_energy_J']-1)<1e-12
 out.update(initial_Y=inp['initial_Y'],mixture_inputs=inp,rate_update=meta,trace_path=str(raw.with_suffix('.npz')),trace_sha256=sha(raw.with_suffix('.npz')),late_outflow_control=inp['late_outflow_control'])
 write(raw.with_suffix('.json'),out)
 for suffix,data in [('.json',out),('.inputs.json',inp),('.initial.json',profile)]:write(P/'tail'/('B2-H2-z256-tail2'+suffix),data)
 print(json.dumps({k:out[k] for k in ['success','termination_reason','elapsed_s','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}),flush=True)
if __name__=='__main__':main()
