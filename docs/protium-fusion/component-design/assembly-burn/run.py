"""Hot nuclear continuation of an exact delivered B1 preparation endpoint."""
from pathlib import Path
import argparse,json,sys,math,hashlib
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent.parent;sys.path[:0]=[str(ROOT),str(ROOT/'b2-completion'),str(ROOT/'mixture-design')]
import radial_completion as radial
import pulse_model as pm
from fast_eos import FastEOS
import updated_network
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--label',choices=['k0','k0p001','k0p01'],required=True);a=ap.parse_args();stem='B1-delivered-'+a.label;ip=P/'inputs'/(stem+'.inputs.json');pp=P/'inputs'/(stem+'.initial.json');inp=json.loads(ip.read_text());profile=json.loads(pp.read_text());y=np.array([inp['initial_Y'].get(s,0.) for s in pm.SPECIES]);_,target,_=radial.initial('B1');radial.initial=lambda stage:(y.copy(),target,float(y[pm.SPECIES.index('n15')]))
 meta=updated_network.install(radial,quantile='median');meta['baseline_record']=inp['preparation_record'];meta['baseline_record_sha256']=inp['preparation_record_sha256'];meta['initial_profile_sha256']=sha(pp);meta['runner_sha256']=sha(__file__)
 eos=FastEOS();rho=inp['rho0_kg_m3']/1000;center=eos.state(rho,inp['T0_K'],float(pm.ZZ@y),float(sum(y)));tu=inp['R0_m']*100/math.sqrt(center.gamma1*center.p/rho)
 raw=Path('/home/tkwa/src/protium-fusion/output/assembly-burn')/stem;inp['rate_update']=meta;inp['runner_sha256']=sha(__file__);inp['bridge_sha256']=sha(P/'bridge.json');write(raw.with_suffix('.inputs.json'),inp);write(raw.with_suffix('.initial.json'),profile)
 result=radial.simulate('B1',inp['T0_K'],inp['rho0_kg_m3'],inp['R0_m'],zones=128,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=profile,final_time_units=inp['duration_s']/tu,rtol=inp['rtol'],atol=inp['atol'],art_linear=.5,art_quad=2.,viscosity_mode='compression_cap',save_trace=raw.with_suffix('.npz'),trace_samples=400)
 assert abs(result['initial_energy_J']/inp['burn_start_energy_J']-1)<1e-12
 dep=result['nuclear_energy_J']-result['weak_energy_escaped_upper_J'];result.update(initial_Y=inp['initial_Y'],mixture_inputs=inp,rate_update=meta,trace_path=str(raw.with_suffix('.npz')),trace_sha256=sha(raw.with_suffix('.npz')),assembly_energy_ledger={'preparation_initial_FastEOS_J':inp['preparation_initial_energy_J'],'preparation_escaped_radiation_J':inp['preparation_escaped_radiation_J'],'burn_start_J':result['initial_energy_J'],'burn_deposited_nuclear_J':dep,'gain_over_preparation_FastEOS_energy':dep/inp['preparation_initial_energy_J'],'gain_over_burn_start_energy':result['deposited_nuclear_over_initial_energy'],'combined_energy_residual_fraction_preparation':(result['final_internal_J']+result['final_kinetic_J']+result['weak_energy_escaped_upper_J']+result['radiation_escaped_J']+inp['preparation_escaped_radiation_J']-inp['preparation_initial_energy_J']-result['nuclear_energy_J'])/inp['preparation_initial_energy_J'],'scope':'Preescaped radiation is kept as energy leaving preparation; no receiver recovery or source physical attainability is credited. FEOS handoff discrepancy remains separately archived in bridge.json.'})
 write(raw.with_suffix('.json'),result);write(P/'records'/(stem+'.json'),result);write(P/'records'/(stem+'.inputs.json'),inp);write(P/'records'/(stem+'.initial.json'),profile)
 print(json.dumps({'case':stem,**{k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','peak_temperature_K','wall_seconds']},'outer_gain':dep/inp['preparation_initial_energy_J']}),flush=True)
if __name__=='__main__':main()
