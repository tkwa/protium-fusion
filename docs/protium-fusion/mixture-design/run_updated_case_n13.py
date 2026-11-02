"""Repeat an exact saved mixture profile with explicit runtime carbon-rate update."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion') if not (HERE.parent/'pulse_model.py').exists() else HERE.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT),str(HERE)]
import radial_completion as radial
import pulse_model as pm
from fast_eos import FastEOS
import updated_network_n13 as update
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
CATS={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--quantile',default='median',choices=['median','minus1sigma','plus1sigma','minus2sigma','plus2sigma','minus3sigma','plus3sigma'])
    ap.add_argument('--mode',choices=['updated','identity'],default='updated');ap.add_argument('--interpolation',choices=['pchip','linear'],default='pchip')
    a=ap.parse_args();src=a.record.resolve();b=json.loads(src.read_text());p=json.loads(src.with_suffix('.initial.json').read_text())
    inp=json.loads(src.with_suffix('.inputs.json').read_text());stage=b['stage'];y=np.array([inp['initial_Y'].get(n,0.) for n in pm.SPECIES])
    assert abs(pm.AA@y-1)<1e-12
    _,target,_=radial.initial(stage);cat=y[pm.SPECIES.index(CATS[stage])]
    radial.initial=lambda s:(y.copy(),target,float(cat))
    meta=update.install(radial,quantile=a.quantile,mode=a.mode,interpolation=a.interpolation)
    meta['baseline_record']=str(src);meta['baseline_record_sha256']=sha(src);meta['initial_profile_sha256']=sha(src.with_suffix('.initial.json'))
    meta['runner_sha256']=sha(__file__)
    inp=dict(inp);inp['rate_update']=meta;inp['runner_source_sha256']=sha(__file__)
    duration=inp.get('duration_s',inp.get('physical_duration_s')); assert duration is not None
    if 'duration_s' in inp and 'physical_duration_s' in inp: assert abs(inp['duration_s']/inp['physical_duration_s']-1)<1e-12
    inp['duration_s']=duration
    assert b.get('viscosity_mode','compression_cap')=='compression_cap'
    eos=FastEOS();rho=b['rho0_kg_m3']/1000
    center=eos.state(rho,b['T0_K'],float(pm.ZZ@y),float(sum(y)))
    tunit=b['R0_m']*100/float(np.sqrt(center.gamma1*center.p/rho))
    a.out.mkdir(parents=True,exist_ok=True)
    suffix=f'-rates-carbon-li2006-{a.mode}-{a.quantile}'
    if a.interpolation!='pchip':suffix+='-'+a.interpolation
    stem=a.out/(src.stem+suffix)
    stem.with_suffix('.initial.json').write_text(json.dumps(p,indent=2)+'\n')
    stem.with_suffix('.inputs.json').write_text(json.dumps(inp,indent=2)+'\n')
    result=radial.simulate(stage,b['T0_K'],b['rho0_kg_m3'],b['R0_m'],zones=b['zones'],nuclear=True,reduced_a=False,
      radiation=b['radiation_diffusion'],kappa=b['opacity_m2_kg'],initial_profile=p,final_time_units=duration/tunit,
      rtol=b['rtol'],atol=b['atol'],viscosity_mode='compression_cap',save_trace=stem.with_suffix('.npz'),trace_samples=400)
    result.update(initial_Y=inp['initial_Y'],mixture_inputs=inp,rate_update=meta,trace_path=str(stem.with_suffix('.npz')))
    stem.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    keys=['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','peak_temperature_K','wall_seconds']
    print(json.dumps({'case':stem.name,**{k:result[k] for k in keys}}),flush=True)
if __name__=='__main__':main()
