"""Independent Lagrangian reference for the first fixed-ion, trapped-radiation model.

Runs the previously repaired radial solver, with two explicit physical changes:
freeze Yi, including its temperature-equation derivative; omit photon diffusion.
The finite-volume Athena++ calculations use the same underlying EOS and rates.
"""
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser()
p.add_argument('--root',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--zones',type=int,nargs='+',default=[64,128,256,512])
p.add_argument('--yi',type=float,default=2/15)
p.add_argument('--physical-ions',action='store_true')
p.add_argument('--radiation',action='store_true')
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.root))
original=a.root/'b2-completion/radial_completion.py'
src=original.read_text()
if not a.physical_ions:
    old='-de_dyi*np.sum(dy,axis=1)-de_dye*(dy@ZZ)'
    assert old in src
    src=src.replace(old,'-de_dye*(dy@ZZ)')
adapted=a.out/'radial_fixed_ions.py';adapted.write_text(src)
spec=importlib.util.spec_from_file_location('radial_fixed_ions',adapted)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
if not a.physical_ions:
    class FrozenEOS(mod.FastEOS):
        def state(self,rho,T,ye,yi):
            return super().state(rho,T,ye,np.zeros_like(np.asarray(yi))+a.yi)
    mod.FastEOS=FrozenEOS
source=a.root/'reassessment-2049/small-stages/scaled-hot-A-64-R30-initial.json'
old=json.loads(source.read_text());oldr=np.r_[0.,old['outer_radius_m']];oldv=np.r_[0.,old['velocity_m_s']]
oldq=np.linspace(0,1,65)**3
for n in a.zones:
    q=np.linspace(0,1,n+1)**3
    rr=np.interp(q,oldq,oldr**3)**(1/3)
    parent=np.searchsorted(oldq,.5*(q[:-1]+q[1:]),side='right')-1
    initial=dict(outer_radius_m=rr[1:].tolist(),velocity_m_s=np.interp(rr,oldr,oldv)[1:].tolist(),
                 T_K=np.asarray(old['T_K'])[parent].tolist(),description='Nested refinement of original A64 hot profile; no discrete energy rescaling')
    begin=time.monotonic()
    r=mod.simulate('A',3e8,1e8,30,zones=n,nuclear=True,reduced_a=True,radiation=a.radiation,
        initial_profile=initial,final_time_units=4,rtol=1e-8,atol=1e-13,
        viscosity_mode='compression_cap',trace_samples=240,save_trace=a.out/f'A-{n}.npz')
    r['experiment']=dict(fixed_ions=not a.physical_ions,yi=a.yi,radiation=a.radiation,
        original_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
        adapted_sha256=hashlib.sha256(adapted.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),wall_s=time.monotonic()-begin)
    (a.out/f'A-{n}.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['zones','success','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','energy_residual_fraction_initial','experiment']}),flush=True)
