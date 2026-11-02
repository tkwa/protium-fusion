"""Uniform fuel-mixture sensitivity using the frozen corrected radial solver.

Each invocation runs in its own process. The only solver override is its
initial-composition function. No archived source or result is modified.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT)]
import radial_completion as radial
import pulse_model as pm
from fast_eos import FastEOS

CATS={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
PROFILES={s:ROOT/f'inputs/pulse-profiles/scaled-hot-{s}-64-R30-initial.json' for s in ['A','B1','C13']}
PROFILES['B2']=ROOT/'b2-completion/B2-64-compression-cap.initial.json'
BASELINES=json.loads((ROOT/'open-loop-redesign/corrected-stage-map.json').read_text())

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def tag(x): return format(x,'.8g').replace('.','p').replace('-','m')
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage',choices=list(CATS),required=True)
    ap.add_argument('--h-per-catalyst',type=float,required=True)
    ap.add_argument('--constraint',choices=['fixed_flow','fixed_energy'],default='fixed_energy')
    ap.add_argument('--zones',type=int,default=64)
    ap.add_argument('--opacity',type=float,default=.001)
    ap.add_argument('--rates',choices=['frozen','carbon-median'],default='frozen')
    ap.add_argument('--duration-factor',type=float,default=1.)
    ap.add_argument('--rtol',type=float,default=1e-7)
    ap.add_argument('--temperature-factor',type=float,default=1.)
    ap.add_argument('--prepared-energy-factor',type=float,default=1.)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.h_per_catalyst<=0: raise ValueError('Positive hydrogen ratio required')
    src=PROFILES[a.stage]; baseline_path=ROOT/BASELINES[CATS[a.stage]]
    p=json.loads(src.read_text()); b=json.loads(baseline_path.read_text())
    base_y,target,_=radial.initial(a.stage)
    y=np.zeros(pm.NN);cat_i=pm.SPECIES.index(CATS[a.stage]);p_i=pm.SPECIES.index('p')
    denom=pm.AA[cat_i]+a.h_per_catalyst
    y[cat_i]=1/denom;y[p_i]=a.h_per_catalyst/denom
    if abs(pm.AA@y-1)>1e-14:raise ValueError('Composition baryons do not sum to one')
    old_m=np.array(p['mass_fraction_of_total']);old_r=np.r_[0.,p['outer_radius_m']]
    old_v=np.r_[0.,p['velocity_m_s']];old_T=np.array(p['T_K'])
    if a.zones!=len(old_m):
        q_old=np.r_[0.,np.cumsum(old_m)];q_old[-1]=1.
        q_new=np.linspace(0,1,a.zones+1)**3
        r=np.interp(q_new,q_old,old_r**3)**(1/3)
        v=np.interp(r,old_r,old_v)[1:]
        centers=.5*(q_new[1:]+q_new[:-1])
        parent=np.clip(np.searchsorted(q_old,centers,side='right')-1,0,len(old_T)-1)
        T=old_T[parent];m=np.diff(q_new)
        old_nm=np.r_[.5*(old_m[:-1]+old_m[1:]),.5*old_m[-1]]
        nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
        v*=np.sqrt((old_nm@old_v[1:]**2)/(nm@v**2))
    else:
        m=old_m;r=old_r;v=old_v[1:];T=old_T
        nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
    rho=b['mass_kg']*m/(4*np.pi/3*np.diff(r**3))
    eos=FastEOS()
    base_state=eos.state(rho/1000,T,float(pm.ZZ@base_y),float(sum(base_y)))
    T=T*a.temperature_factor
    new_state=eos.state(rho/1000,T,float(pm.ZZ@y),float(sum(y)))
    U0=float(m@base_state.e)*1e-4;U1=float(m@new_state.e)*1e-4
    K0=float(.5*nm@v**2);vscale=1.
    if a.constraint=='fixed_energy':
        if U1>=(K0+U0)*a.prepared_energy_factor:raise ValueError('Thermal energy exceeds available prepared energy')
        vscale=float(np.sqrt(((K0+U0)*a.prepared_energy_factor-U1)/K0));v=v*vscale
    elif a.prepared_energy_factor!=1.:raise ValueError('Energy factor requires fixed_energy comparison')
    p={'outer_radius_m':r[1:].tolist(),'velocity_m_s':v.tolist(),
       'T_K':T.tolist(),'mass_fraction_of_total':m.tolist(),
       'description':f'Mixture sensitivity {a.stage}, H/catalyst={a.h_per_catalyst}; {a.constraint}; fixed pulse mass, density and temperature profiles.'}
    initial_Y={n:float(vv) for n,vv in zip(pm.SPECIES,y) if vv!=0}
    # A process-local composition override leaves all frozen solver physics intact.
    radial.initial=lambda stage:(y.copy(),target,float(y[cat_i]))
    center=eos.state(b['rho0_kg_m3']/1000,b['T0_K'],float(pm.ZZ@y),float(sum(y)))
    sound_cm_s=float(np.sqrt(center.gamma1*center.p/(b['rho0_kg_m3']/1000)))
    tunit=b['R0_m']*100/sound_cm_s
    duration=b['elapsed_s']*a.duration_factor
    a.out.mkdir(parents=True,exist_ok=True)
    name=f'{a.stage}-H{tag(a.h_per_catalyst)}-{a.constraint}-z{a.zones}-k{tag(a.opacity)}-t{tag(a.duration_factor)}'
    if a.temperature_factor!=1.:name+=f'-T{tag(a.temperature_factor)}'
    if a.prepared_energy_factor!=1.:name+=f'-E{tag(a.prepared_energy_factor)}'
    if a.rtol!=1e-7:name+=f'-rtol{tag(a.rtol)}'
    rate_update=None
    if a.rates=='carbon-median':
        import updated_network
        rate_update=updated_network.install(radial,quantile='median')
        name+='-rates-carbon-median'
    stem=a.out/name
    inputs={'stage':a.stage,'hydrogen_per_catalyst':a.h_per_catalyst,'initial_Y':initial_Y,
      'initial_counts':{CATS[a.stage]:1.,'p':a.h_per_catalyst},'constraint':a.constraint,
      'held_fixed':['total pulse mass','enclosed-mass density profile','temperature profile','physical final time unless temperature guard'],
      'velocity_scale':vscale,'temperature_factor':a.temperature_factor,
      'prepared_energy_factor':a.prepared_energy_factor,'reference_specific_internal_J_kg':U0,
      'reference_specific_kinetic_J_kg':K0,'new_specific_internal_J_kg':U1,
      'specific_prepared_energy_J_kg':U1+float(.5*nm@v**2),
      'initial_max_pressure_Pa':float(np.max(new_state.p))*.1,
      'initial_mass_weighted_pressure_Pa':float(m@new_state.p)*.1,
      'initial_electrons_per_baryon':float(pm.ZZ@y),'initial_ions_per_baryon':float(sum(y)),
      'duration_s':duration,'zones':a.zones,'opacity_m2_kg':a.opacity,
      'profile_source':str(src.relative_to(ROOT)),'profile_source_sha256':sha(src),
      'baseline_source':str(baseline_path.relative_to(ROOT)),'baseline_source_sha256':sha(baseline_path),
      'solver_source_sha256':sha(ROOT/'b2-completion/radial_completion.py'),
      'network_source_sha256':sha(ROOT/'pulse_model.py'),
      'runner_source_sha256':sha(Path(__file__)),
      'scope':'Full 27-isotope spherical hydro sensitivity; prescribed hot inflow; constant grey opacity. Source preparation, 3D merger and apparatus not simulated.'}
    if rate_update is not None:inputs['rate_update']=rate_update
    if a.temperature_factor!=1.:inputs['held_fixed'].remove('temperature profile')
    if a.constraint=='fixed_flow':inputs['held_fixed'].append('radial velocity profile')
    else:inputs['held_fixed'].append('total prepared energy at stated reference multiplier')
    stem.with_suffix('.initial.json').write_text(json.dumps(p,indent=2)+'\n')
    stem.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    result=radial.simulate(a.stage,b['T0_K'],b['rho0_kg_m3'],b['R0_m'],zones=a.zones,
      nuclear=True,reduced_a=False,radiation=True,kappa=a.opacity,initial_profile=p,
      final_time_units=duration/tunit,rtol=a.rtol,atol=1e-12,viscosity_mode='compression_cap',
      save_trace=stem.with_suffix('.npz'),trace_samples=400)
    result['initial_Y']=initial_Y;result['mixture_inputs']=inputs
    if rate_update is not None:result['rate_update']=rate_update
    result['trace_path']=str(stem.with_suffix('.npz'))
    stem.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'case':name,**{k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','peak_temperature_K','wall_seconds']}}),flush=True)
if __name__=='__main__':main()
