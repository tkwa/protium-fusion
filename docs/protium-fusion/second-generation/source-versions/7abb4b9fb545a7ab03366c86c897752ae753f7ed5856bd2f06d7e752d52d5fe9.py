#!/usr/bin/env python3
"""Paired spherical N15 size/speed controls, using the latest delivered precursor family.
Prescribed hot-flow assembly, 27 isotopes, current carbon medians, gray LTE.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

P=Path(__file__).resolve().parent; F=P.parent
sys.path[:0]=[str(F/'b2-completion'),str(F),str(F/'mixture-design')]
import pulse_model as pm
import radial_completion as radial
from fast_eos import FastEOS
import updated_network

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--mass-factor',type=float,default=100.)
    ap.add_argument('--speed-factor',type=float,default=1.);ap.add_argument('--zones',type=int,default=32)
    ap.add_argument('--time-factor',type=float,default=1.);a=ap.parse_args()
    source=P/'probe-input.json';saved=json.loads(source.read_text())
    d=saved['record'];old=saved['initial_profile']
    scale=a.mass_factor**(1/3);n=a.zones;oldx=np.linspace(0,1,len(old['T_K'])+1);x=np.linspace(0,1,n+1)
    oldmid=(oldx[1:]**3+oldx[:-1]**3)/2;mid=(x[1:]**3+x[:-1]**3)/2
    profile=dict(outer_radius_m=(np.interp(x[1:],oldx,np.r_[0.,old['outer_radius_m']])*scale).tolist(),
        velocity_m_s=(np.interp(x[1:],oldx,np.r_[0.,old['velocity_m_s']])*a.speed_factor).tolist(),
        T_K=np.interp(mid,oldmid,old['T_K']).tolist(),
        description='Delivered B1 hot source, radius scaled at fixed density; temperatures interpolated in enclosed mass, velocities/radii at boundary nodes. Coarse paired control.')
    y=np.array([d['initial_Y'].get(n,0.) for n in pm.SPECIES]);radial.initial=lambda stage:(y.copy(),'c12',float(y[pm.SPECIES.index('n15')]))
    meta=updated_network.install(radial)
    eos=FastEOS();rho=d['rho0_kg_m3']/1000.;center=eos.state(rho,d['T0_K'],float(pm.ZZ@y),float(sum(y)))
    R=d['R0_m']*scale;tunit=R*100/np.sqrt(center.gamma1*center.p/rho)
    duration=d['elapsed_s']*scale*a.time_factor
    stem=P/(f"pulse-n15-m{a.mass_factor:g}-v{a.speed_factor:g}-z{n}-t{a.time_factor:g}".replace('.','p'))
    trace_dir=P.parents[2]/'output/protium-second-generation-2026-11-08';trace_dir.mkdir(parents=True,exist_ok=True)
    trace=trace_dir/(stem.name+'.npz')
    result=radial.simulate('B1',d['T0_K'],d['rho0_kg_m3'],R,zones=n,nuclear=True,
        reduced_a=False,radiation=True,kappa=.001,initial_profile=profile,
        final_time_units=duration/tunit,rtol=1e-6,atol=1e-11,
        viscosity_mode='compression_cap',save_trace=trace,trace_samples=400)
    result.update(initial_Y=d['initial_Y'],rate_update=meta,probe=vars(a),
        input_source=str(source),input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        qualification=__doc__,trace_path=str(trace),trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest())
    stem.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy',
        'desired_product_per_initial_catalyst','initial_energy_J','radiation_escaped_J','energy_residual_fraction_initial','wall_seconds']}),flush=True)

if __name__=='__main__':main()
