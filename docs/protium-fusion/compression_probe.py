"""Controlled convergent-flow experiment derived from an adiabatic expansion.

The velocity-reversed profile is a specified initial condition, not a demonstrated
driver. Its thermal and kinetic energy both count as input. Compare nonreacting
return with nuclear burning during the same inflow. Numerical viscosity makes
the discrete evolution imperfectly reversible; resolution matters.
"""
import argparse
import json
from pathlib import Path
from radial_pulse import simulate

ROOT=Path(__file__).parent
ap=argparse.ArgumentParser();ap.add_argument('--zones',type=int,default=128)
ap.add_argument('--stage',default='A',choices=['A','B1','B2','C13','He'])
ap.add_argument('--full-network',action='store_true')
ap.add_argument('--temp',type=float);ap.add_argument('--rho',type=float)
ap.add_argument('--radius',type=float);ap.add_argument('--tag',default='')
ap.add_argument('--carbon13',type=float,default=0.)
a=ap.parse_args();nz=a.zones
reduced=a.stage=='A' and not a.full_network
params={'A':(3e8,1e8,300),'B1':(2e8,1e6,1500),'B2':(3e8,1e8,300),
        'C13':(3e8,1e8,100),'He':(1e9,1e10,100)}[a.stage]
params=tuple(new if new is not None else old for new,old in zip([a.temp,a.rho,a.radius],params))
stem=f'compression-{a.stage}-{nz}'+('-full' if a.full_network else '')
if a.tag: stem+='-'+a.tag

def save(label,data):
    (ROOT/f'{stem}-{label}.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(label,{k:data.get(k) for k in ['nuclear_energy_J','initial_energy_J',
        'desired_product_per_initial_catalyst','peak_temperature_K',
        'energy_residual_fraction_initial','termination_reason','wall_seconds']},flush=True)

ref=simulate(a.stage,*params,zones=nz,nuclear=False,reduced_a=reduced,carbon13_fraction=a.carbon13)
save('reference',ref)
profile=ref['profile'].copy();profile['velocity_m_s']=[-v for v in profile['velocity_m_s']]
profile['description']=f'Velocity reversal of nonreacting expansion after two sound times, {nz} cells; no modeled preparation.'
(ROOT/f'{stem}-initial.json').write_text(json.dumps(profile,indent=2)+'\n')
back=simulate(a.stage,*params,zones=nz,nuclear=False,reduced_a=reduced,
              initial_profile=profile,save_trace=ROOT/f'{stem}-nonreacting.npz',carbon13_fraction=a.carbon13)
save('nonreacting',back)
burn=simulate(a.stage,*params,zones=nz,nuclear=True,reduced_a=reduced,radiation=True,
              initial_profile=profile,final_time_units=4,save_trace=ROOT/f'{stem}-burning.npz',carbon13_fraction=a.carbon13)
save('burning',burn)
