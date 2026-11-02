"""Finite-length and free-space controls for the prescribed electron process."""
from pathlib import Path
import argparse,json,gzip,time
import numpy as np
from lepton_escape_mc import (Cylinder,local_state,prepare_species,ensemble,summary,
    EPS0,ME,C,EV,pi)
ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,default=Path(__file__).parent)
ap.add_argument('--outdir',type=Path,default=Path(__file__).parent)
args=ap.parse_args()
ROOT=args.root;OUT=args.outdir
plant=next(p for p in json.loads((ROOT/'network-plant.json').read_text())['cases']
    if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
c=Cylinder(plant['selected_source_profiles']['22049'])
cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])
r=float(np.interp(.5,cdf,c.eq['x'])*c.scale);local=local_state(c,r)
table=prepare_species(local['species'],15.)
tau=3*pi*EPS0*ME**3*C**3/(EV**4*c.B**2)
columns=['time_s','radiated_eV','collision_eV','final_eV','energy_time_eV_s',
    'perp_energy_time_eV_s','parallel_energy_time_eV_s','occupation_time_s',
    'radial_force_displacement_m','steps','escaped','maximum_eV','maximum_step_s','final_z_m','negative_energy_proposals','initial_vz_m_s','final_vz_m_s']
for name,factor,n,seed,step,tmax,closed,burn in [
    ('half',.5,2048,19371,.02,.1,False,0.),
    ('double',2.,2048,19379,.02,.4,False,0.),
    ('closed-fine',1.,2048,28537,.01,.006,True,.001)]:
    print('START',name,flush=True);start=time.monotonic()
    L=factor*plant['active_length_m']
    rows,hist,pitch=ensemble(n,seed,table,local['relative_azimuthal_m_s'],tau,L,tmax,step,1.,burn,not closed,c.B)
    out=summary(rows,hist,pitch,1.,tau)
    out.update({'scope':__doc__,'reaction_fraction_inside':.5,'radius_m':r,'log_lambda':15.,
        'birth_energy_eV':1.,'step_fraction':step,'seed':seed,'time_limit_s':tmax,'burn_in_s':burn,
        'absorbing_ends':not closed,'runtime_s':time.monotonic()-start,'active_length_m':L,
        'length_factor':factor,'local_state':{k:v for k,v in local.items() if k!='species'}})
    if closed:
        estimates=rows[:,13]**2/(2*(tmax-burn))
        out['finite_time_axial_diffusion_m2_s']={'mean':float(estimates.mean()),'standard_error':float(estimates.std(ddof=1)/np.sqrt(n))}
    (OUT/f'lepton-escape-{name}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    payload=json.dumps({'columns':columns,'rows':rows.tolist()},separators=(',',':'),allow_nan=False).encode()
    (OUT/f'lepton-escape-{name}-trajectories.json.gz').write_bytes(gzip.compress(payload,mtime=0))
    print('DONE',name,{k:out[k] for k in ['runtime_s','residence_s','radiated_eV','occupation_mean_energy_eV','negative_radial_energy_proposals','censored']},flush=True)
