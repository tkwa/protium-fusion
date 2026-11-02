"""Explicit first-contact inflow probe with an unvalidated cold EOS extension.

Cold fully ionized matter is a surrogate, not a model of solid N/H pellets.
This deliberately exposes rather than hides the initial-condition gap.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))


def build_row(task):
    d,temps=task
    import pynucastro as p
    from pynucastro.eos.electron_eos import ElectronEOS
    eos=ElectronEOS();comp=p.Composition(['p']);comp.set_all(1.)
    pressure=[];energy=[]
    for T in temps:
        e,pos=eos.pe_state(d,T,comp,compute_derivs=False,eta_guess_max=1e10)
        pressure.append(e.p+pos.p)
        energy.append(d*(e.e+pos.e))
    return np.log(pressure),np.log(energy)


def build_table():
    original=np.load(ROOT/'electron-eos.npz')
    # Preserve every original hot-table node exactly. Add 24 nodes per decade.
    lt_cold=np.linspace(np.log(100.),original['log_temperature'][0],73)[:-1]
    ld=original['log_density']
    pp=[];ee=[];started=time.monotonic()
    with ProcessPoolExecutor(max_workers=2) as pool:
        for i,(p,e) in enumerate(pool.map(build_row,[(d,np.exp(lt_cold)) for d in np.exp(ld)])):
            pp.append(p);ee.append(e)
            if i%10==0:print('cold EOS rows',i+1,'/',len(ld),'seconds',time.monotonic()-started,flush=True)
    np.savez_compressed(HERE/'electron-eos-cold-surrogate.npz',log_density=ld,
        log_temperature=np.r_[lt_cold,original['log_temperature']],
        log_pressure=np.c_[np.array(pp),original['log_pressure']],
        log_energy_density=np.c_[np.array(ee),original['log_energy_density']],
        source='pynucastro ElectronEOS; fully ionized unvalidated cold-material surrogate; original hot nodes preserved')
    print('table completed',time.monotonic()-started,flush=True)


def profile(zones,temperature,contact_density):
    source=json.loads((ROOT/'compression-A-128-initial.json').read_text())
    summary=json.loads((ROOT/'compression-A-128-burning.json').read_text())
    if 128%zones:raise ValueError('zones must divide 128')
    ids=np.arange(128//zones-1,128,128//zones)
    r=np.array(source['outer_radius_m'])[ids]
    v=-np.array(source['velocity_m_s'])[ids]
    mass=summary['mass_kg']
    m=mass*np.diff(np.linspace(0,1,zones+1)**3)
    node_m=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
    rescale=np.sqrt(summary['initial_kinetic_energy_J']/(.5*np.dot(node_m,v*v)))
    v*=rescale
    def density(dt):return m/(4*np.pi/3*np.diff(np.r_[0.,r+v*dt]**3))
    hi=.01
    while density(hi).max()>contact_density:hi*=2
    back=brentq(lambda dt:density(dt).max()-contact_density,0.,hi,xtol=1e-15)
    rr=r+v*back
    return {'outer_radius_m':rr.tolist(),'velocity_m_s':(-v).tolist(),
        'rho_kg_m3':density(back).tolist(),'T_K':np.full(zones,temperature).tolist(),
        'mass_fraction_of_total':(m/mass).tolist(),
        'ballistic_lookback_s':back,'velocity_rescale_for_kinetic_energy':rescale,
        'description':f'First-contact granular-flow surrogate: maximum density {contact_density} kg/m3, uniform T={temperature} K, remotely deliverable radial velocity schedule; no merger or cold-material EOS validation.'}


def run(zones,temperature,contact_density,nuclear,final_units,art_linear=.5,art_quad=2.):
    from fast_eos import FastEOS
    import radial_pulse as radial
    radial.FastEOS=lambda:FastEOS(HERE/'electron-eos-cold-surrogate.npz')
    p=profile(zones,temperature,contact_density)
    stem=f'cold-driver-A-{zones}-T{temperature:g}-rho{contact_density:g}-'+('burn' if nuclear else 'nonreacting')
    if art_linear != .5 or art_quad != 2.: stem+=f'-visc{art_linear:g}_{art_quad:g}'
    (HERE/(stem+'-initial.json')).write_text(json.dumps(p,indent=2)+'\n')
    started=time.monotonic()
    try:
        result=radial.simulate('A',3e8,1e8,300.,zones=zones,nuclear=nuclear,reduced_a=True,
            radiation=True,initial_profile=p,final_time_units=final_units,
            save_trace=HERE/(stem+'.npz'),rtol=2e-6,atol=1e-11,art_linear=art_linear,art_quad=art_quad)
    except Exception as exc:
        failure={'record_kind':'numerical failure; no physical outcome established',
                 'exception_type':type(exc).__name__,'message':str(exc),
                 'zones':zones,'initial_temperature_K':temperature,
                 'contact_density_kg_m3':contact_density,'nuclear_enabled':nuclear,
                 'artificial_viscosity_linear':art_linear,'artificial_viscosity_quadratic':art_quad,
                 'wall_seconds':time.monotonic()-started}
        (HERE/(stem+'-failure.json')).write_text(json.dumps(failure,indent=2)+'\n')
        raise
    result['cold_driver_surrogate']=True
    result['cold_driver_limitations']=[
        'No discrete-grain contact, ionization, chemical binding, or Coulomb EOS.',
        'The inward velocity schedule is a deliverable kinematic requirement, not a constructed accelerator.',
        'Spherical symmetry and constant opacity; no conduction, screening, or multidimensional mixing.',
        'Initial temperature is a sensitivity input, not an externally delivered heating design.']
    (HERE/(stem+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('profile','output_mass_fractions')},indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--build-table',action='store_true')
    ap.add_argument('--zones',type=int,default=32);ap.add_argument('--temperature',type=float,default=1000.)
    ap.add_argument('--contact-density',type=float,default=2000.)
    ap.add_argument('--nuclear',action='store_true');ap.add_argument('--final-units',type=float,default=32.)
    ap.add_argument('--art-linear',type=float,default=.5);ap.add_argument('--art-quad',type=float,default=2.)
    args=ap.parse_args()
    if args.build_table:build_table()
    else:run(args.zones,args.temperature,args.contact_density,args.nuclear,args.final_units,args.art_linear,args.art_quad)
