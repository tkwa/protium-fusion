"""Ballistic pulse arrival and a deliberately simple thermal-collector screen.

The radial hydro result supplies a monotone velocity distribution. Subsequent
free expansion is approximated as ballistic; its remaining internal energy is
distributed with kinetic energy as a separate sensitivity. No wall survival,
catalyst retention, magnetic capture, or multidimensional stability is claimed.
"""
import json
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
from nuclear_screen import MU,MEV,LIB
from cycle_ledger import nuc,MEAN_NU

ROOT=Path(__file__).parent
FILES=['compression-A-128-full-burning','compression-A-512-burning',
       'compression-B1-64-burning','compression-B2-64-burning',
       'compression-B1-128-burning','compression-B2-128-burning',
       'compression-C13-64-burning','compression-He-64-burning',
       'compression-He-64-r30-burning']
PROPS={'illustrative_carbon':{'k_W_m_K':100.,'rho_kg_m3':2000.,'cp_J_kg_K':1000.},
       'tungsten_paper_constants':{'k_W_m_K':130.,'rho_kg_m3':19300.,'cp_J_kg_K':140.}}

def temperature_rise(flux,dt,props):
    # Exact convolution weights for a piecewise-constant heat flux on a
    # semi-infinite solid with constant properties and no surface cooling.
    n=np.arange(len(flux),dtype=float)
    kernel=2*np.sqrt(dt)*(np.sqrt(n+1)-np.sqrt(n))
    eff=np.sqrt(np.pi*props['k_W_m_K']*props['rho_kg_m3']*props['cp_J_kg_K'])
    return fftconvolve(flux,kernel)[:len(flux)]/eff

def ballistic(d,radius):
    v=np.r_[0.,d['profile']['velocity_m_s']]
    dm=d['mass_kg']*np.array(d['profile']['mass_fraction_of_total'])
    assert np.all(np.diff(v)>0), 'Shell collision invalidates ballistic screen'
    slope=dm/np.diff(v)
    e=slope*(v[1:]**3-v[:-1]**3)/6
    correction=d['final_kinetic_J']/sum(e)
    slope*=correction;e*=correction
    # Keep the slow tail until less than 1e-6 of kinetic energy remains.
    e_cum=np.r_[0.,np.cumsum(e)]
    vlow=np.interp(e_cum[-1]*1e-6,e_cum,v)
    tmax=radius/vlow
    nt=60000;dt=tmax/nt;t=(np.arange(nt)+.5)*dt
    speed=radius/t
    ids=np.searchsorted(v,speed,side='right')-1
    valid=(ids>=0)&(ids<len(dm))
    flux=np.zeros(nt)
    flux[valid]=.5*speed[valid]**4/radius*slope[ids[valid]]/(4*np.pi*radius**2)
    # Boundary interpolation gives a small quadrature discrepancy; report it.
    deposited=flux.sum()*dt*4*np.pi*radius**2
    tq=[radius/np.interp(e_cum[-1]*f,e_cum,v) for f in [.99,.5,.01]]
    return t,flux,dt,{'arrival_1_50_99_percent_energy_s':tq,
       'energy_quadrature_relative_error':deposited/d['final_kinetic_J']-1,
       'continuous_velocity_energy_normalization':correction,
       'truncated_slow_tail_energy_fraction':1e-6}

def prompt(d,stem,radius):
    archive=ROOT/'compression-traces.json'
    compact=json.loads(archive.read_text())['traces'].get(stem) if archive.exists() else None
    if compact is not None:
        tt=np.array(compact['time_s']);cumulative=np.array(compact['cumulative_radiation_escaped_J'])
    else:
        tr=np.load(ROOT/(stem+'.npz'));tt=tr['time_s']
        nz=int(tr['zones']);nv=int(tr['variables_per_zone']);nn=nv-5
        z=tr['state'].reshape(nz,nv,len(tt))
        weights=np.array(d['profile']['mass_fraction_of_total'])
        cumulative=weights@z[:,nn+4,:]
        cumulative*=d['radiation_escaped_J']/cumulative[-1]
    assert np.min(np.diff(cumulative))>-1e-7*max(cumulative)
    grid=np.linspace(0,tt[-1],12001);dt=grid[1]-grid[0]
    flux=np.diff(np.interp(grid,tt,cumulative))/(dt*4*np.pi*radius**2)
    return grid[1:],flux,dt

def cold_decay(d,radius,t):
    area=4*np.pi*radius**2;flux=np.zeros_like(t);energies={}
    for parent,daughter in [('n13','c13'),('o14','n14'),('o15','n15')]:
        count=d['mass_kg']/MU*d['output_mass_fractions'].get(parent,0.)/nuc(parent).A
        heat=(nuc(parent).mass-nuc(daughter).mass-MEAN_NU[parent])*MEV*count
        rate=LIB.get_rate_by_name(parent+'(,)'+daughter).eval(1e8)
        flux+=heat*rate*np.exp(-rate*t)/area
        energies[parent]=heat
    return flux,energies

def run(stem,radius):
    d=json.loads((ROOT/(stem+'.json')).read_text())
    t,q,dt,details=ballistic(d,radius)
    qp,energies=cold_decay(d,radius,t)
    tp,pq,pdt=prompt(d,stem,radius)
    # The photons/positrons have an earlier light-flight delay, omitted from
    # this separate temperature calculation. Their spatial source is treated
    # as central; this is adequate only when v*t_beta << collector radius.
    result={'source':stem+'.json','collector_radius_m':radius,
      'collector_area_m2':4*np.pi*radius**2,**details,
      'kinetic_fluence_J_m2':d['final_kinetic_J']/(4*np.pi*radius**2),
      'peak_ballistic_kinetic_flux_W_m2':float(max(q)),
      'peak_prompt_radiation_flux_W_m2':float(max(pq)),
      'prompt_radiation_fluence_J_m2':float(sum(pq)*pdt),
      'cold_decay_energy_J':energies,'thermal_results':{}}
    heat_factor=1+d['final_internal_J']/d['final_kinetic_J']
    for name,props in PROPS.items():
        main=temperature_rise(q*heat_factor,dt,props)
        cold=temperature_rise(qp,dt,props)
        direct=temperature_rise(pq,pdt,props)
        result['thermal_results'][name]={
          'kinetic_plus_residual_heat_peak_rise_K':float(max(main)),
          'cold_decay_all_surface_absorption_peak_rise_K':float(max(cold)),
          'prompt_all_surface_absorption_peak_rise_K':float(max(direct)),
          'sum_of_separate_peak_rises_K_upper_screen':float(max(main)+max(cold)+max(direct)),
          'peak_main_time_s':float(t[np.argmax(main)]),
          'main_equivalent_heat_flux_factor_MJ_m2_s_half':float(max(main)*np.sqrt(np.pi*props['k_W_m_K']*props['rho_kg_m3']*props['cp_J_kg_K'])/2e6)}
    return result

if __name__=='__main__':
    # A rectangular-pulse check verifies the convolution normalization.
    p=PROPS['illustrative_carbon'];dt=.01;q=2e5*np.ones(1000)
    got=temperature_rise(q,dt,p)[-1]
    expected=2*q[0]*np.sqrt(10/(np.pi*p['k_W_m_K']*p['rho_kg_m3']*p['cp_J_kg_K']))
    assert abs(got/expected-1)<1e-12
    cases=[run(f,r) for r in [1e10,2e10] for f in FILES]
    out={'scope':'Ballistic and constant-property semi-infinite thermal screen; no validated collector design.',
      'material_properties':PROPS,'rectangular_pulse_relative_error':got/expected-1,
      'limitations':['No continuing hydrodynamic acceleration or angular structure',
      'No nonlinear material response, thermal radiation, active cooling, erosion or implantation',
      'Gamma energy treated as surface heat for a conservative temperature screen; no photon transport',
      'Escaped prompt energy and timing inherit the grey diffusion model and its resolution limits',
      'Sum of separate peak temperatures is a conservative superposition screen, not a temperature history'],
      'cases':cases}
    (ROOT/'collector-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for d in cases:
        print(d['source'],d['collector_radius_m'],
          'MW/m2',d['peak_ballistic_kinetic_flux_W_m2']/1e6,
          'carbon rise K',d['thermal_results']['illustrative_carbon']['sum_of_separate_peak_rises_K_upper_screen'])
