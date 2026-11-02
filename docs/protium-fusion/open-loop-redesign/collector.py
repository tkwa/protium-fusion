"""Receiver, source and finite-temperature screens for corrected open-loop pulses.

Reproduce from saved small inputs: python collector.py
Refresh raw trace reduction: python collector.py --refresh-inputs
No hydro chain is run. The cached records, initial profiles and radiation
integrals preserve exactly the information this calculation consumes.
"""
from pathlib import Path
import argparse, hashlib, json, math, sys
import numpy as np
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import collector_screen as old
from fast_eos import FastEOS
RREF=1e9
AREA=4*math.pi*RREF**2
SB=5.670374419e-8
PROPS=old.PROPS['illustrative_carbon']

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def cache_inputs(process):
    result={'scope':'Reduced corrected pulse inputs; no hydro rerun required. Radiation integral explicitly selects radiated_energy/uunit, not NV-5.', 'stages':{}}
    for cat,s in process['selected']['stage_pulse_data'].items():
        path=Path(s['source']);d=json.loads(path.read_text());stage=d['stage']
        initial=Path(d.get('initial_profile_path',ROOT/'inputs/pulse-profiles'/f'scaled-hot-{stage}-64-R30-initial.json'))
        trpath=path.with_suffix('.npz')
        if not trpath.exists():
            directory='protium-b2-completion' if stage=='B2' else 'protium-b2-regression'
            trpath=Path('/var/tmp/tkwa-chatgpt')/directory/trpath.name
        tr=np.load(trpath);t=tr['time_s'];nz=int(tr['zones']);nv=int(tr['variables_per_zone'])
        z=tr['state'].reshape(nz,nv,len(t));index=d['state_variables'].index('radiated_energy/uunit')
        w=np.asarray(d['profile']['mass_fraction_of_total'])
        cumulative=w@z[:,index,:]
        cumulative*=d['radiation_escaped_J']/cumulative[-1]
        assert np.min(np.diff(cumulative))>=-1e-7*cumulative[-1]
        result['stages'][cat]={'record':d,'initial':json.loads(initial.read_text()),
           'prompt_time_s':t.tolist(),'prompt_cumulative_J':cumulative.tolist(),
           'source_sha256':{'record':sha(path),'initial':sha(initial),'raw_trace':sha(trpath)},
           'source_paths':{'record':str(path),'initial':str(initial),'raw_trace':str(trpath)}}
    (HERE/'collector-inputs.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result

def distribution(profile,mass):
    v=np.r_[0.,np.abs(profile['velocity_m_s'])]
    dm=mass*np.asarray(profile['mass_fraction_of_total'])
    assert np.all(np.diff(v)>0), 'Nonmonotone velocities invalidate this ballistic distribution'
    slope=dm/np.diff(v)
    energy=slope*(v[1:]**3-v[:-1]**3)/6
    return v,dm,slope,energy

def velocity_quantile(profile,mass,fraction,weight):
    v,dm,slope,energy=distribution(profile,mass)
    bins=dm if weight=='mass' else energy
    c=np.r_[0.,np.cumsum(bins)];target=fraction*c[-1]
    i=min(len(dm)-1,max(0,np.searchsorted(c,target,side='right')-1))
    if weight=='mass':return float(v[i]+(target-c[i])/slope[i])
    return float((v[i]**3+6*(target-c[i])/slope[i])**(1/3))

def representatives(profile):
    r=np.asarray(profile['outer_radius_m']);ri=np.r_[0.,r[:-1]]
    v=np.abs(profile['velocity_m_s']);vi=np.r_[0.,v[:-1]]
    rc=.75*(r**4-ri**4)/(r**3-ri**3)
    vc=vi+(v-vi)*(rc-ri)/(r-ri)
    return rc,vc

def source_rating(d,initial):
    m=d['mass_kg']*np.asarray(initial['mass_fraction_of_total'])
    rc,vc=representatives(initial);v=np.abs(initial['velocity_m_s'])
    launch=-(RREF-rc)/vc
    edge=np.r_[launch[0]-.5*(launch[1]-launch[0]),.5*(launch[1:]+launch[:-1]),launch[-1]+.5*(launch[-1]-launch[-2])]
    assert np.all(np.diff(edge)>0)
    nodem=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
    ke=.5*np.dot(nodem,v*v);cell=.5*m*vc**2;cell*=ke/cell.sum()
    # Exact initial per-cell internal energy from the same EOS and reconstructed
    # Lagrangian density, placed in the same preparation bins as kinetic energy.
    ye,yi={'B2':(8/14,3/14),'C13':(7/14,2/14),'A':(8/15,2/15),'B1':(8/16,2/16)}[d['stage']]
    bounds=np.r_[0.,np.asarray(initial['outer_radius_m'])]
    density=m/(4*math.pi/3*np.diff(bounds**3))/1000
    internal=m*FastEOS().state(density,np.asarray(initial['T_K']),ye,yi).e*1e-4
    internal_error=sum(internal)/d['initial_internal_energy_J']-1
    assert abs(internal_error)<1e-10,internal_error
    total_peak=float(np.max((cell+internal)/np.diff(edge)))/.6
    return {'peak_kinetic_electric_W':float(np.max(cell/np.diff(edge)))/.6,
       'peak_total_prepared_electric_W':total_peak,
       'initial_internal_reconstruction_relative_error':float(internal_error),
       'electric_storage_J':d['initial_energy_J']/.6,
       'initial_kinetic_reconstruction_relative_error':float(ke/d['initial_kinetic_energy_J']-1),
       'scope':'Total prepared kinetic plus exact per-cell internal energy in launch-time bins at 60% efficiency; switching/focusing and realizing the thermodynamic path remain engineering requirements.'}

def ballistic(d):
    v,dm,slope,e=distribution(d['profile'],d['mass_kg'])
    correction=d['final_kinetic_J']/sum(e);slope*=correction
    vmin=velocity_quantile(d['profile'],d['mass_kg'],1e-6,'energy')
    nt=240000;dt=RREF/vmin/nt;t=(np.arange(nt)+.5)*dt
    speed=RREF/t;ii=np.searchsorted(v,speed,side='right')-1
    good=(ii>=0)&(ii<len(dm));q=np.zeros(nt)
    q[good]=.5*speed[good]**4/RREF*slope[ii[good]]/AREA
    residual=sum(q)*dt*AREA/d['final_kinetic_J']-1
    assert abs(residual)<2e-3, residual
    rise=old.temperature_rise(q*(1+d['final_internal_J']/d['final_kinetic_J']),dt,PROPS)
    return {'peak_rise_K':float(max(rise)), 'peak_flux_W_m2':float(max(q)),
       'peak_time_s':float(t[np.argmax(rise)]),'energy_quadrature_relative_error':float(residual),
       'velocity_energy_normalization':float(correction),
       'arrival_1_50_99_percent_energy_s':[RREF/velocity_quantile(d['profile'],d['mass_kg'],f,'energy') for f in [.99,.5,.01]]}

def stage_screen(item):
    d=item['record'];p=item['initial'];main=ballistic(d)
    time=np.asarray(item['prompt_time_s']);cum=np.asarray(item['prompt_cumulative_J'])
    grid=np.linspace(0,time[-1],12001);dt=grid[1]-grid[0]
    q=np.diff(np.interp(grid,time,cum))/(dt*AREA)
    prompt=float(max(old.temperature_rise(q,dt,PROPS)))
    # Fixed beta lifetimes do not scale with receiver radius. This isolated
    # uniform-area screen assigns all decay heat to the receiving surface.
    cold_dt=.1;cold_t=(np.arange(200000)+.5)*cold_dt
    cold_q,energies=old.cold_decay(d,RREF,cold_t)
    cold=float(max(old.temperature_rise(cold_q,cold_dt,PROPS)))
    timing={}
    for name,f,weight in [('mass999',.001,'mass'),('energy99',.01,'energy')]:
        incoming=RREF/velocity_quantile(p,d['mass_kg'],f,weight)
        outgoing=RREF/velocity_quantile(d['profile'],d['mass_kg'],f,weight)
        timing[name]={'incoming_s':incoming,'outgoing_s':outgoing,'occupied_s':incoming+outgoing}
    _,vi=representatives(p);_,vo=representatives(d['profile'])
    timing['full_resolved']={'incoming_s':float(RREF/min(vi)),'outgoing_s':float(RREF/min(vo)),
        'occupied_s':float(RREF/min(vi)+RREF/min(vo))}
    received=d['final_kinetic_J']+d['final_internal_J']+d['radiation_escaped_J']+sum(energies.values())
    return {'stage':d['stage'],'source_sha256':item['source_sha256'],
        'source_success':d['success'],'source_termination':d['termination_reason'],'source_zones':d['zones'],
        'material_peak_rise_K':main['peak_rise_K'],
        'prompt_peak_rise_K':prompt+cold,
        'prompt_only_peak_rise_K':prompt,'delayed_decay_peak_rise_K':cold,
        'sum_separate_maxima_K':main['peak_rise_K']+prompt+cold,
        'received_energy_J':received, 'prompt_radiation_J':d['radiation_escaped_J'],
        'delayed_decay_energy_J':energies, 'ballistic':main, 'timing_at_R1e9':timing,
        'source_rating':source_rating(d,p),
        'central_initial_shell_mass_fraction':p['mass_fraction_of_total'][0],
        'central_outgoing_shell_mass_fraction':d['profile']['mass_fraction_of_total'][0],
        'maximum_surface_radius_m':d['profile']['outer_radius_m'][-1]}

def engine(th,tr):
    eta=.4;phi=.8;ratio=1-eta/phi
    if ratio<=tr/th:return {'hot_K':th,'radiator_K':tr,'feasible':False}
    b=math.sqrt((1-eta)*ratio);tx=(tr+b*th)/(ratio+b);tc=ratio*tx
    return {'hot_K':th,'radiator_K':tr,'feasible':True,'engine_hot_K':tx,'engine_cold_K':tc,
        'total_conductance_per_W_hot_K_inverse':1/(th-tx)+(1-eta)/(tc-tr)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--refresh-inputs',action='store_true');args=ap.parse_args()
    process=json.loads((HERE/'process.json').read_text());cycle=process['selected']
    inputs=cache_inputs(process) if args.refresh_inputs else json.loads((HERE/'collector-inputs.json').read_text())
    assert all(cycle['stage_pulse_data'][k]['sha256']==v['source_sha256']['record'] for k,v in inputs['stages'].items())
    stages={k:stage_screen(v) for k,v in inputs['stages'].items()}
    occupied={mode:sum(stages[k]['timing_at_R1e9'][mode]['occupied_s']*cycle['stage_pulse_data'][k]['pulses_per_He'] for k in stages) for mode in ['mass999','energy99','full_resolved']}
    # k=100 is an illustrative bulk carbon property, not qualified irradiated material.
    mat=[]
    for sigma in [100.,1000.,5000.,10000.]:
        thickness=sigma/2000
        mat.append({'sigma_kg_m2':sigma,'thickness_m':thickness,
          'max_pulse_s_for_L_gt_4sqrt_alpha_t_at_k100':(thickness/4)**2/(100/2000/1000),
          'max_pulse_s_for_L_gt_4sqrt_alpha_t_at_k200':(thickness/4)**2/(200/2000/1000),
          'uncollided_10MeV_fraction':math.exp(-.001959*sigma),
          'uncollided_20MeV_fraction':math.exp(-.001575*sigma)})
    scans=[]
    for rise in [500.,1000.,1500.,2000.]:
        radii={k:brentq(lambda r:s['material_peak_rise_K']*(RREF/r)**2.5+s['prompt_peak_rise_K']*(RREF/r)**2-rise,1e8,1e10) for k,s in stages.items()}
        scans.append({'allowed_rise_K':rise,'minimum_R_m_in_surface_deposition_screen':max(radii.values()),'per_stage_R_m':radii})
    out={'scope':'Corrected 64-zone pulses; ballistic free flight and illustrative carbon conduction only. A conditional engineering envelope, not a qualified receiver or source.',
       'resolution_status':'Provisional 64-zone baseline; B2 nested128 increases nuclear yield by about5-6%, with256 convergence pending. No mesh-converged receiver claim.',
       'radius_reference_m':RREF,'material_properties':PROPS,
       'stages':stages,'minimum_screened_absorber_sigma':5000.,
       'peak_source_electric_W_at_R1e9':max(s['source_rating']['peak_total_prepared_electric_W'] for s in stages.values()),
       'peak_kinetic_source_electric_W_at_R1e9':max(s['source_rating']['peak_kinetic_electric_W'] for s in stages.values()),
       'occupied_s_per_He_at_R1e9':occupied,
       'He_rate_at_80pct_occupancy_per_s':{k:.8/v for k,v in occupied.items()},
       'timing_scope':{'mass999':'Incoming plus outgoing 99.9% mass velocity quantile, leaving 0.1% slow mass on each leg for explicitly protected staging/cleanup; it cannot be discarded within the fuel-return budget.',
          'energy99':'Incoming plus outgoing 99% kinetic-energy quantile; not a fuel-retention criterion.',
          'full_resolved':'All resolved mass-weighted shell representative speeds. Central-shell subcell distribution reaches zero; this finite diagnostic also requires an active slow-core pickup.',
          'literal_100_percent_ballistic_clearance':'No finite time for a continuous central profile whose velocity reaches zero.'},
       'radius_scaling':{'material_peak_rise':'R^-2.5 for ballistic kinetic plus residual internal heat',
         'prompt_peak_rise':'R^-2; this interface field includes both prompt photons and the separate fixed-lifetime beta-decay heat maximum',
         'occupied_time':'R','kinetic_source_peak':'R^-1','storage_energy':'independent of R','absorber_mass':'R^2'},
       'thermal_radius_screens':scans,'absorber_screens':mat,
       'radiator_capacity_per_outer_m2_at5AU_W':{str(t):.9*SB*t**4-1361/25/4 for t in [330,450,600]},
       'finite_engine_40pct_at80pct_Carnot':[engine(h,t) for h in [1000,1500,2000] for t in [330,450,600]],
       'sources':[{'url':'https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html','use':'Graphite uncollided-photon attenuation: mu/rho .001959 and .001575 m2/kg at 10 and 20 MeV. Not a full energy-containment calculation.'}],
       'limitations':['Hydro endpoint may retain pressure-driven acceleration; monotone terminal speeds alone do not prove ballistic convergence.',
         'Surface deposition of photons is a conservative thermal approximation; real gamma penetration may lower skin temperature while demanding a thick absorber.',
         'Fixed constant properties, no ablation, implantation, fatigue, radiative surface cooling or active cooling validation.',
         'Delayed heat is uniformly spread over receiver area; decay location, charged-particle collection and gamma transport remain design tasks.',
         'Peak and mean loads require distributed surface hardware; a central heat-engine/pumping layout is not granted.',
         'Warm receiver cavity radiation requires protected cold feed trajectories; cold pellet preparation is not thermally automatic.'],
       'input_sha256':{'process.json':sha(HERE/'process.json'),'collector-inputs.json':sha(HERE/'collector-inputs.json')},
       'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'collector_screen.py',ROOT/'cycle_ledger.py',ROOT/'nuclear_screen.py',ROOT/'fast_eos.py',ROOT/'electron-eos.npz']}}
    year=365.25*86400
    rate=out['He_rate_at_80pct_occupancy_per_s']['mass999']
    wall=[]
    for share in [.1,.2]:
        loss=share*.007*cycle['fresh_hydrogen_kg_per_He']*rate/AREA*year
        wall.append({'fraction_of_total_fresh_carbon_allowance':share,
          'reference_H_flow_kg_s':cycle['fresh_hydrogen_kg_per_He']*rate,
          'maximum_net_wall_carbon_loss_kg_m2_year':loss,
          'maximum_fraction_lost_from_retired_Sigma5000_at_1e_4_per_year':loss/.5,
          'carbon_thickness_equivalent_m_year':loss/2000})
    out['wall_carbon_budget_at_R1e9_mass999_cap']=wall
    out['wall_carbon_budget_scaling']='Net wall loss per area = reserved_fraction*0.007*fresh_H_kg_per_He*He_rate*year/(4*pi*R^2). Capital stock and gross internally recycled substrate are not fresh operating feed.'
    out['annual_surface_fuel_loading_at_R1e9_mass999_cap_kg_m2']={k:rate*cycle['stage_pulse_data'][k]['pulses_per_He']*year*inputs['stages'][k]['record']['mass_kg']/AREA for k in stages}
    if 'B2_eight_sound_time_sensitivity' in inputs:
        item=dict(inputs['stages']['c12']);item['record']=inputs['B2_eight_sound_time_sensitivity']['record']
        item['source_sha256']=inputs['B2_eight_sound_time_sensitivity']['source_sha256']
        out['B2_eight_sound_time_sensitivity']=stage_screen(item)
        out['B2_eight_sound_time_sensitivity']['scope']='Later hydrodynamic endpoint; prompt waveform retained from the better early-time baseline trace. Baseline process/receiver unchanged.'
    (HERE/'collector.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'stages':{k:{x:s[x] for x in ['material_peak_rise_K','prompt_only_peak_rise_K','delayed_decay_peak_rise_K','sum_separate_maxima_K']} for k,s in stages.items()},'occupied':occupied,'rates':out['He_rate_at_80pct_occupancy_per_s'],'radius':scans,'source_W':out['peak_source_electric_W_at_R1e9']},indent=2))

if __name__=='__main__':main()
