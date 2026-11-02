"""Independent receiver and source-rating checks for compact opaque pulses.

This wraps the existing ballistic/conduction model with the new pulse outputs.
It does not turn a successful prepared-state hydro trajectory into a driver.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, math, sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'phase-two/plant-audit'))
import collector_screen as collector
from heat_transport import exchanger
SB=5.670374419e-8
STEMS=['scaled-hot-A-64-R30','scaled-hot-B1-64-R30','scaled-hot-C13-64-R30',
       'scaled-hot-He-64-R3','scaled-hot-B2-64-R30-norad','scaled-hot-B2-64-R30-tight']
RADIUS=1e9
ETA_DRIVER=.6

def rating(stem,pulse_root,d):
    p=json.loads((pulse_root/(stem+'-initial.json')).read_text())
    m=d['mass_kg']*np.asarray(p['mass_fraction_of_total']);r=np.asarray(p['outer_radius_m'])
    v=-np.asarray(p['velocity_m_s']);ri=np.r_[0.,r[:-1]];vi=np.r_[0.,v[:-1]]
    rc=.75*(r**4-ri**4)/(r**3-ri**3);vc=vi+(v-vi)*(rc-ri)/(r-ri)
    launch=-(RADIUS-rc)/vc
    edge=np.r_[launch[0]-.5*(launch[1]-launch[0]),.5*(launch[1:]+launch[:-1]),
               launch[-1]+.5*(launch[-1]-launch[-2])]
    node_m=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
    ke=.5*np.dot(node_m,v*v);cell=.5*m*vc**2;cell*=ke/cell.sum()
    assert np.all(np.diff(edge)>0)
    peak=float(np.max(cell/np.diff(edge)))/ETA_DRIVER
    return dict(driver_efficiency=ETA_DRIVER,initial_electric_energy_J=d['initial_energy_J']/ETA_DRIVER,
        installed_storage_at_100kJ_kg=d['initial_energy_J']/ETA_DRIVER/1e5,
        peak_kinetic_launch_electric_W=peak,peak_kinetic_launch_W_per_collector_m2=peak/(4*math.pi*RADIUS**2),
        launch_range_s=[float(launch.min()),float(launch.max())],
        initial_kinetic_reconstruction_relative_error=float(ke/d['initial_kinetic_energy_J']-1),
        source_rating_scope='Ballistic initial-cell source timing; kinetic electrical peak excludes separately timed internal-energy preparation and switching/field hardware.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pulse-root',type=Path,default=HERE/'small-stages');args=ap.parse_args()
    pulse_root=args.pulse_root;collector.ROOT=pulse_root;out=[];hashes={}
    cache_path=HERE/'opaque-receiver-prompt-traces.json.gz'
    cache=json.loads(gzip.decompress(cache_path.read_bytes())) if cache_path.exists() else {'traces':{}}
    raw_prompt=collector.prompt
    def saved_prompt(d,stem,radius):
        if stem not in STEMS:return raw_prompt(d,stem,radius)
        if stem not in cache['traces']:
            t,q,dt=raw_prompt(d,stem,RADIUS)
            npz=pulse_root/(stem+'.npz')
            cache['traces'][stem]={'time_s':t.tolist(),'flux_at_R1e9_W_m2':q.tolist(),'dt_s':dt,
                'npz_sha256':hashlib.sha256(npz.read_bytes()).hexdigest(),
                'record_sha256':hashlib.sha256((pulse_root/(stem+'.json')).read_bytes()).hexdigest()}
        z=cache['traces'][stem]
        return np.asarray(z['time_s']),np.asarray(z['flux_at_R1e9_W_m2'])*(RADIUS/radius)**2,z['dt_s']
    collector.prompt=saved_prompt
    for stem in STEMS:
        path=pulse_root/(stem+'.json')
        if not path.exists():continue
        d=json.loads(path.read_text());r=dict(stem=stem,source_success=d['success'],
            termination=d.get('effective_termination_reason',d.get('termination_reason')),
            mass_kg=d['mass_kg'],initial_energy_J=d['initial_energy_J'],
            final_kinetic_J=d['final_kinetic_J'],final_internal_J=d['final_internal_J'],
            source_energy_error_fraction_initial=d['energy_residual_fraction_initial'])
        for ext in ['.json','.npz','-initial.json']:
            f=pulse_root/(stem+ext)
            if f.exists():hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
        r['source_rating']=rating(stem,pulse_root,d)
        r['uniform_bulk_temperature_rise_K_at_Sigma10000_cp1000']=(d['final_kinetic_J']+d['final_internal_J'])/(4*math.pi*RADIUS**2*1e7)
        try:r['receiver']=collector.run(stem,RADIUS)
        except AssertionError as ex:r['receiver_failure']=str(ex)
        if not d['success']:r['qualification']='Partial failed numerical trajectory; source and energy scales only. Receiver history unqualified.'
        else:r['qualification']='Completed prescribed-flow source trajectory; collector remains ballistic and constant-property only.'
        out.append(r)
    # Direct check of the old large factory catches a previously incorrect R scaling.
    collector.ROOT=ROOT
    old=collector.run('compression-He-64-r30-burning',RADIUS)
    heat=exchanger(1.,eta=.4,relative_Carnot=.8,Th=1000.,Tr=330.)
    result=dict(scope=__doc__,radius_m=RADIUS,area_m2=4*math.pi*RADIUS**2,
        absorber_areal_mass_kg_m2=1e4,absorber_mass_kg=4*math.pi*RADIUS**2*1e4,
        illustrative_carbon_properties=collector.PROPS['illustrative_carbon'],cases=out,
        old_30m_factory_receiver=old,
        finite_engine_per_W_hot_input=heat,
        radiator_330K_emissivity_point9_flux_W_m2=.9*SB*330**4,
        radiator_maximum_W_on_same_sphere=.9*SB*330**4*4*math.pi*RADIUS**2,
        limitations=['No nonlinear material properties, erosion, implantation or catalyst retention calculation',
          'Summed independent thermal maxima are an upper screen only within constant-property conduction assumptions',
          'B2 partial profiles require continued hydrodynamics and can have colliding shells',
          'Surface radiation is counted in receiver loading even where cycle ledger subtracts escaped prompt photons; the final plant heat ledger must choose one collection convention',
          'Peak source powers are ratings, not manufactured components or proof of pulse preparation'],
        input_sha256=hashes,
        calculation_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),ROOT/'collector_screen.py',ROOT/'phase-two/plant-audit/heat_transport.py']})
    cache_path.write_bytes(gzip.compress(json.dumps(cache,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    result['prompt_trace_cache_sha256']=hashlib.sha256(cache_path.read_bytes()).hexdigest()
    (HERE/'opaque-receiver-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in out:print(r['stem'],r['source_success'],'peak_source',r['source_rating']['peak_kinetic_launch_electric_W'],
        'stored_J',r['source_rating']['initial_electric_energy_J'],
        'surface_sum_K',r.get('receiver',{}).get('thermal_results',{}).get('illustrative_carbon',{}).get('sum_of_separate_peak_rises_K_upper_screen'),r.get('receiver_failure',''))
    print('330K sphere radiation W',result['radiator_maximum_W_on_same_sphere'])
    print('finite exchanger per W',heat)
if __name__=='__main__':main()
