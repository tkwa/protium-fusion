"""Bounded periodic accumulation in a finite cooled carbon slab.

This is an explicit additional scheduling/cooling requirement: at least one
specified guard interval between any two large A/B2/C13 pulses, and periodic
B1 pulses separated by at least their complete cavity slot. The front is heated, the back held
at1000K. Constant properties and linear superposition remain assumptions.
"""
from pathlib import Path
import argparse,hashlib,json,math,sys
import numpy as np
from scipy.special import zeta
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import collector_screen as physical
RADIUS=1.1e9;K=100.;RHO=2000.;CP=1000.;SIGMA=5000.;L=SIGMA/RHO
ALPHA=K/(RHO*CP);AREA=4*math.pi*RADIUS**2;YEAR=365.25*86400
GUARD=86400.;EARLY=20000.;SINK=1000.

def history_impulse_kernel(period,last_heat_time):
    """Sum h(nP-D), n>=1, exactly per finite-slab thermal eigenmode.

    h(t)=2alpha/(kL)*sum exp[-alpha*((n+.5)pi/L)^2*t].
    Old pulse heat is conservatively concentrated at the latest time D.
    """
    assert 0<last_heat_time<period
    eigen=ALPHA*((np.arange(1000)+.5)*math.pi/L)**2
    return float(2*ALPHA/(K*L)*sum(np.exp(-eigen*(period-last_heat_time))/(-np.expm1(-eigen*period))))

def tail_flux_bound(cat,record,screen,period,cut):
    # After cut, v<R/cut. Bound dm/dv by its largest value in any
    # intersected remaining velocity bin, giving q(t)<=C*t^-4.
    v=np.r_[0.,record['profile']['velocity_m_s']]
    dm=record['mass_kg']*np.asarray(record['profile']['mass_fraction_of_total'])
    slope=dm/np.diff(v)*screen['ballistic']['velocity_energy_normalization']
    remaining=slope[v[:-1]<=RADIUS/cut]
    coeff=.5*max(remaining)*RADIUS**3/AREA*(1+record['final_internal_J']/record['final_kinetic_J'])
    ballistic=coeff/period**4*float(zeta(4,cut/period))
    decay=0.;terms={}
    for parent,daughter in [('n13','c13'),('o14','n14'),('o15','n15')]:
        rate=physical.LIB.get_rate_by_name(parent+'(,)'+daughter).eval(1e8)
        e=screen['delayed_decay_energy_J'][parent]/AREA
        term=e*rate*math.exp(-rate*cut)/(-math.expm1(-rate*period))
        terms[parent]=term;decay+=term
    return {'ballistic_W_m2':ballistic,'decay_W_m2':decay,'decay_terms':terms,
            'temperature_upper_K':(ballistic+decay)*L/K}

def main():
    global RADIUS,AREA
    ap=argparse.ArgumentParser()
    ap.add_argument('--station-inputs',type=Path,default=HERE/'station-inputs.json')
    ap.add_argument('--collector-inputs',type=Path,default=HERE/'collector-inputs.json')
    ap.add_argument('--output',type=Path,default=HERE/'collector-periodic.json')
    ap.add_argument('--radius',type=float,default=RADIUS)
    args=ap.parse_args();RADIUS=args.radius;AREA=4*math.pi*RADIUS**2
    station=json.loads(args.station_inputs.read_text());cycle=station['cycle'];receiver=station['receiver']
    inputs=json.loads(args.collector_inputs.read_text())
    rate=.8/(receiver['occupied_s_per_He_at_R1e9']['mass999']*RADIUS/1e9)
    intervals={k:1/(rate*p['pulses_per_He']) for k,p in cycle['stage_pulse_data'].items()}
    big=['c12','c13','n14'];small='n15'
    big_rate=sum(1/intervals[k] for k in big)
    assert big_rate*GUARD<1
    fluence={k:s['received_energy_J']/AREA for k,s in receiver['stages'].items()}
    x=1e9/RADIUS
    isolated={k:s['material_peak_rise_K']*x**2.5+s['prompt_peak_rise_K']*x**2 for k,s in receiver['stages'].items()}
    big_history=max(fluence[k] for k in big)*history_impulse_kernel(GUARD,EARLY)
    # Stage-by-stage late envelopes, each allowed the entire large-shot rate,
    # are summed rather than taking their max: intentionally conservative.
    big_tail={k:tail_flux_bound(k,inputs['stages'][k]['record'],receiver['stages'][k],GUARD,EARLY) for k in big}
    small_minimum=receiver['stages'][small]['timing_at_R1e9']['mass999']['occupied_s']*RADIUS/1e9
    small_cut=min(EARLY,small_minimum/2)
    small_history=fluence[small]*history_impulse_kernel(small_minimum,small_cut)
    small_tail=tail_flux_bound(small,inputs['stages'][small]['record'],receiver['stages'][small],small_minimum,small_cut)
    residual=big_history+small_history+sum(z['temperature_upper_K'] for z in big_tail.values())+small_tail['temperature_upper_K']
    phase_independent=max(isolated[k] for k in big)+isolated[small]+residual
    mean_flux=sum(fluence[k]/intervals[k] for k in fluence)
    all_phase_loose=sum(isolated.values())
    # Mode summation for the steady unit-flux thermal resistance approachesL/k.
    eigen=ALPHA*((np.arange(100000)+.5)*math.pi/L)**2
    resistance=float(sum(2*ALPHA/(K*L)/eigen))
    assert abs(resistance/(L/K)-1)<1e-5
    out={'scope':__doc__,'radius_m':RADIUS,'sigma_kg_m2':SIGMA,'thickness_m':L,
       'thermal_diffusivity_m2_s':ALPHA,'slowest_thermal_decay_time_s':4*L**2/(ALPHA*math.pi**2),
       'back_face_sink_K':SINK,'He_rate_s':rate,'per_stage_intervals_s':intervals,
       'large_shots_per_year':big_rate*YEAR,'large_shot_guard_s':GUARD,
       'large_guard_capacity_used':big_rate*GUARD,'early_pulse_bound_cut_s':EARLY,
       'B1_minimum_spacing_bound_s':small_minimum,'B1_mean_interval_s':intervals[small],
       'received_heat_mean_W_m2':mean_flux,'mean_front_to_sink_offset_K':mean_flux*L/K,
       'mean_front_temperature_K':SINK+mean_flux*L/K,
       'isolated_stage_peak_bounds_K':isolated,
       'sum_all_stage_isolated_maxima_with_arbitrary_coincidence_K':all_phase_loose,
       'old_large_pulse_background_upper_K':big_history,'old_B1_background_upper_K':small_history,
       'large_late_tail_envelopes':big_tail,'B1_late_tail_envelope':small_tail,
       'accumulated_background_upper_K':residual,
       'phase_independent_peak_rise_bound_with_guard_K':phase_independent,
       'phase_independent_front_temperature_bound_K':SINK+phase_independent,
       'extra_cold_inventory_upper_at_one_day_gross_flow_kg':cycle['gross_fuel_kg_per_He']*rate*GUARD,
       'steady_resistance_mode_sum_relative_error':resistance/(L/K)-1,
       'bound_logic':['At most one large pulse contributes its full isolated maximum; its predecessors are at least86400s apart.',
         'Each old large pulse is bounded by the maximum full-pulse fluence deposited entirely at20000s, plus a separate late-tail bound.',
         'B1 may have its isolated maximum at the same instant as the large-pulse maximum, with all previous pulses at their minimum allowable spacing included; irregular cadence can only reduce this bound.',
         'Positive finite-slab impulse response decreases with age; the late ballistic t^-4 and beta-exponential flux envelopes are bounded by a constant flux timesL/k.',
         'The sink is held at1000K, so mean and residual heating are included once. No free reset of absorber temperature is assumed.'],
       'limitations':['Requires the explicit large-pulse guard schedule, cold stock routing and distributed backside heat extraction. Mean count alone is only a scheduling-capacity check.',
         'No nonlinear material, erosion, radiation cooling or multidimensional receiver calculation.',
         'All incoming pulse heat is placed in the slab, conservative relative to99% captured-energy assumption.',
         'Provisional64 pulse family; source refinement may change quantitative bounds.'],
       'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.station_inputs,args.collector_inputs]},
       'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['large_shots_per_year','large_guard_capacity_used','received_heat_mean_W_m2','mean_front_to_sink_offset_K','isolated_stage_peak_bounds_K','accumulated_background_upper_K','phase_independent_front_temperature_bound_K','extra_cold_inventory_upper_at_one_day_gross_flow_kg']},indent=2))
if __name__=='__main__':main()
