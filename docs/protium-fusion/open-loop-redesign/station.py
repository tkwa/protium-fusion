"""Coupled station mass, power, radiation and pulse-schedule requirements.

A finite family of conditional engineering allocations, not a global optimum.
Pure arithmetic: inputs summarize completed stage and receiver calculations.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

MEV=1.602176634e-13
G=6.67430e-11
SB=5.670374419e-8
YEAR=365.25*86400
RREF=1e9
MSUN=1.98847e30
AU=149597870700.


def evaluate(rate, radius, cycle, receiver, *, radiator_T=450., hot_T=1000.,
             absorber_sigma=10000., driver_peak_W_kg=2., service_sigma=23873.2414637843,
             reserve_fraction=.2, maximum_surface_T=2000., timing='mass999',
             source_energy_J_kg=1e5, driver_efficiency=.6, heat_efficiency=.4,
             energy_capture=.99, auxiliary_fraction=.1, replacement_fraction=1e-4,
             replacement_J_kg=1e8, feed_buffer_days=30., process_residence_s=86400.,
             total_supply_ratio=.007, surface_model="surface", soft_surface_fraction=0.,
             photon_depth_K=None):
    """Return one coupled allocation; source/receiver ratings are explicit.

    service_sigma retains the previous processing, structural and thermal
    allowance per receiver area. Peak driver, average conditioning and heat
    conversion are additional equipment. A fractional installed reserve is
    included, not omitted from the reported denominator.
    """
    A=4*math.pi*radius**2
    if heat_efficiency >= .8*(1-radiator_T/hot_T):
        return {'feasible':False,'reason':'finite_engine_temperature_margin'}
    # Finite hot/cold approaches, with an engine at 80% of internal Carnot.
    c=1-heat_efficiency/.8
    approach_ratio=math.sqrt((1-heat_efficiency)*c)
    engine_hot=(radiator_T+approach_ratio*hot_T)/(c+approach_ratio)
    engine_cold=c*engine_hot
    conductance_per_thermal_W=1/(hot_T-engine_hot)+(1-heat_efficiency)/(engine_cold-radiator_T)
    # MeV per net He converted to joules per net He.
    h=cycle['H_nonnu_MeV_per_He']*MEV
    u=cycle['U_prepared_MeV_per_He']*MEV
    eg=heat_efficiency*energy_capture*(h+u)
    ed=u/driver_efficiency
    en=eg-ed
    eq=energy_capture*(h+u)-u
    eexport=(1-auxiliary_fraction)*en
    if eexport<=0:
        return {'feasible':False,'reason':'nonpositive_process_margin'}
    flux=.9*SB*radiator_T**4
    solar=1361/(5**2)/4  # full absorption, sphere mean, 5 AU
    rep=replacement_fraction*replacement_J_kg/YEAR
    speed=1e4
    support_coeff=.001*G*speed/(8*math.pi*radius**2)
    stream_coeff=G/(2*radius*speed**2)
    # Conservative scalar load envelope: maximum Hill-frame solar tide
    # 3 GM_sun R/a^3 relative to shell self-gravity GM/(2 R^2).
    tidal_equivalent_mass=6*MSUN*(radius/(5*AU))**3
    tidal_support_linear=support_coeff*tidal_equivalent_mass
    peak_per_reference=receiver['peak_source_electric_W_at_R1e9']
    peak=peak_per_reference*RREF/radius*(.6/driver_efficiency)
    stored=2*max(x['initial_energy_J'] for x in cycle['stage_pulse_data'].values())/driver_efficiency
    circulation=cycle['gross_fuel_kg_per_He']*rate
    working=circulation*process_residence_s
    fresh=(1+total_supply_ratio)*cycle['fresh_hydrogen_kg_per_He']*rate
    feed_buffer=fresh*feed_buffer_days*86400
    # M = base + coeff1*M + coeff2*M^2. Positive small root is physical.
    base={
        'receiver':absorber_sigma*A,
        'insulation':100*A,
        'processing_structure_thermal_service':service_sigma*A,
        'peak_driver_hardware':peak/driver_peak_W_kg,
        'two_pulse_storage':stored/source_energy_J_kg,
        'average_driver_conditioning':.1*ed*rate,
        'gross_electric_conversion':.1*eg*rate,
        'working_inventory':working,
        'inflight_fuel_inventory':sum(v['gross_fuel_kg_per_He']*rate*receiver['stages'][k]['timing_at_R1e9'][timing]['occupied_s']*radius/RREF for k,v in cycle['stage_pulse_data'].items()),
        'fresh_feed_buffer':feed_buffer,
    }
    constant=sum(base.values())+20*((eq-eexport)*rate+solar*A)/flux
    linear=reserve_fraction+replacement_fraction*10+rep+20*(rep+tidal_support_linear)/flux+stream_coeff*tidal_equivalent_mass
    quadratic=stream_coeff+20*support_coeff/flux
    b=1-linear
    discriminant=b*b-4*quadratic*constant
    if b<=0 or discriminant<0:
        return {'feasible':False,'reason':'mass_support_feedback'}
    mass=2*constant/(b+math.sqrt(discriminant))
    replacement=rep*mass
    support=support_coeff*mass*(mass+tidal_equivalent_mass)
    power=eexport*rate-replacement-support
    waste=eq*rate-power
    radiator_area=(waste+solar*A)/flux
    base.update({
        'radiator_panels':20*radiator_area,
        'replacement_manufacturing':replacement,
        'ten_year_spares':replacement_fraction*10*mass,
        'support_streams':stream_coeff*mass*(mass+tidal_equivalent_mass),
        'installed_reserve':reserve_fraction*mass,
    })
    occupancy=rate*receiver['occupied_s_per_He_at_R1e9'][timing]*radius/RREF
    # The thermal screen supplies an upper sum of surface maxima, at RREF.
    # Prompt and material arrival contributions have different R scaling.
    if surface_model=='surface':
        prompt_rise={k:x['prompt_peak_rise_K'] for k,x in receiver['stages'].items()}
    elif surface_model=='finite_depth':
        assert photon_depth_K is not None
        prompt_rise={k:soft_surface_fraction*x['prompt_only_peak_rise_K']+(1-soft_surface_fraction)*photon_depth_K[k]+x['delayed_decay_peak_rise_K'] for k,x in receiver['stages'].items()}
    else:
        raise ValueError(surface_model)
    rise=max(x['material_peak_rise_K']*(RREF/radius)**2.5+prompt_rise[k]*(RREF/radius)**2 for k,x in receiver['stages'].items())
    bulk_rise=max(x['received_energy_J']/(A*absorber_sigma*1000)
                  for x in receiver['stages'].values())
    reasons=[]
    if power<=0:reasons.append('nonpositive_net_power')
    if radiator_area>A:reasons.append('radiator_area')
    if occupancy>.8*(1+1e-12):reasons.append('sequential_pulse_occupancy')
    if peak<ed*rate:reasons.append('source_peak_below_average')
    if cycle['fresh_C_H_mass_ratio']>total_supply_ratio+1e-12:reasons.append('carbon_supply_budget')
    if hot_T+rise>maximum_surface_T:reasons.append('surface_temperature')
    # The supplied receiver screen sets a thickness floor for the declared
    # thermal/attenuation approximation; actual transport remains separate.
    if absorber_sigma<receiver.get('minimum_screened_absorber_sigma',5000.):
        reasons.append('absorber_thickness_scope')
    electrical_residual=eg*rate-ed*rate-auxiliary_fraction*en*rate-replacement-support-power
    radiation_escape=(1-energy_capture)*(h+u)*rate
    global_residual=h*rate-power-waste-radiation_escape
    return {
        'feasible':not reasons,'reasons':reasons,
        'helium_per_s':rate,'radius_m':radius,'mass_kg':mass,
        'net_W':power,'net_W_kg':power/mass,
        'gross_W':eg*rate,'driver_W':ed*rate,'auxiliary_W':auxiliary_fraction*en*rate,
        'supply_boundary':'Externally provided feed; extraction/delivery apparatus and costs are outside this station boundary.','replacement_W':replacement,'support_W':support,
        'waste_W':waste,'absorbed_solar_W':solar*A,'radiator_area_fraction':radiator_area/A,
        'surface_model':surface_model,'soft_surface_fraction':soft_surface_fraction,
        'engine_internal_hot_K':engine_hot,'engine_internal_cold_K':engine_cold,
        'heat_exchanger_conductance_W_K':conductance_per_thermal_W*eg*rate/heat_efficiency,
        'radiator_K':radiator_T,'receiver_base_K':hot_T,'peak_surface_K':hot_T+rise,
        'single_pulse_bulk_rise_K':bulk_rise,
        'absorber_sigma_kg_m2':absorber_sigma,'service_sigma_kg_m2':service_sigma,
        'driver_peak_W_kg':driver_peak_W_kg,'source_peak_W':peak,'stored_electric_J':stored,
        'timing_basis':timing,'sequential_occupancy':occupancy,
        'feed_C_kg_s':cycle['fresh_carbon_kg_per_He']*rate,
        'feed_H_kg_s':cycle['fresh_hydrogen_kg_per_He']*rate,
        'fresh_C_H_mass_ratio':cycle['fresh_C_H_mass_ratio'],
        'total_C_H_supply_budget':total_supply_ratio,
        'net_material_C_allowance_kg_s':(total_supply_ratio-cycle['fresh_C_H_mass_ratio'])*cycle['fresh_hydrogen_kg_per_He']*rate,
        'net_material_C_allowance_kg_m2_year':(total_supply_ratio-cycle['fresh_C_H_mass_ratio'])*cycle['fresh_hydrogen_kg_per_He']*rate/A*YEAR,
        'receiver_refurbishment_C_return_requirement':1-(total_supply_ratio-cycle['fresh_C_H_mass_ratio'])*cycle['fresh_hydrogen_kg_per_He']*rate/(replacement_fraction*absorber_sigma*A/YEAR),
        'all_stage_flow_kg_s':circulation,'cold_process_buffer_s':process_residence_s,
        'pulse_intervals_s':{k:1/(rate*v['pulses_per_He']) for k,v in cycle['stage_pulse_data'].items()},
        'mass_breakdown_kg':base,'mass_relative_residual':sum(base.values())/mass-1,
        'electrical_relative_residual':electrical_residual/max(1.,eg*rate),
        'global_energy_relative_residual':global_residual/max(1.,h*rate),
        'tidal_equivalent_mass_kg':tidal_equivalent_mass,
        'maximum_tidal_to_self_load_ratio':tidal_equivalent_mass/mass,
        'replacement_fraction_per_year':replacement_fraction,
        'replacement_energy_J_kg':replacement_J_kg,
        'fabrication_repayment_years_at100MJkg':1e8/(power/mass)/YEAR if power>0 else None,
        'assumed_capital_per_W_not_economic_cost':mass/power if power>0 else None,
    }


def scan(cycle,receiver):
    rows=[]
    thermal=[(330.,1000.),(450.,1000.),(450.,1500.),(600.,1500.)]
    radii=[10**(8.55+i*(9.65-8.55)/60) for i in range(61)]
    for timing in ['mass999','full_resolved']:
      for service in [23873.2414637843,5000.,1000.]:
        for rating in [2.,10.,100.]:
          for sigma in [5000.,10000.]:
            best=None;count=0
            for tr,th in thermal:
                for radius in radii:
                    fmax=.8/(receiver['occupied_s_per_He_at_R1e9'][timing]*radius/RREF)
                    for j in range(20):
                        f=fmax*10**(-2+2*j/19)
                        d=evaluate(f,radius,cycle,receiver,radiator_T=tr,hot_T=th,
                                   absorber_sigma=sigma,driver_peak_W_kg=rating,
                                   service_sigma=service,timing=timing)
                        if d.get('feasible'):
                            count+=1
                            if best is None or d['net_W_kg']>best['net_W_kg']:best=d
            # Locally refine radius and throughput, while preserving discrete
            # engineering choices. This is a finite requirements family.
            if best:
              for width in [.06,.006,.0006]:
                center=best['radius_m']
                for i in range(61):
                  radius=center*10**(width*(i/30-1))
                  fmax=.8/(receiver['occupied_s_per_He_at_R1e9'][timing]*radius/RREF)
                  for factor in [1.,.98,.95,.9,.8,.6,.4,.2]:
                    d=evaluate(fmax*factor,radius,cycle,receiver,radiator_T=best['radiator_K'],hot_T=best['receiver_base_K'],absorber_sigma=sigma,driver_peak_W_kg=rating,service_sigma=service,timing=timing)
                    if d.get('feasible') and d['net_W_kg']>best['net_W_kg']:best=d
            rows.append({'peak_driver_rating_W_kg':rating,'absorber_sigma_kg_m2':sigma,
                         'service_sigma_kg_m2':service,'timing_basis':timing,
                         'feasible_tested_points':count,'best':best})
    return rows


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('station.json'))
    args=ap.parse_args();inputs=json.loads(args.inputs.read_text())
    rows=scan(inputs['cycle'],inputs['receiver'])
    good=[x['best'] for x in rows if x['best']]
    for d in good:
        assert abs(d['mass_relative_residual'])<1e-10
        assert abs(d['electrical_relative_residual'])<1e-10
        assert abs(d['global_energy_relative_residual'])<1e-10
    out={'scope':'Finite conditional station family; hardware ratings, service mass and qualified pulse-tail handling are required engineering achievements. Supplied-feed station boundary.',
         'inputs_sha256':hashlib.sha256(args.inputs.read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'family':rows,'best_within_tested_family':max(good,key=lambda x:x['net_W_kg']) if good else None,
         'checks':{'mass_power_global_energy':True,'feasible_family_count':len(good)}}
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{k:d[k] for k in ['peak_driver_rating_W_kg','absorber_sigma_kg_m2','service_sigma_kg_m2','timing_basis']}|{'best_W_kg':d['best']['net_W_kg'] if d['best'] else None} for d in rows],indent=2))
