# Catcher follow-up copy of station_generalized.py; archived source remains unchanged.
# Adds explicitly returned auxiliary dissipation to the hot loop, with cold/electrical closure.
# Generalized copy of frozen cold-sink-redesign/station.py.
# Original SHA256: 8437bd3df1c66b53940639d2ae98d71b8be9c913bee6c67b3bacb77b841ebea5
# Only capture and driver efficiency are exposed as parameters; their defaults
# preserve the original algebra. The original file remains unchanged.
"""Cold-sink thermal station accounting, with explicit finite exchangers.

Pulse rate, radius and source equipment follow the supplied frozen operating
point. Factory hardware is proportional to repaired throughput, independently
of its energy cost. Direct conversion is a separate marginal ledger below.
"""
import math

MEV=1.602176634e-13
G=6.67430e-11
SB=5.670374419e-8
YEAR=365.25*86400


def finite_engine(eta, phi, hot, radiator):
    """Minimum-UA temperature allocation at a chosen fraction phi of Carnot.

    eta is gross heat-to-electric efficiency. phi is a selected engineering
    quality, not a thermodynamic constant. The two exchanger conductances are
    added; heat-engine cold rejection is (1-eta) times its hot input.
    """
    if not (0<eta<1 and 0<phi<=1 and 0<radiator<hot):
        return {'valid':False,'reason':'invalid_engine_parameters'}
    limit=phi*(1-radiator/hot)
    if eta>=limit-1e-12*max(1.,limit):
        return {'valid':False,'reason':'no_finite_temperature_margin','external_quality_limit':limit}
    ratio=1-eta/phi
    a=math.sqrt((1-eta)*ratio)
    th=(radiator+a*hot)/(ratio+a)
    tc=ratio*th
    dh=hot-th;dc=tc-radiator
    if min(dh,dc)<=1e-6:
        return {'valid':False,'reason':'unresolved_exchanger_temperature_margin'}
    return {'valid':True,'internal_hot_K':th,'internal_cold_K':tc,
            'hot_approach_K':dh,'cold_approach_K':dc,
            'UA_per_hot_W_per_K':1/dh+(1-eta)/dc,
            'external_carnot':1-radiator/hot,'internal_carnot':1-tc/th,
            'fraction_of_internal_carnot':phi,
            'fraction_of_external_carnot':eta/(1-radiator/hot)}


def evaluate(baseline,cycle,*,environment_absorbed_W_m2,energy_capture_fraction=.99,driver_efficiency=.6,heat_efficiency=.6,carnot_quality=.8,
             radiator_K=200.,receiver_K=1000.,routine_MJ_kg=2.,
             direct_solar_absorbed_W_m2=0.,
             environment_basis='whole_shell',emissivity=.9,escape_view_factor=1.,
             radiator_available_fraction=1.,bypass_W_m2=0.,hot_eligible_fraction=1.,
             hx_W_K_kg=10.,thermal_converter_We_kg=10.,pump_fraction_hot=0.,
             hx_internal_UA_per_hot_W_K=.03,
             pump_We_kg=10.,radiator_kg_m2=20.,auxiliary_fraction=.1,
             include_incremental_hx=True, internal_hot_return_W_m2=0.):
    """One conditional requirement case, with all absorbed environment counted.

    environment_absorbed_W_m2 is a *net absorbed external* IR input, after its
    source view factors and absorptivity, not an incident blackbody irradiance.
    whole_shell multiplies that input by all 4piR²; radiator multiplies by Arad.
    escape_view_factor is the cold-opening cosine-weighted view. Its complement
    views identical gray radiators, giving an effective emission coefficient.
    It must not be applied to the already-absorbed environment input again.
    """
    params=dict(internal_hot_return_W_m2=internal_hot_return_W_m2,energy_capture_fraction=energy_capture_fraction,driver_efficiency=driver_efficiency,heat_efficiency=heat_efficiency,carnot_quality=carnot_quality,
        radiator_K=radiator_K,receiver_K=receiver_K,routine_MJ_kg=routine_MJ_kg,
        direct_solar_absorbed_W_m2=direct_solar_absorbed_W_m2,
        environment_absorbed_W_m2=environment_absorbed_W_m2,environment_basis=environment_basis,
        emissivity=emissivity,escape_view_factor=escape_view_factor,
        radiator_available_fraction=radiator_available_fraction,bypass_W_m2=bypass_W_m2,
        hot_eligible_fraction=hot_eligible_fraction,hx_W_K_kg=hx_W_K_kg,
        thermal_converter_We_kg=thermal_converter_We_kg,pump_fraction_hot=pump_fraction_hot,
        pump_We_kg=pump_We_kg,radiator_kg_m2=radiator_kg_m2,auxiliary_fraction=auxiliary_fraction,
        include_incremental_hx=include_incremental_hx,
        hx_internal_UA_per_hot_W_K=hx_internal_UA_per_hot_W_K)
    engine=finite_engine(heat_efficiency,carnot_quality,receiver_K,radiator_K)
    if not engine['valid']:
        return {'feasible':False,'reasons':[engine['reason']],'parameters':params,'engine':engine}
    assert 0<emissivity<=1 and 0<escape_view_factor<=1
    assert 0<radiator_available_fraction<=1 and 0<=hot_eligible_fraction<=1
    assert routine_MJ_kg>=0 and bypass_W_m2>=0 and environment_absorbed_W_m2>=0
    assert thermal_converter_We_kg>0 and pump_We_kg>0 and hx_W_K_kg>0
    assert hx_internal_UA_per_hot_W_K>=0
    assert 0<=pump_fraction_hot<1 and 0<=auxiliary_fraction<1
    assert environment_basis in ['whole_shell','radiator']
    R=baseline['radius_m'];A=4*math.pi*R*R;rate=baseline['helium_per_s']
    H=cycle['H_nonnu_MeV_per_He']*MEV*rate
    nuclear_total=cycle['nuclear_Q_MeV_per_He']*MEV*rate
    neutrino=cycle['neutrino_MeV_per_He']*MEV*rate
    U=cycle['U_prepared_MeV_per_He']*MEV*rate
    assert 0<energy_capture_fraction<=1 and 0<driver_efficiency<=1
    capture=energy_capture_fraction
    captured=capture*(H+U);escape=(1-capture)*(H+U)
    bypass=A*bypass_W_m2
    if bypass>=hot_eligible_fraction*captured:
        return {'feasible':False,'reasons':['bypass_exceeds_hot_supply'],'parameters':params}
    returned=A*internal_hot_return_W_m2
    assert returned>=0
    hot=hot_eligible_fraction*captured-bypass+returned
    cold_source=(1-hot_eligible_fraction)*captured
    gross=heat_efficiency*hot
    driver=U/driver_efficiency
    process_margin=gross-driver
    if process_margin<=0:
        return {'feasible':False,'reasons':['nonpositive_process_margin'],'parameters':params}
    aux=auxiliary_fraction*process_margin
    pump=pump_fraction_hot*hot
    before_massloads=process_margin-aux-pump
    external_UA_min=engine['UA_per_hot_W_per_K']*hot
    UA=(engine['UA_per_hot_W_per_K']+hx_internal_UA_per_hot_W_K)*hot
    hx_mass=UA/hx_W_K_kg if include_incremental_hx else 0.
    eta_emit=emissivity*escape_view_factor/(1-(1-emissivity)*(1-escape_view_factor))
    emitted_flux=eta_emit*SB*radiator_K**4
    full_env=direct_solar_absorbed_W_m2
    rad_env=0.
    if environment_basis=='whole_shell':full_env+=environment_absorbed_W_m2
    else:rad_env=environment_absorbed_W_m2
    net_capacity=emitted_flux-rad_env
    if net_capacity<=0:
        return {'feasible':False,'reasons':['environment_exceeds_radiator_capacity'],'parameters':params,
                'radiator_emission_W_m2':emitted_flux}
    omitted={'gross_electric_conversion','radiator_panels','replacement_manufacturing',
             'ten_year_spares','support_streams','installed_reserve'}
    masses={k:v for k,v in baseline['mass_breakdown_kg'].items() if k not in omitted}
    masses.update({'thermal_electric_converter':gross/thermal_converter_We_kg,
                   'incremental_heat_exchangers':hx_mass,'incremental_pump_hardware':pump/pump_We_kg})
    ref_flow=baseline['mass_kg']*baseline['replacement_fraction_per_year']/YEAR
    factory_mass_per_flow=baseline['mass_breakdown_kg']['replacement_manufacturing']/ref_flow
    repair_mass_fraction=factory_mass_per_flow*baseline['replacement_fraction_per_year']/YEAR
    rep=baseline['replacement_fraction_per_year']*routine_MJ_kg*1e6/YEAR
    spare_fraction=10*baseline['replacement_fraction_per_year'];reserve=.2
    speed=1e4
    ks=.001*G*speed/(8*math.pi*R*R)
    km=G/(2*R*speed*speed)
    Mt=baseline['tidal_equivalent_mass_kg']
    eq=captured-U
    # P=before_massloads-rep*M-ks*M(M+Mt), Arad=(eq-P+full_env*A)/capacity.
    constant=sum(masses.values())+radiator_kg_m2*(eq-before_massloads+full_env*A)/net_capacity
    linear=reserve+spare_fraction+repair_mass_fraction+km*Mt+radioratio(radiator_kg_m2,net_capacity)*(rep+ks*Mt)
    quadratic=km+radioratio(radiator_kg_m2,net_capacity)*ks
    b=1-linear;disc=b*b-4*quadratic*constant
    if b<=0 or disc<0:
        return {'feasible':False,'reasons':['mass_support_feedback'],'parameters':params}
    mass=2*constant/(b+math.sqrt(disc))
    support=ks*mass*(mass+Mt);repair=rep*mass
    power=before_massloads-support-repair
    waste=eq-power
    Arad=(waste+full_env*A)/net_capacity
    environment=full_env*A+rad_env*Arad
    masses.update({'radiator_panels':radiator_kg_m2*Arad,
                   'replacement_factory_capacity':repair_mass_fraction*mass,
                   'ten_year_spares':spare_fraction*mass,'support_streams':km*mass*(mass+Mt),
                   'installed_reserve':reserve*mass})
    engine_cold=(1-heat_efficiency)*hot
    driver_loss=driver-U
    heat_components={'engine_cold_rejection':engine_cold,'hot_to_cold_bypass':bypass,
        'captured_source_ineligible_for_hot_engine':cold_source,'driver_preparation_loss':driver_loss,
        'ordinary_auxiliaries_after_hot_return':aux-returned,'external_heat_loop_pumps':pump,
        'routine_refurbishment':repair,'support_dissipation':support}
    assert returned<=aux+1e-10*max(1.,aux)
    reasons=[]
    if power<=0:reasons.append('nonpositive_net_power')
    if Arad/A>radiator_available_fraction*(1+1e-12):reasons.append('radiator_area')
    residuals={'mass':sum(masses.values())/mass-1,
        'electrical':(gross-driver-aux-pump-repair-support-power)/max(1,gross),
        'global_nuclear':(H-escape-power-waste)/max(1,H),
        'global_including_neutrinos':(nuclear_total-neutrino-escape-power-waste)/max(1,nuclear_total),
        'cold_heat_components':(sum(heat_components.values())-waste)/max(1,waste),
        'radiator':(Arad*emitted_flux-environment-waste)/max(1,waste)}
    assert max(abs(v) for v in residuals.values())<1e-10
    return {'feasible':not reasons,'reasons':reasons,'parameters':params,'engine':engine,
        'mass_kg':mass,'net_W':power,'net_W_kg':power/mass,'radius_m':R,'area_m2':A,'helium_per_s':rate,
        'nuclear_non_neutrino_W':H,'nuclear_total_W':nuclear_total,'neutrino_W':neutrino,
        'net_fraction_of_nuclear_energy':power/nuclear_total,
        'net_fraction_of_non_neutrino_nuclear_energy':power/H,'prepared_pulse_W':U,'captured_pulse_W':captured,
        'escaping_source_W':escape,'engine_hot_W':hot,'gross_electric_W':gross,
        'driver_W':driver,'ordinary_auxiliary_W':aux,'pump_W':pump,'replacement_W':repair,'support_W':support,
        'cold_rejection_W':waste,'environment_absorbed_W':environment,'internal_hot_return_W':returned,
        'external_exchanger_minimum_UA_W_K':external_UA_min,
        'additional_internal_exchanger_UA_W_K':hx_internal_UA_per_hot_W_K*hot,
        'heat_exchanger_UA_W_K':UA,'heat_exchanger_UA_W_K_per_shell_m2':UA/A,
        'radiator_effective_emissivity_to_escape':eta_emit,'radiator_emission_W_m2':emitted_flux,
        'radiator_area_m2':Arad,'radiator_area_fraction':Arad/A,
        'radiator_net_capacity_W_m2':net_capacity,
        'max_environment_W_m2_at_available_area':(radiator_available_fraction*A*emitted_flux-waste)/(A if environment_basis=='whole_shell' else radiator_available_fraction*A),
        'hot_engine_input_W_m2':hot/A,'net_export_W_m2':power/A,'cold_rejection_W_m2':waste/A,
        'factory_mass_per_repaired_kg_s':factory_mass_per_flow,
        'routine_repaired_kg_s':mass*baseline['replacement_fraction_per_year']/YEAR,
        'mass_breakdown_kg':masses,'cold_heat_components_W':heat_components,'closure_relative_residuals':residuals,
        'scope':'Conditional stationary heat/power/mass requirements. Hot eligibility, insulation, converter/HX/pump ratings and cold-sky view require hardware qualification. Original generic services retained; named HX and pump hardware are incremental allocations.'}


def radioratio(mass_per_area,net_flux):
    return mass_per_area/net_flux


def direct_marginal_bound(thermal_case,kinetic_captured_W,recovered_fraction,
                          extra_direct_aux_W=0.,cold_dump_fraction=0.):
    """Energy-only marginal ledger; no peak equipment or storage assumed free.

    Recovered fraction converts eligible captured kinetic energy to bus work.
    The unconverted remainder is hot unless cold_dump_fraction sends it directly
    to the cold sink. This bound presumes that same energy was hot-eligible in
    the reference thermal case; otherwise reject/use a new routing model.
    """
    p=thermal_case['parameters'];assert p['hot_eligible_fraction']==1
    assert 0<=recovered_fraction<=1 and 0<=cold_dump_fraction<=1
    D=recovered_fraction*kinetic_captured_W
    cold=cold_dump_fraction*(1-recovered_fraction)*kinetic_captured_W
    removed=D+cold
    eta=p['heat_efficiency'];a=p['auxiliary_fraction'];pump=p['pump_fraction_hot']
    gross_delta=D-eta*removed
    delta=(1-a)*gross_delta+pump*removed-extra_direct_aux_W
    oldpower=thermal_case['net_W'];M=thermal_case['mass_kg']
    return {'captured_kinetic_W':kinetic_captured_W,'recovered_fraction':recovered_fraction,
        'direct_bus_W':D,'kinetic_cold_dump_W':cold,'removed_from_hot_engine_W':removed,
        'thermal_gross_electric_change_W':-eta*removed,'gross_electric_change_W':gross_delta,
        'extra_direct_aux_W':extra_direct_aux_W,'ideal_net_increment_W_before_new_massloads':delta,
        'same_specific_power_extra_mass_limit_kg_ignoring_massloads':delta/(oldpower/M),
        'one_mWkg_total_extra_mass_limit_kg_ignoring_massloads':(oldpower+delta)/.001-M,
        'scope':'Marginal ideal bound only. Additional peak-rated converter, buffer/export, support and maintenance mass/power must be debited. Existing driver storage is already allocated. Reduced exchanger/engine mass is not credited.'}
