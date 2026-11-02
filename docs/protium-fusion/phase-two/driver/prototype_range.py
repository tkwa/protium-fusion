"""Arithmetic for a chosen fictional open test range, not a performance prediction.

Q=2.6 is a narrative scenario input. All consequences are elementary accounting;
material survival, accelerator construction, nuclear yield and schedule are unproved.
"""
from pathlib import Path
import hashlib
import json
from math import pi, sqrt

HERE = Path(__file__).resolve().parent
SIGMA = 5.670374419e-8
EPS0 = 8.8541878128e-12
G = 6.67430e-11
GM_SUN = 1.3271244e20  # IAU 2015 nominal solar mass parameter, m3/s2
AU = 149597870700.0   # IAU exact astronomical unit, m
C = 299792458.0


def main():
    hot_path = HERE/'scaled-hot-A-64-R30.json'
    hot = json.loads(hot_path.read_text())
    mass, delivered = hot['mass_kg'], hot['initial_energy_J']
    speed = sqrt(2*hot['initial_kinetic_energy_J']/mass)
    R, distance_AU, Q = 1e9, 120., 2.6
    distance = distance_AU*AU
    eta_d, eta_h = .3, .4
    flight = R/speed
    nuclear = Q*delivered
    outward = delivered+nuclear
    solar_flux = 1361./distance_AU**2
    omega = sqrt(GM_SUN/distance**3)
    tide = 2*GM_SUN*R/distance**3
    grain_density, surface_E, applied_E, accelerator_L = 1000., 1e9, 1e7, 1e8
    pellet_a = 3*EPS0*surface_E*applied_E*accelerator_L/(grain_density*.5*speed**2)
    accel_time = 2*accelerator_L/speed
    # Unidirectional sunlight: intercepted area pi*a^2, mass 4*pi*a^3*rho/3.
    solar_heat = 3*solar_flux*flight/(4*grain_density*pellet_a)
    pressure = .5*EPS0*surface_E**2
    charge_mass = 3*EPS0*surface_E/(grain_density*pellet_a)
    total_charge = charge_mass*mass
    kinetic = hot['initial_kinetic_energy_J']
    # Uniform sphere: U = 3 Q^2 / (20 pi epsilon0 R). Geometry is illustrative.
    coulomb_barrier_radius = 3*total_charge**2/(20*pi*EPS0*kinetic)
    output = {
        'scope': 'Chosen fictional open-range scenario: arithmetic only, no measured or predicted deployment performance.',
        'source_hot_diagnostic_sha256': hashlib.sha256(hot_path.read_bytes()).hexdigest(),
        'scenario_inputs': {'range_launch_radius_m': R, 'heliocentric_distance_AU': distance_AU,
            'pulse_mass_kg': mass, 'delivered_energy_J': delivered,
            'chosen_fictional_nuclear_over_delivered_Q': Q, 'driver_efficiency': eta_d,
            'thermal_conversion_efficiency': eta_h, 'solar_irradiance_at_1AU_W_m2': 1361.},
        'energy_accounting': {
            'convention': 'Q is deposited nuclear energy / delivered pulse energy. The fluence uses immediate capture heat; O15 decay, neutrinos and later processing are separate. Recover the delivered input heat once.',
            'chosen_nuclear_energy_J': nuclear, 'delivered_plus_nuclear_J': outward,
            'electric_driver_input_J': delivered/eta_d,
            'electric_return_at_complete_collection_J': eta_h*outward,
            'net_electric_at_complete_collection_before_auxiliaries_J': eta_h*outward-delivered/eta_d,
            'breakeven_Q_at_complete_collection_before_auxiliaries': 1/(eta_d*eta_h)-1,
            'breakeven_Q_at_99percent_collection_before_auxiliaries': 1/(eta_d*eta_h*.99)-1,
            'isotropic_energy_fluence_at_range_J_m2': outward/(4*pi*R**2),
            'uniform_mass_column_at_range_kg_m2': mass/(4*pi*R**2),
            'fluence_over_assumed_10MJ_per_kg_sacrificial_enthalpy_kg_m2': outward/(4*pi*R**2)/1e7,
            'shield_qualification': 'The last value is a toy energy/enthalpy quotient, not an ablation depth or a validated shield. Spectrum, ion stopping, impulse, nonuniformity and surface transients remain unknown.',
            'maximum_isotropic_energy_intercept_fraction_for_collector_projected_area_1e12m2': 1e12/(4*pi*R**2),
        },
        'flight_and_source_requirements': {'rms_inflow_speed_m_s': speed,
            'rms_speed_flight_time_s': flight, 'photon_flight_time_s': R/C,
            'all_outward_energy_as_kinetic_characteristic_speed_m_s': sqrt(2*outward/mass),
            'corresponding_debris_time_s': R/sqrt(2*outward/mass),
            'pointing_30m_angular_error_rad': 30/R,
            'pointing_30m_uncorrelated_transverse_velocity_error_m_s': 30/flight,
            'pointing_30m_arrival_error_s': 30/speed,
            'qualification': 'Characteristic flight and ballistic sensitivities only; all-kinetic outward speed is illustrative. Arrival distributions and useful diagnostic/shield timing require a computed outgoing spectrum.'},
        'open_space_thermal_environment': {'solar_irradiance_W_m2': solar_flux,
            'black_isothermal_sphere_equilibrium_K': (solar_flux/(4*SIGMA))**.25,
            'assumed_radiation_geometry': 'Direct sunlight only. Nearby equipment and prior shots can dominate local radiation and residual gas.',
            'solar_absorbed_energy_per_kg_during_flight_for_black_maximum_grain_J_kg': solar_heat},
        'accelerator_scenario': {'grain_density_kg_m3': grain_density, 'surface_field_V_m': surface_E,
            'applied_gradient_V_m': applied_E, 'length_m': accelerator_L,
            'maximum_grain_radius_at_rms_speed_m': pellet_a,
            'constant_acceleration_time_from_rest_s': accel_time,
            'integrated_accelerating_voltage_V': applied_E*accelerator_L,
            'surface_maxwell_pressure_Pa': pressure,
            'thin_shell_mass_over_core_at_1GPa_and_shell_density2000_kg_m3': 3*2000*pressure/(2*grain_density*1e9),
            'qualification': 'An example necessary requirement; no demonstrated material or accelerator. Acceleration may occur outside the range launch radius. Film mass, core strength, hydrogen retention and composition must be counted.'},
        'charge_neutralization_requirement': {
            'initial_grain_charge_per_mass_C_kg': charge_mass,
            'total_charge_if_all_grains_remain_at_accelerating_charge_C': total_charge,
            'uniform_sphere_radius_where_Coulomb_energy_equals_initial_kinetic_m': coulomb_barrier_radius,
            'uniform_30m_sphere_Coulomb_over_initial_kinetic': coulomb_barrier_radius/30.,
            'maximum_residual_net_charge_fraction_for_1percent_kinetic_penalty_at_30m': sqrt(.01*30/coulomb_barrier_radius),
            'isolated_grain_field_self_energy_per_kg_J_kg': 1.5*EPS0*surface_E**2/grain_density,
            'charge_compensating_electron_mass_kg': total_charge/1.602176634e-19*9.1093837015e-31,
            'qualification': 'Unscreened uniform-sphere electrostatic example, not the real beam geometry. Charged grains cannot coast inward unchanged; a discharge or charge-neutralization system must act before convergence. Electron mass is small, but its transport, currents, grain heating and final trajectory corrections are unmodeled.'},
        'gravity': {'solar_tidal_acceleration_across_range_m_s2': tide,
            'tidal_displacement_during_rms_flight_m': .5*tide*flight**2,
            'pulse_self_gravity_acceleration_at_range_m_s2': G*mass/R**2,
            'uniform_30m_sphere_binding_J_per_kg': 3*G*mass/(5*30.),
            'binding_over_delivered_energy': (3*G*mass**2/(5*30.))/delivered,
            'circular_orbit_angular_frequency_s_inverse': omega,
            'orbital_frame_rotation_length_during_rms_flight_m': omega*flight*R,
            'qualification': 'Tidal and self-gravity perturbations are small over flight. Rotation of orbit-following nodes requires predictable aiming corrections; no gravity-supported full collector shell is assumed.'},
        'checks': {'uniform_sphere_Coulomb_field_integral': abs((total_charge**2/(40*pi*EPS0*30.)+total_charge**2/(8*pi*EPS0*30.))/(kinetic*coulomb_barrier_radius/30.)-1)<1e-14,
            'breakeven_reconstructs_zero': abs(eta_h*(1+(1/(eta_d*eta_h)-1))-1/eta_d)<1e-14,
            'full_capture_scenario_loses_electricity': eta_h*outward<delivered/eta_d,
            'fluence_integrates_to_energy': abs(outward/(4*pi*R**2)*(4*pi*R**2)/outward-1)<1e-14,
            'accelerator_work_matches_kinetic_energy': abs((3*EPS0*surface_E/(grain_density*pellet_a))*applied_E*accelerator_L/(.5*speed**2)-1)<1e-14},
        'primary_constant_source': 'https://arxiv.org/abs/1510.07674',
        'remaining_requirements': ['Measured partially ionized cold-material EOS and shock Hugoniot.',
            'Discrete pellet collision, vaporization, ionization and multidimensional liner formation.',
            'Full reaction network and delayed O15 decay heat/neutrino accounting.',
            'Source hardware, pellet coating/strength/retention, field emission and finite accelerator geometry.',
            'Open-range diagnostic capture, spectrum-dependent shielding and safe shot cadence.',
            'Other catalytic stages, separation and recovery before any closed hydrogen-to-helium cycle claim.']}
    assert all(output['checks'].values())
    (HERE/'prototype-range.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ['energy_accounting','gravity','open_space_thermal_environment']},indent=2))


if __name__ == '__main__': main()
