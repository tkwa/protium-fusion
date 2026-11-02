"""Necessary conditions for remote delivery of the saved opaque-pulse inflows.

No delivery efficiency, material lifetime, or successful fluid assembly is inferred.
Units are SI unless a field explicitly says otherwise. See REPORT.md.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from fast_eos import FastEOS

EPS0 = 8.8541878128e-12
SIGMA = 5.670374419e-8
KB = 1.380649e-23
MU = 1.66053906660e-27
R_DRIVER = 2e10
FILES = [('A', 128, ''), ('A', 512, ''), ('B1', 128, ''),
         ('B2', 128, ''), ('C13', 64, ''), ('He', 64, '-r30')]
COMPOSITION = {'A': (8/15, 2/15), 'B1': (8/16, 2/16),
               'B2': (8/14, 3/14), 'C13': (7/14, 2/14), 'He': (2/4, 1/4)}


def wquant(x, w, quantiles=(.01, .05, .5, .95, .99)):
    order = np.argsort(x)
    x, w = np.asarray(x)[order], np.asarray(w)[order]
    c = (np.cumsum(w) - .5*w) / np.sum(w)
    return dict(zip(map(str, quantiles), np.interp(quantiles, c, x).tolist()))


def charge_accelerator(speed, density, radius, surface_E, applied_E, length):
    qm = 3*EPS0*surface_E/(density*radius)
    specific_ke = .5*speed**2
    return {'maximum_charge_per_mass_C_kg': qm,
            'surface_maxwell_pressure_Pa': .5*EPS0*surface_E**2,
            'acceleration_m_s2': qm*applied_E,
            'minimum_length_m': specific_ke/(qm*applied_E),
            'minimum_voltage_V': specific_ke/qm,
            'maximum_pellet_radius_m': 3*EPS0*surface_E*applied_E*length/(density*specific_ke)}


def survival(speed, density, radius, temperature, absorptivity, heat_allowance,
             path=R_DRIVER):
    travel = path/speed
    radiative_heat = 3*absorptivity*SIGMA*temperature**4*travel/(density*radius)
    column = 8*density*radius*heat_allowance/(3*speed**2)
    return {'travel_s': travel,
            'radiative_absorbed_energy_J_kg': radiative_heat,
            'minimum_pellet_radius_m': 3*absorptivity*SIGMA*temperature**4*travel/(density*heat_allowance),
            'maximum_absorptivity': density*radius*heat_allowance/(3*SIGMA*temperature**4*travel),
            'maximum_gas_column_kg_m2_at_full_energy_deposition': column,
            'maximum_uniform_cavity_gas_mass_kg': 4*np.pi*path**2*column/3}


def isentrope(rho_start, T_start, ye, yi, rho_goal):
    eos = FastEOS()
    # EOS is fully ionized. The integral is a diagnostic only at grain-contact
    # densities, where bound electrons and Coulomb terms may matter strongly.
    def fun(logrho, logT):
        rho = np.exp(logrho)/1000  # cgs
        T = np.exp(logT[0])
        st = eos.state(rho, max(T, np.exp(eos.lt[0])*1.000001), ye, yi)
        dT_dlogrho = (st.p/rho - rho*st.de_drho)/st.c_v
        return [float(dT_dlogrho/T)]
    floor = float(np.exp(eos.lt[0])*1.00001)
    def event(logrho, logT): return logT[0]-np.log(floor)
    event.terminal = True
    sol = solve_ivp(fun, (np.log(rho_start), np.log(rho_goal)), [np.log(T_start)],
                    rtol=2e-9, atol=1e-10, events=event, max_step=.1)
    rho_end, T_end = float(np.exp(sol.t[-1])), float(np.exp(sol.y[0,-1]))
    st = eos.state(rho_end/1000, T_end, ye, yi)
    return {'requested_density_kg_m3': rho_goal, 'density_kg_m3': rho_end,
            'temperature_K': T_end, 'specific_internal_J_kg': float(st.e)*1e-4,
            'pressure_Pa': float(st.p)*.1,
            'sound_speed_m_s': float(np.sqrt(st.gamma1*st.p/(rho_end/1000)))*.01,
            'reached_requested_density': len(sol.t_events[0]) == 0,
            'qualification': 'Fully ionized ideal-ion/electron/radiation EOS only; not a physical cold-material isentrope.'}


def screen(stage, zones, tag):
    stem = f'compression-{stage}-{zones}{tag}'
    initial_path = ROOT/(stem+'-initial.json')
    burning_path = ROOT/(stem+'-burning.json')
    p = json.loads(initial_path.read_text())
    b = json.loads(burning_path.read_text())
    mass = b['mass_kg']
    w = np.asarray(p['mass_fraction_of_total'])
    m = mass*w
    r = np.asarray(p['outer_radius_m'])
    v = -np.asarray(p['velocity_m_s'])
    node_m = np.r_[.5*(m[:-1]+m[1:]), .5*m[-1]]
    kinetic = float(.5*np.dot(node_m, v*v))
    rho = np.asarray(p['rho_kg_m3'])
    temp = np.asarray(p['T_K'])
    # Cell representative positions and velocities, for source scheduling only.
    ri = np.r_[0., r[:-1]]
    vi = np.r_[0., v[:-1]]
    rc = .75*(r**4-ri**4)/(r**3-ri**3)
    vc = vi+(v-vi)*(rc-ri)/(r-ri)
    launch = -(R_DRIVER-rc)/vc
    # Finite-volume midpoint schedule excludes the measure-zero central point.
    # Each cell is spread over neighboring midpoint launch-time boundaries.
    edge = np.r_[launch[0]-.5*(launch[1]-launch[0]),
                 .5*(launch[1:]+launch[:-1]),
                 launch[-1]+.5*(launch[-1]-launch[-2])]
    dt = np.diff(edge)
    cell_K = .5*m*vc**2
    cell_K *= kinetic/cell_K.sum()  # preserve the original node kinetic ledger
    launch_power = cell_K/dt
    flight = -launch
    v_rms = np.sqrt(2*kinetic/mass)
    contact = []
    for solid_density in (100., 1000., 2000., 3000.):
        def earlier_density(t):
            rr = np.r_[0., r+v*t]
            return m/(4*np.pi/3*np.diff(rr**3))
        f = lambda t: np.max(earlier_density(t))-solid_density
        hi = max(r[-1]/v[-1], 1e-3)
        while f(hi) > 0: hi *= 2
        back = brentq(f, 0., hi, xtol=1e-15)
        rr = r+v*back
        densities = earlier_density(back)
        contact.append({'material_density_kg_m3': solid_density,
                        'ballistic_lookback_s': back,
                        'outer_radius_m': float(rr[-1]),
                        'uniform_solid_ball_radius_m': float((3*mass/(4*np.pi*solid_density))**(1/3)),
                        'minimum_cell_density_kg_m3': float(densities.min()),
                        'maximum_cell_density_kg_m3': float(densities.max()),
                        'mean_density_kg_m3': float(mass/(4*np.pi/3*rr[-1]**3)),
                        'qualification': 'Constant-velocity backward map; pressure acceleration omitted. First contact of densest cells, not a delivered fluid state.'})
    ye, yi = COMPOSITION[stage]
    isos = [isentrope(float(rho[0]), float(temp[0]), ye, yi, rho_goal)
            for rho_goal in (1000., 2000., 3000.)]
    scenarios = []
    for wall_T in (300., 500., 1000.):
        for E_applied, length in ((1e6, 1e8), (1e7, 1e9)):
            density, surface_E, heat_allowance = 1000., 1e9, 1e6
            acc = charge_accelerator(v_rms,density,1e-3,surface_E,E_applied,length)
            amax = acc['maximum_pellet_radius_m']
            surv = survival(v_rms,density,amax,wall_T,1.,heat_allowance)
            fast_acc = charge_accelerator(float(v.max()),density,1e-3,surface_E,E_applied,length)
            fast_surv = survival(float(v.max()),density,fast_acc['maximum_pellet_radius_m'],wall_T,1.,heat_allowance)
            scenarios.append({'wall_temperature_K': wall_T,
                              'maximum_speed_pellet_radius_m': fast_acc['maximum_pellet_radius_m'],
                              'maximum_speed_Planck_absorptivity_limit': fast_surv['maximum_absorptivity'],
                              'pellet_density_kg_m3': density,
                              'surface_field_V_m': surface_E,
                              'accelerator_gradient_V_m': E_applied,
                              'accelerator_length_m': length,
                              'allowable_absorbed_heat_J_kg': heat_allowance,
                              'maximum_accelerable_pellet_radius_m': amax,
                              'maximum_Planck_absorptivity_for_cold_flight': surv['maximum_absorptivity'],
                              'surface_maxwell_pressure_Pa': acc['surface_maxwell_pressure_Pa'],
                              'maximum_uniform_cavity_gas_mass_kg': surv['maximum_uniform_cavity_gas_mass_kg'],
                              'gas_mass_fraction_of_this_pulse': surv['maximum_uniform_cavity_gas_mass_kg']/mass})
    granular_preheat = []
    contact_state = isos[1]
    if contact_state['reached_requested_density']:
        specific_contact = contact_state['specific_internal_J_kg']
        sound_contact = contact_state['sound_speed_m_s']
        theta_max = float(np.arcsin(min(1., np.sqrt(specific_contact/(.5*v_rms*v_rms)))))
        for grain_radius in (1e-3, .01, .1, 1.):
            grain_density = 1000.
            total_cross_section = 3*mass/(4*grain_density*grain_radius)
            opaque_radius = np.sqrt(total_cross_section/np.pi)
            granular_preheat.append({
                'pellet_radius_m': grain_radius,
                'pellet_material_density_kg_m3': grain_density,
                'mean_projected_grain_optical_depth_one_radius_m': float(opaque_radius),
                'expansion_length_if_vaporized_at_that_radius_m': float(sound_contact*opaque_radius/v_rms),
                'maximum_vaporization_distance_for_300m_sound_expansion_m': float(300*v_rms/sound_contact),
                'contact_thermal_to_inflow_kinetic_energy_ratio': float(specific_contact/(.5*v_rms*v_rms)),
                'maximum_symmetric_jet_merge_half_angle_rad_if_all_contact_heat_is_from_merger': theta_max,
                'qualification': 'Opaque grains; independently expanding vapor at diagnostic contact sound speed. Collisional collective focusing, recombination, and radiation cooling omitted; not a universal focus bound.'})
    pointing=[]
    for footprint in (30., 300., 3000.):
        for last_correction_distance in (R_DRIVER, 1e8, 1e6):
            t = last_correction_distance/v_rms
            dv = footprint/t
            pointing.append({'footprint_tolerance_m': footprint,
                             'distance_after_last_correction_m':last_correction_distance,
                             'angular_error_rad': footprint/last_correction_distance,
                             'uncorrelated_transverse_velocity_error_m_s':dv,
                             'uncorrelated_longitudinal_velocity_error_m_s':dv,
                             'fractional_longitudinal_velocity_error':dv/v_rms,
                             'arrival_time_error_s':footprint/v_rms,
                             'equivalent_unfocused_proton_temperature_K':MU*dv*dv/KB})
    return {'stage':stage,'zones':zones,'tag':tag,
            'inputs':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest()
                      for x in (initial_path,burning_path)},
            'mass_kg':mass,'input_energy_J':b['initial_energy_J'],
            'kinetic_energy_J':kinetic,'kinetic_ledger_relative_error':kinetic/b['initial_kinetic_energy_J']-1,
            'input_internal_energy_J':b['initial_internal_energy_J'],
            'energy_rms_speed_m_s':v_rms,'maximum_speed_m_s':float(v.max()),
            'launch_schedule':{'source_radius_m':R_DRIVER,'mass_weighted_flight_time_quantiles_s':wquant(flight,m),
                               'energy_weighted_flight_time_quantiles_s':wquant(flight,cell_K),
                               'peak_kinetic_source_power_W':float(launch_power.max()),
                               'electrical_peak_at_assumed_30pct_conversion_W':float(launch_power.max()/.3),
                               'largest_source_bin_duration_s':float(dt.max()),
                               'smallest_source_bin_duration_s':float(dt.min()),
                               'qualification':'Kinematic time-of-flight schedule; centroid cell quadrature. No source or accelerator design; kinetic power only.'},
            'contact_backmaps':contact,'ideal_fully_ionized_contact_isentropes':isos,
            'cold_pellet_acceleration_survival_trade':scenarios,
            'granular_preheat_requirements':granular_preheat,'uncorrected_pointing_requirements':pointing}


def main():
    cases=[screen(*spec) for spec in FILES]
    shell_cases=[]
    for strength in (1e9,1e10,1e11):
        pressure=.5*EPS0*1e9**2
        density_core,density_shell=1000.,2000.
        thickness_ratio=pressure/(2*strength)
        shell_cases.append({'grain_surface_field_V_m':1e9,'electrostatic_pressure_Pa':pressure,
                            'shell_allowable_tensile_stress_Pa':strength,
                            'core_density_kg_m3':density_core,'shell_density_kg_m3':density_shell,
                            'minimum_shell_thickness_over_radius':thickness_ratio,
                            'thin_shell_mass_over_core_mass':3*density_shell*thickness_ratio/density_core,
                            'qualification':'Spherical membrane carrying the entire electric stress; core support, adhesion, field evaporation, defects and other pressure loads excluded.'})
    result={'scope':'Necessary remote driver delivery constraints, not a realizable driver or a burn result.',
            'load_bearing_conductive_shell_scenarios':shell_cases,
            'constants':{'epsilon0_F_m':EPS0,'sigma_W_m2_K4':SIGMA,'source_radius_m':R_DRIVER},
            'assumptions':['Straight line constant-speed flight after source acceleration.',
                           'Spherical rigid pellets; surface field sets their charge.',
                           'Radiation is isotropic blackbody cavity irradiation with specified Planck-mean absorptivity.',
                           'Cold-pellet emission neglected; heat allowance is an explicit scenario, not a measured latent heat.',
                           'Full deposition of intercepted stationary-gas kinetic energy; transparency would relax the gas bound.',
                           'No collisions, gravity, charged-pellet space charge, ablation, or active corrections in the trajectory map.'],
            'cases':cases}
    (HERE/'delivery-screen.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for c in cases:
        s=c['launch_schedule'];p=c['cold_pellet_acceleration_survival_trade'][-1]
        print(c['stage'],c['zones'],'v_rms',c['energy_rms_speed_m_s'],
              'flight50',s['mass_weighted_flight_time_quantiles_s']['0.5'],
              'peak_electric_W',s['electrical_peak_at_assumed_30pct_conversion_W'],
              'amax',p['maximum_accelerable_pellet_radius_m'],'epsmax_1000K',p['maximum_Planck_absorptivity_for_cold_flight'],flush=True)

if __name__=='__main__':main()
