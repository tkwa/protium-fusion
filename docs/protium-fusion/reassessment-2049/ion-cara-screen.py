"""Parametric CARA clearing screen; not a self-consistent accelerator design.

Run with Python 3. The sole output is ion-cara-screen.json beside this file.
Cold plane-wave kinematics, first-order radiation, and prescribed RF/capture
hardware are kept separate. No measured capture or future performance is input.
"""
from pathlib import Path
from math import pi, sqrt, exp, erf
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
E_CHARGE = 1.602176634e-19
ME = 9.1093837139e-31
EPS0 = 8.8541878188e-12
C = 299792458.0
MEV = E_CHARGE * 1e6
MC2_MEV = ME * C * C / MEV


def kinematics(B, radius, kinetic_MeV, electric_V_m):
    omega = E_CHARGE * B / ME
    acceleration = E_CHARGE * electric_V_m / (ME * C)
    u = sqrt(2 * kinetic_MeV / MC2_MEV)
    gamma = 1 + kinetic_MeV / MC2_MEV
    tau = 3 * pi * EPS0 * ME**3 * C**3 / (E_CHARGE**4 * B**2)
    pulse_s = u / acceleration
    time_s = (u + u**3 / 6) / acceleration
    length_m = C * u**3 / (6 * acceleration)
    energy_per_length = EPS0 * pi * radius**2 * electric_V_m**2
    # Radiation on the unperturbed orbit. This does not solve resonance
    # detuning: a chirp/taper must still be constructed with radiation reaction.
    radiation = MC2_MEV / (2 * tau * acceleration) * (u**3 / 3 + u**5 / 10)
    wavelength = 2 * pi * C / omega
    return dict(B_T=B, radius_m=radius, target_kinetic_MeV=kinetic_MeV,
        E_V_m=electric_V_m, cyclotron_Hz=omega/(2*pi), wavelength_m=wavelength,
        gamma=gamma, u_perp=u, u_parallel=gamma-1,
        gyro_radius_m=u*C/omega, transverse_pitch_tangent=u/(gamma-1),
        nonrelativistic_perpendicular_cooling_s=tau,
        particle_transit_s=time_s, retarded_pulse_s=pulse_s,
        acceleration_length_m=length_m,
        first_order_synchrotron_MeV=radiation,
        uniform_bore_mode_energy_J_m=energy_per_length,
        uniform_bore_peak_wave_power_W=energy_per_length*C,
        empty_pulse_J=energy_per_length*C*pulse_s,
        gaussian_Rayleigh_length_for_waist_equal_radius_m=pi*radius**2/wavelength,
        identity_time_residual_s=time_s-length_m/C-pulse_s)


def clearing(k, primary_per_m_s, *, length_m=1000, attenuation_m_inv=.01,
             capture_fraction=.1, pair_per_primary=.000906,
             lepton_radiation_MeV_s=100, RF_generation=.5,
             RF_pulses_per_attempt=2, total_RF_work_MeV_per_removed=120,
             heater=.95, heat=.6, network_guide_factor=1.):
    rate_unit = primary_per_m_s * MEV
    pulse = k['empty_pulse_J']
    # Prescribed maintained field: distributed feeds replace local loss.
    # Count full unrecovered wave output once at each physical cell exit.
    # Two pulses allow both directions; p is removal per COMPLETE attempt.
    rf_energy_J_m = network_guide_factor * RF_pulses_per_attempt * pulse * (attenuation_m_inv + 1/length_m) / RF_generation
    A = rf_energy_J_m / rate_unit
    cost_factor = 1/heater - heat
    # Exact mean population for constant births and independent removal p at
    # equally spaced attempts. A permanently inaccessible tail is separate.
    B = pair_per_primary * lepton_radiation_MeV_s * cost_factor * (1/capture_fraction-.5)
    f = sqrt(B/A)
    particle_work = pair_per_primary * total_RF_work_MeV_per_removed / RF_generation
    source_residence = (1/capture_fraction-.5)/f
    return dict(source_cell_length_m=length_m, network_guide_factor=network_guide_factor,
        power_attenuation_m_inv=attenuation_m_inv,
        attenuation_dB_m=attenuation_m_inv*10/2.302585092994046,
        one_end_only_power_transmission=exp(-attenuation_m_inv*length_m),
        capture_fraction_per_complete_attempt=capture_fraction,
        pair_supply_per_primary=pair_per_primary,
        positron_radiation_MeV_s=lepton_radiation_MeV_s,
        RF_generation_efficiency=RF_generation,
        RF_pulses_per_attempt=RF_pulses_per_attempt,
        total_RF_work_MeV_per_removed=total_RF_work_MeV_per_removed,
        extra_RF_J_per_active_m_per_attempt=rf_energy_J_m,
        A_MeV_s=A, B_MeV_per_s=B, optimized_attempts_s=f,
        mean_susceptible_positron_residence_s=source_residence,
        local_RF_field_duty=RF_pulses_per_attempt*k['retarded_pulse_s']*f,
        RF_empty_structure_MeV_per_primary=A*f,
        accumulated_positron_MeV_per_primary=B/f,
        injected_RF_particle_work_MeV_per_primary=particle_work,
        total_added_MeV_per_primary=2*sqrt(A*B)+particle_work,
        removed_particles_per_m_per_attempt=pair_per_primary*primary_per_m_s/f,
        successful_capture_RF_work_J_per_m_per_attempt=pair_per_primary*primary_per_m_s/f*total_RF_work_MeV_per_removed*MEV,
        requirement_cell_half_length_exceeds_acceleration=length_m/2 > k['acceleration_length_m'])


def main():
    p = json.loads((ROOT/'network-plant.json').read_text())
    # Radiator temperature is downstream of source-family selection; the saved
    # 60%-heat/30-K case has the same selected source as the 60%/50-K layout.
    ref = next(x for x in p['cases'] if x['assumptions']['B']==140 and x['assumptions']['heat']==.6 and x['assumptions']['end_in']==.999)
    ti = ref['selected_source_profiles']['22049']
    B = ti['B_T']; radius = ti['radius_m']; rate = ti['reaction_rate_per_m_s']
    # All fifteen processors need clearing. Equal amplitudes, schedules and cell
    # lengths give guide power proportional to total bore volume per primary.
    volume_per_primary = 0.
    for allocation in ref['channel_allocations']:
        key=str(1000*allocation['ZA'][0]+allocation['ZA'][1])
        s=ref['selected_source_profiles'][key]
        flow=allocation['pn_rate_s']/ref['primary_pn_rate_s']
        volume_per_primary += flow*pi*s['radius_m']**2/s['reaction_rate_per_m_s']
    network_factor = volume_per_primary/(pi*radius**2/rate)
    kin = [kinematics(B, radius, 60, E) for E in [1e5, 1e6, 1e7]]
    base = kin[1]
    cases = [clearing(base,rate,attenuation_m_inv=a,capture_fraction=p,network_guide_factor=network_factor)
             for a in [.001,.01,.1] for p in [.001,.01,.05,.1,.3]]
    selected = next(x for x in cases if x['power_attenuation_m_inv']==.01 and x['capture_fraction_per_complete_attempt']==.1)
    # Exact ideal monochromatic, detuned plane-wave bound: the driven transverse
    # momentum is at most2a/|delta|. Two non-overlapping resonance bands per
    # attempt are generously permitted. These are illustrative thermal states.
    bandwidth=[]
    omega=E_CHARGE*B/ME
    acceleration=E_CHARGE*base['E_V_m']/(ME*C)
    width=2*acceleration/(omega*base['u_perp'])
    for T_eV in [1e2,1e3,1e4]:
        sigma=sqrt(T_eV/(MC2_MEV*1e6))
        fraction=min(1.,2*erf(width/(sqrt(2)*sigma)))
        budget=clearing(base,rate,capture_fraction=fraction,network_guide_factor=network_factor)
        bandwidth.append(dict(illustrative_axial_temperature_eV=T_eV,
            fractional_resonance_halfwidth=width,
            maximum_fraction_two_optimally_placed_bands=fraction,
            corrected_budget=budget))
    # OLD70% reference only: rare inaccessible tail and full-reset tradeoff.
    # The new60% family needs a new inventory/reset converter closure before use.
    reset_A = 2055.3950472623355
    tail=[]
    for fraction in [1e-5,1e-4,1e-3]:
        rad_B=fraction*.000906*100*(1/.95-.7)
        period=sqrt(2*reset_A/rad_B)
        tail.append(dict(permanently_uncaptured_birth_fraction=fraction,
            full_reset_period_s=period,
            loss_only_minimum_added_MeV_per_primary=sqrt(2*reset_A*rad_B),
            scope='old70%-heat reference only, loss-only; NOT an allocation for the selected60% plant'))
    input_paths=[ROOT/'network-plant.json',ROOT/'phase-two/periodic-reset/periodic-reset.json',Path(__file__)]
    out=dict(scope=__doc__,constants=dict(e=E_CHARGE,me=ME,epsilon0=EPS0,c=C),
        input_sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in input_paths},
        source_heat_conversion=ref['assumptions']['heat'],
        source_profile_note='Saved60%-heat profile; radiator temperature does not enter source-family selection.',
        equal_field_full_network_guide_factor=network_factor,
        primary_rate_per_active_Ti_m_s=rate,
        primary_MeV_to_W_per_active_Ti_m=rate*MEV,
        kinematics=kin,clearing_cases=cases,selected_parametric_case=selected,
        monochromatic_bandwidth_discriminator=bandwidth,
        inaccessible_tail_screen=tail,
        warnings=[
          'Plane-wave equations assume initially resting leptons, uniform B and no source radial field. Mode dispersion, radiation reaction, source fields and capture need a kinetic Maxwell solve.',
          'The selected10% capture is an ACCEPTANCE TARGET for a new nondegenerate capture stage. The monochromatic bandwidth discriminator shows that the simple pulse does not achieve it in the illustrative thermal populations.',
          'Effective total RF work of120MeV per removed positron is a specified allowance for60MeV final energy, synchrotron radiation and unsuccessful excitation. It does not prove failed kicks fit the allowance.',
          'End launches alone cannot maintain the field at these attenuation lengths. Distributed injection must replace guide dissipation and beam loading; unabsorbed output pulse energy is explicitly discarded.',
          'The1km cell gives more than the ideal306m acceleration length to the farther end. Direction selection and acceleration near finite end fields are not modeled.',
          'No RF-to-electric heat recovery or recovered positron energy credit is assigned to RF losses. Accumulated lepton radiation gets the existing heater minus heat-engine credit.',
          'Peak source/coupler/switch ratings, waveform shaping, gas, residual electrons and selective heavy-product removal are outside this clearing sum.',
          'The particle supply is the recorded prompt-pair scenario, not a complete pair-production bound.',
          'Guide cost is extended to all processors by equal-field/equal-duty bore volume; stage-specific schedules may change this only after their particle handling is supplied.',
          'The permanently inaccessible fraction cannot be represented by independent repeated removal p; its separate reset screen prevents silently hiding a persistent population.'
        ])
    assert all(abs(k['identity_time_residual_s'])<1e-18 for k in kin)
    assert all(abs(x['RF_empty_structure_MeV_per_primary']-x['accumulated_positron_MeV_per_primary'])<1e-12 for x in cases)
    assert selected['requirement_cell_half_length_exceeds_acceleration']
    (HERE/'ion-cara-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(selected,indent=2));print(json.dumps(tail,indent=2))


if __name__=='__main__':main()
