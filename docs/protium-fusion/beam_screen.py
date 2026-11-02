"""Measured 14N(p,gamma) resonance versus neutral-target stopping.

Inputs: Gyurky et al., arXiv:1907.07446 (2019), abstract and section III.
Optimistically credit the entire CNO cycle energy to its nitrogen capture.
This rejects neither a hot plasma target nor recoverable nonneutral beams.
"""
import json
from pathlib import Path
from math import pi

HBARC_MEV_FM = 197.3269804
MP_MEV = 938.27208816
MEV_J = 1.602176634e-13


def case(energy_lab_keV, width_lab_keV, strength_eV, stop_eV_per_1e15atoms_cm2):
    reduced_mass = MP_MEV * 14/15
    energy_cm_MeV = energy_lab_keV / 1000 * 14/15
    # Optimistic convention while the quoted natural-width frame is checked:
    # treat it as a lab-energy width. If it is a CM width, peak sigma and
    # peak credit/stopping are smaller by 14/15; the conclusion is unchanged.
    width_cm_eV = width_lab_keV * 1000 * 14/15
    pi_over_k2_barn = pi*HBARC_MEV_FM**2/(2*reduced_mass*energy_cm_MeV)/100
    peak_sigma_barn = pi_over_k2_barn * 4*strength_eV/width_cm_eV
    # Upper bound: monoenergetic beam fixed precisely at line center,
    # no energy spread, and all later captures and decays supplied free.
    epsilon_eV_cm2 = stop_eV_per_1e15atoms_cm2 / 1e15
    full_cycle_credit_eV = 26.73e6
    ratio = peak_sigma_barn*1e-24*full_cycle_credit_eV/epsilon_eV_cm2
    # Thick-target resonance integral uses dE_lab, not dE_cm.
    integral_barn_eV_lab = 2*pi*pi_over_k2_barn*strength_eV*15/14
    yield_thick = integral_barn_eV_lab*1e-24/epsilon_eV_cm2
    return {'lab_energy_keV':energy_lab_keV,
            'strength_eV':strength_eV, 'width_lab_keV':width_lab_keV,
            'quoted_width_frame_status':'unverified; lab convention maximizes peak credit',
            'if_quoted_width_is_CM_multiply_peak_values_by':14/15,
            'peak_cross_section_barn':peak_sigma_barn,
            'optimistic_full_cycle_heat_over_stopping':ratio,
            'stopping_over_full_cycle_credit':1/ratio,
            'required_stopping_energy_recovery_fraction':1-ratio,
            'resonance_thick_target_yield_per_proton':yield_thick,
            'resonance_full_cycle_credit_over_incident_beam_energy':yield_thick*full_cycle_credit_eV/(energy_lab_keV*1000),
            'assumptions':['cold elemental nitrogen stopping estimate',
                           'zero beam energy spread at resonance for peak ratio',
                           'all later cycle steps free',
                           'total 26.73 MeV credited including otherwise escaping neutrinos',
                           'narrow isolated Breit-Wigner resonance',
                           'proton/N14 mass ratio approximated as 1/14'],
            'not_evaluated':['other energies and cycles','hot plasma targets',
                             'nonneutral ion beams','reversible recovery of scattered particles']}


if __name__ == '__main__':
    result={'source':'https://arxiv.org/abs/1907.07446',
            'cases':[case(278,1.12,.0134,10.72),case(1058,3.8,.442,4.733)]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
