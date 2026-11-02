"""Closed CN-I/hot-CN-I stationary screen; SI public interface.

Requires pynucastro==3.1.0, numpy and scipy. The library's ReacLib rates
are in cm^3 mol^-1 s^-1 for captures and s^-1 for beta decay.

This is an equilibrium catalytic throughput calculation at fixed T, rho,
and hydrogen mass fraction. It omits screening, NO branches, hot breakout,
finite fuel depletion and hydrodynamics. It is not an ignition calculation.
"""
import json
from pathlib import Path
import numpy as np
import pynucastro as pyna

MU = 1.66053906660e-27
KB = 1.380649e-23
C = 299792458.0
SIGMA_SB = 5.670374419e-8
ARAD = 4 * SIGMA_SB / C
MEV = 1.602176634e-13
SIGMA_T = 6.6524587321e-29
ME = 9.1093837015e-31
HBAR = 1.054571817e-34
MJ = 1.89813e27

LIB = pyna.ReacLibLibrary()
NAMES = {
    'cp12': 'c12(p,g)n13', 'bn13': 'n13(,)c13',
    'np13': 'n13(p,g)o14', 'bo14': 'o14(,)n14',
    'cp13': 'c13(p,g)n14', 'np14': 'n14(p,g)o15',
    'bo15': 'o15(,)n15', 'npa15': 'n15(p,a)c12',
}
RATES = {key: LIB.get_rate_by_name(name) for key, name in NAMES.items()}
STATES = ['c12', 'n13', 'c13', 'o14', 'n14', 'o15', 'n15']
A = np.array([12, 13, 13, 14, 14, 15, 15.0])
Z = np.array([6, 7, 6, 8, 7, 8, 7.0])


def screen(temp, rho, x_h=0.7, x_cat=0.1):
    """Return model diagnostics. rho is kg/m^3; T is kelvin."""
    x_he = 1 - x_h - x_cat
    if min(x_h, x_cat) <= 0 or x_he < 0:
        raise ValueError('Invalid mass fractions')
    k = {name: rate.eval(temp) * (1 if name.startswith('b') else rho / 1000 * x_h)
         for name, rate in RATES.items()}
    hot = k['np13'] / (k['bn13'] + k['np13'])
    dwell = np.array([1/k['cp12'], 1/(k['bn13']+k['np13']),
                      (1-hot)/k['cp13'], hot/k['bo14'],
                      1/k['np14'], 1/k['bo15'], 1/k['npa15']])
    cycle = dwell.sum()
    fractions = dwell / cycle
    abar = fractions @ A
    zbar = fractions @ Z
    z2bar = fractions @ (Z**2)
    catalysts_per_kg = x_cat / (abar * MU)
    # Total atomic mass-defect energy, including positron annihilation.
    # Approximate mean neutrino energies: N13 .706, O14 1.05, O15 1.00 MeV.
    # Verify against full beta spectra before treating these as precision data.
    q_total = 26.73 * MEV
    q_nu = ((1-hot)*0.706 + hot*1.05 + 1.00) * MEV
    eps = catalysts_per_kg / cycle * (q_total-q_nu)
    ye = x_h + x_he/2 + x_cat*zbar/abar
    ni_per_mass = (x_h + x_he/4 + x_cat/abar)/MU
    ne_per_mass = ye/MU
    pgas = rho * KB * temp * (ni_per_mass + ne_per_mass)
    prad = ARAD * temp**4 / 3
    ugas = 1.5 * pgas / rho
    urad = 3 * prad / rho
    # Nonrelativistic, nondegenerate optically thin free-free diagnostic,
    # bolometric Gaunt factor 1.2. 1.426e-27 cgs -> 1.426e-40 SI.
    ne = rho*ne_per_mass
    sum_ni_z2 = rho/MU*(x_h+x_he+x_cat*z2bar/abar)
    brem = 1.426e-40 * 1.2 * temp**0.5 * ne * sum_ni_z2 / rho
    ef = HBAR**2/(2*ME) * (3*np.pi**2*ne)**(2/3)
    return {
        'T_K': temp, 'rho_kg_m3': rho, 'X_H': x_h, 'X_CNO': x_cat,
        'cycle_s': float(cycle), 'hot_branch_fraction': float(hot),
        'catalyst_A_mean': float(abar), 'stationary_number_fractions': dict(zip(STATES, fractions.tolist())),
        'nuclear_heat_W_kg_active': float(eps),
        'neutrino_fraction': float(q_nu/q_total),
        'brem_thin_W_kg_active': float(brem), 'brem_over_fusion': float(brem/eps),
        'gas_pressure_Pa': float(pgas), 'LTE_radiation_pressure_Pa': float(prad),
        'gas_plus_LTE_radiation_energy_J_kg': float(ugas+urad),
        'thermal_energy_over_nuclear_heat_s': float((ugas+urad)/eps),
        'hydrogen_depletion_local_s': float(x_h*cycle/(4*MU*catalysts_per_kg)),
        'electron_scattering_opacity_m2_kg': float(SIGMA_T*ne_per_mass),
        'Fermi_energy_over_kT': float(ef/(KB*temp)),
        'max_total_over_active_mass_at_40pct_1mWkg': float(.4*eps/.001),
        'flags': ['fixed background steady catalyst distribution', 'closed CN-I only',
                  'no screening', 'LTE radiation assumption', 'no confinement or transport solution',
                  'no fuel feed/ash removal budget', 'ideal nondegenerate gas approximation'],
    }


def main():
    records = [screen(t, rho) for rho in [1, 1e3, 1e5, 1e7]
               for t in [1.5e7, 2e7, 3e7, 5e7, 1e8, 2e8, 3e8]]
    payload = {'model':'closed CN-I/hot-CN-I screening calculation',
               'pynucastro_version':pyna.__version__,
               'rate_ids':{key:str(rate.id) for key,rate in RATES.items()},
               'records':records}
    Path(__file__).with_name('nuclear_screen.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('rho kg/m3 | T MK | heat W/kg active | thin brem/fusion | cycle s | u/epsilon s | EF/kT')
    for r in records:
        print(f"{r['rho_kg_m3']:9.1g} {r['T_K']/1e6:6.0f} {r['nuclear_heat_W_kg_active']:13.3g} "
              f"{r['brem_over_fusion']:13.3g} {r['cycle_s']:12.3g} "
              f"{r['thermal_energy_over_nuclear_heat_s']:12.3g} {r['Fermi_energy_over_kT']:8.3g}")


if __name__ == '__main__':
    main()
