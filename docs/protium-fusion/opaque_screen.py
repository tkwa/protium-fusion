"""Radiation/confinement scales for the nuclear screen, not equilibria.

Frequency-dependent free-free absorption plus Thomson scattering is
Rosseland averaged. This is nonrelativistic, classical, fully ionized LTE.
The size estimate assumes a homogeneous sphere with dT/dr ~ -T/R.
Its numerical coefficient is a geometry convention, not a solved profile.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from nuclear_screen import screen, MU, KB, C, ARAD, MJ, Z, A

H = 6.62607015e-34
MU0 = 1.25663706212e-6
G = 6.67430e-11


def opacity(temp, rho, nuclear):
    f = np.array(list(nuclear['stationary_number_fractions'].values()))
    abar = float(f@A)
    zbar = float(f@Z)
    z2bar = float(f@(Z**2))
    xh, zcat = nuclear['X_H'], nuclear['X_CNO']
    yhe = 1-xh-zcat
    ne = rho/MU*(xh+yhe/2+zcat*zbar/abar)
    ni_z2 = rho/MU*(xh+yhe+zcat*z2bar/abar)
    x = np.geomspace(1e-4, 80, 4000)
    nu = x*KB*temp/H
    # Rybicki-Lightman alpha_ff: 3.692e8 in cgs; 3.692e-2 in SI.
    # A constant g_ff=1.2 is only an order-of-magnitude approximation.
    alpha_abs = 3.692e-2*temp**-0.5*ne*ni_z2*nu**-3*(-np.expm1(-x))*1.2
    kabs = alpha_abs/rho
    kes = nuclear['electron_scattering_opacity_m2_kg']
    weight = x**4*np.exp(-x)/(-np.expm1(-x))**2
    norm = simpson(weight,x=x)
    kr = norm/simpson(weight/(kes+kabs),x=x)
    ka4 = float(np.interp(4,x,kabs))
    return kr, ka4, kes


def scales(temp, rho, xh=.7, zcat=.1):
    n=screen(temp,rho,xh,zcat)
    kr,ka4,kes=opacity(temp,rho,n)
    eps=n['nuclear_heat_W_kg_active']
    # L_diff/M ~= 4 a c T^4/(kappa rho^2 R^2), with T-gradient T/R.
    r=np.sqrt(4*ARAD*C*temp**4/(kr*rho**2*eps))
    m=4*np.pi/3*rho*r**3
    pg=n['gas_pressure_Pa']; pr=n['LTE_radiation_pressure_Pa']
    pressure=pg+pr
    v=4*np.pi/3*r**3
    # A thin sphere has membrane mass 3PV/(2*specific_strength).
    # This is a virial-style mass floor, not a realizable pressure vessel:
    # the membrane must be thin and sustain local stress and temperature.
    strength=1e8
    membrane=1.5*pressure*v/strength
    # Both adiabatic gas and radiation contribute to restoring pressure.
    # Exact mixture gamma depends on gas/radiation energy exchange.
    cs_scale=np.sqrt((5*pg/3+4*pr/3)/rho)
    return {**n,
        'Rosseland_opacity_m2_kg':float(kr),
        'absorption_opacity_at_hnu_4kT_m2_kg':float(ka4),
        'diffusion_balance_radius_scale_m':float(r),
        'sphere_mass_at_radius_scale_kg':float(m),
        'sphere_mass_fraction_Jupiter':float(m/MJ),
        'radial_transport_optical_depth':float(kr*rho*r),
        'effective_thermalization_depth_at_hnu_4kT':float(np.sqrt(3*ka4*(ka4+kes))*rho*r),
        'diffusion_time_scale_s':float(3*kr*rho*r*r/C),
        'sound_crossing_scale_s':float(r/cs_scale),
        'heat_replacement_over_sound_crossing':float(n['thermal_energy_over_nuclear_heat_s']/(r/cs_scale)),
        'equivalent_magnetic_pressure_field_T':float(np.sqrt(2*MU0*pressure)),
        'thin_membrane_mass_floor_at_1e8Jkg_strength_kg':float(membrane),
        'membrane_floor_fraction_Jupiter':float(membrane/MJ),
        'self_gravity_binding_over_thermal_energy':float(.6*G*m/r/n['gas_plus_LTE_radiation_energy_J_kg']),
        'extra_flags':['size is a gradient scale, not a solution of stellar structure',
                       'constant Gaunt factor', 'no degenerate opacity or screening',
                       'no conductive or convective transport',
                       'membrane formula invalid as a physical design at extreme pressure',
                       'magnetic pressure equivalence does not imply stable confinement']}


if __name__=='__main__':
    records=[scales(t,rho) for rho in [1e3,1e5,1e7]
             for t in [3e7,5e7,1e8,2e8,3e8]]
    Path(__file__).with_suffix('.json').write_text(json.dumps({'records':records},indent=2)+'\n')
    print('rho | T MK | R km | M/Mj | tau | B T | Mmem/Mj | heat/sound time')
    for r in records:
        print(f"{r['rho_kg_m3']:7.1g} {r['T_K']/1e6:5.0f} {r['diffusion_balance_radius_scale_m']/1e3:10.3g} "
              f"{r['sphere_mass_fraction_Jupiter']:10.3g} {r['radial_transport_optical_depth']:10.3g} "
              f"{r['equivalent_magnetic_pressure_field_T']:10.3g} {r['membrane_floor_fraction_Jupiter']:10.3g} "
              f"{r['heat_replacement_over_sound_crossing']:10.3g}")
