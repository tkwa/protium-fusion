"""Weak straight-ion pulse response of an initially stationary electron.

Exact linear response to a prescribed unscreened Coulomb force pulse in uniform
B. The rate screen adds a static outer cutoff and classical/quantum inner scale.
It is a cold-electron diagnostic, not a finite-energy collision operator.
"""
from pathlib import Path
from math import pi,sqrt
import argparse,json
import numpy as np
from scipy.special import k0,k1,erf
from numpy.polynomial.legendre import leggauss
from electron_wall_orbits import ME,E_CHARGE as E,EPS0
from lepton_collision_screen import ion_coefficients

ROOT=Path(__file__).parent
HBAR=1.054571817e-34


def pulse_fractions(a,mu2):
    """Impact-azimuth-averaged energy / total zero-field impulse energy.

    a=|omega_e|*b/u. mu2 is (u_z/u)^2, or its conditional average.
    The parallel impulse uses the zero-frequency force integral. The rotating
    transverse velocity uses the force Fourier transform at the gyrofrequency.
    """
    a=np.asarray(a);safe=np.maximum(a,1e-150)
    ak1=np.where(a<1e-8,1.,safe*k1(safe))
    ak0=np.where(a<1e-8,0.,safe*k0(safe))
    parallel=np.broadcast_to((1-mu2)/2,a.shape)
    perp=ak1**2*(1+mu2)/2+ak0**2*(1-mu2)
    return parallel,perp


def speed_density(u,vth,U):
    U=abs(U);u=np.asarray(u)
    if U<1e-10*vth:
        return 4/sqrt(pi)*u*u/vth**3*np.exp(-(u/vth)**2)
    return u/(sqrt(pi)*vth*U)*(np.exp(-((u-U)/vth)**2)-np.exp(-((u+U)/vth)**2))


def conditional_mu2(u,vth,U):
    """Mean squared angle to B, for a Maxwellian drifting perpendicular to B."""
    a=2*np.asarray(u)*abs(U)/vth**2
    safe=np.maximum(a,1e-100)
    return np.where(a<1e-3,1/3-a*a/45+2*a**4/945,
                    (1/np.tanh(safe)-1/safe)/safe)


def rate(mass,n,z,temp,U,B,bmax,nv=256,nb=128,fixed_log=None,inner_factor=1.):
    vth=sqrt(2*temp*E/mass);mu=ME*mass/(ME+mass)
    K=z*E*E/(4*pi*EPS0);omega=abs(E*B/ME)
    # Integrate speed to eight thermal widths beyond the mean; the omitted
    # Gaussian tail is negligible for these moment screens.
    nodes,weights=leggauss(nv);umax=abs(U)+8*vth
    u=(nodes+1)*umax/2;wu=weights*umax/2
    density=speed_density(u,vth,U);mu2=conditional_mu2(u,vth,U)
    if fixed_log is None:
        bmin=inner_factor*np.maximum(abs(K)/(mu*u*u),HBAR/(mu*u))
    else:
        bmin=np.full_like(u,bmax*np.exp(-fixed_log))
    logarithm=np.maximum(0.,np.log(bmax/bmin))
    nodesb,weightsb=leggauss(nb)
    # Values with bmin>bmax have zero integration weight.
    b=np.exp(np.log(np.minimum(bmin,bmax))[:,None]+(nodesb[None,:]+1)*logarithm[:,None]/2)
    a=omega*b/u[:,None]
    parallel,perp=pulse_fractions(a,mu2[:,None])
    lp=np.sum(parallel*weightsb[None,:],axis=1)*logarithm/2
    lt=np.sum(perp*weightsb[None,:],axis=1)*logarithm/2
    pref=4*pi*n*K*K/(ME*E)*wu*density/u
    unmag=float(np.sum(pref*logarithm));par=float(np.sum(pref*lp));trans=float(np.sum(pref*lt))
    nu,dp,dt=ion_coefficients(np.array([abs(U)]),mass,n,z,temp,15.)
    nrl=float(.5*ME*(dp[0]+2*dt[0])/E)
    return {'unmagnetized_heating_eV_s':unmag,'parallel_heating_eV_s':par,
        'perpendicular_heating_eV_s':trans,'total_heating_eV_s':par+trans,
        'magnetized_over_unmagnetized':(par+trans)/unmag if unmag else None,
        'parallel_fraction_of_magnetized':par/(par+trans) if par+trans else None,
        'NRL_fixed_log_15_cold_heating_eV_s':nrl,
        'speed_normalization':float(np.sum(wu*density)),
        'cutoff_excluded_speed_probability':float(np.sum(wu*density*(bmin>=bmax))),
        'B_T':B,'outer_cutoff_m':bmax,'inner_scale_multiplier':inner_factor,
        'fixed_logarithm':fixed_log,'speed_nodes':nv,'impact_log_nodes':nb}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    scales=json.loads((args.outdir/'magnetized-lepton-scales.json').read_text())
    cases=[]
    for row in scales['cases']:
        local=row['local_state'];B=row['electron_omega_rad_s']*ME/E
        for species in row['species']:
            for outer_fraction in [.1,1.]:
                for inner_factor in [1.,10.]:
                    result=rate(species['mass_kg'],species['density_m-3'],species['charge'],local['ion_temperature_eV'],
                        local['relative_azimuthal_m_s'],B,local['ion_static_screening_m']*outer_fraction,inner_factor=inner_factor)
                    result.update({'reaction_fraction_inside':row['reaction_fraction_inside'],'species':species['name'],
                        'outer_fraction_of_static_screening':outer_fraction})
                    cases.append(result)
    out={'scope':__doc__,'cases':cases,
        'limitations':['Initial electron is at rest in the local electric-drift frame; it does not retain the finite gyroradius or energy distribution of the escape simulation.',
            'Ion paths are straight and prescribed; finite recoil, gyroturn, drift/gradient evolution, collective dielectric screening, nonlinear slow encounters and radiative feedback during a pulse are omitted.',
            'Impact cutoffs are conventions. Multiplying the inner scale or changing the outer range probes that choice and is not a physical uncertainty interval.',
            'Parallel heating denotes the increase of electron axial kinetic energy from the weak pulse; it is not an equilibrium temperature relaxation rate.']}
    (args.outdir/'cold-magnetized-response.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in cases:
        if r['reaction_fraction_inside']==.5:print({k:v for k,v in r.items() if k in ['species','outer_fraction_of_static_screening','inner_scale_multiplier','unmagnetized_heating_eV_s','parallel_heating_eV_s','perpendicular_heating_eV_s','magnetized_over_unmagnetized','cutoff_excluded_speed_probability']})
