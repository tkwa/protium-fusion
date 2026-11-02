"""Diagnostic trace-lepton collision/radiation balance, not a reactor model.

Uses the unmagnetized Landau test-particle coefficients, a prescribed Coulomb
logarithm, an isotropic energy distribution in the local E x B frame, and
free-space cyclotron emission. The actual source is strongly magnetized;
neither this collision operator nor the isotropic closure is qualified there.
No lepton abundance, residence time or plant power is inferred by this script.
"""
from pathlib import Path
from math import pi, sqrt
import argparse,json,gzip
import numpy as np
from scipy.special import gammainc
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import solve_banded
from numpy.polynomial.legendre import leggauss
from electron_wall_orbits import Cylinder,ME,MP,E_CHARGE,EPS0,C
from nonneutral_screen import MEV

ROOT=Path(__file__).parent
EV=E_CHARGE


def ion_coefficients(speed,mass,density,Z,T_eV,log_lambda=15.):
    """Velocity drift frequency and parallel/per-component transverse diffusion.

    Diffusion entries mean <dv_i dv_j>/dt (no factor of two absorbed).
    The parallel direction is the relative velocity to one Maxwellian species.
    """
    speed=np.asarray(speed,dtype=float)
    vth2=2*T_eV*EV/mass
    x=speed*speed/vth2
    psi=gammainc(1.5,x)
    dpsi=2/np.sqrt(pi)*np.sqrt(x)*np.exp(-x)
    G=4*pi*(Z*E_CHARGE**2/(4*pi*EPS0))**2*density*log_lambda/ME**2
    safe=np.maximum(speed,1e-80)
    xs=np.maximum(x,1e-160)
    nu=G*(1+ME/mass)*psi/safe**3
    parallel=G*psi/(xs*safe)
    perp=G*((1-1/(2*xs))*psi+dpsi)/safe
    small=x<1e-7
    # Series limits avoid 0/0 at zero relative speed. First two corrections
    # retain continuity and are tested against the exact functions.
    a=4/(3*sqrt(pi))
    nu0=G*(1+ME/mass)*a/vth2**1.5
    D0=G*a/sqrt(vth2)
    nu=np.where(small,nu0*(1-3*x/5+3*x*x/14),nu)
    parallel=np.where(small,D0*(1-3*x/5+3*x*x/14),parallel)
    perp=np.where(small,D0*(1-x/5+3*x*x/70),perp)
    return nu,parallel,perp


def coefficients(energy_eV,species,flow_m_s=0.,log_lambda=15.,nangle=64):
    """Angle-averaged energy drift/diffusion and momentum transfer.

    Each species supplies (mass_kg, number_density_m-3, charge, temperature_eV).
    All ions have the same mean velocity, flow_m_s along an arbitrary direction.
    Returns coefficients in eV/s and eV^2/s, plus force along that direction.
    """
    energy=np.atleast_1d(energy_eV)
    v=np.sqrt(2*energy*EV/ME)[:,None]
    mu,weights=leggauss(nangle);mu=mu[None,:];weights=weights/2
    U=flow_m_s
    w2=v*v+U*U-2*v*U*mu
    w=np.sqrt(np.maximum(w2,0.))
    vw=v*v-v*U*mu
    nu=np.zeros_like(w);dp=nu.copy();dt=nu.copy()
    for mass,density,z,temp in species:
        n,p,t=ion_coefficients(w,mass,density,z,temp,log_lambda)
        nu+=n;dp+=p;dt+=t
    projected=dt*v*v+(dp-dt)*vw*vw/np.maximum(w2,1e-160)
    A=(-ME*nu*vw+ME/2*(dp+2*dt))/EV
    B=ME**2*projected/EV**2
    force=-ME*nu*(v*mu-U)
    # A direction-decorrelation proxy, exact for fixed-speed Lorentz scattering.
    angular=(dp+2*dt-projected/(v*v))/(2*v*v)
    mean=lambda a:a@weights
    return {'A_eV_s':mean(A),'B_eV2_s':mean(B),'force_N':mean(force),
            'direction_rate_s-1':mean(angular)}


def stationary(energy_eV,coeff,tau_perp_s=np.inf):
    E=np.asarray(energy_eV)
    rad=2*E/(3*tau_perp_s)
    A=coeff['A_eV_s']-rad;B=coeff['B_eV2_s']
    # Integrate in log E, resolving the regular E^(1/2) low-energy density.
    integral=cumulative_trapezoid(2*A/B*E,np.log(E),initial=0.)
    logf=integral-np.log(B);logf-=np.max(logf)
    f=np.exp(logf);f/=np.trapezoid(f,E)
    avg=lambda a:float(np.trapezoid(f*a,E))
    cdf=cumulative_trapezoid(f,E,initial=0.)
    return {'density_per_eV':f,'cdf':cdf,'mean_energy_eV':avg(E),
        'radiation_eV_s':avg(rad),'collision_eV_s':avg(coeff['A_eV_s']),
        'force_N':avg(coeff['force_N']),
        'mean_direction_rate_s-1':avg(coeff['direction_rate_s-1']),
        'fixed_speed_spatial_diffusion_proxy_m2_s':avg(2*E*EV/ME/(3*coeff['direction_rate_s-1'])),
        'probability_above_10keV':float(1-np.interp(1e4,E,cdf)),
        'upper_edge_probability_density_per_eV':float(f[-1])}


def green_kubo_diffusion(E,coeff,f):
    """Diffusion in the isotropized homogeneous velocity process.

    Solve (nu_direction - L_energy) g = v. The stationary energy generator
    is L_energy g = (B f g')'/(2 f). Then D_parallel = <v g>/3.
    A finite-volume discretization preserves the positive quadratic form.
    Energy-angle correlations in the real drifting plasma are not retained.
    """
    E=np.asarray(E);f=np.asarray(f)
    use=f>f.max()*1e-24
    last=np.where(use)[0][-1]+1
    E=E[:last];f=f[:last];B=coeff['B_eV2_s'][:last]
    nu=coeff['direction_rate_s-1'][:last]
    # Geometric cell faces on the logarithmic energy mesh.
    faces=np.r_[E[0]*sqrt(E[0]/E[1]),np.sqrt(E[:-1]*E[1:]),E[-1]*sqrt(E[-1]/E[-2])]
    widths=np.diff(faces);weight=f*widths;weight/=weight.sum()
    conduct=.5*np.sqrt(B[:-1]*f[:-1]*B[1:]*f[1:])/np.diff(E)
    conduct/=np.sum(f*widths)
    diagonal=nu*weight+np.r_[0.,conduct]+np.r_[conduct,0.]
    bands=np.zeros((3,len(E)));bands[1]=diagonal
    bands[0,1:]=-conduct;bands[2,:-1]=-conduct
    v=np.sqrt(2*E*EV/ME)
    g=solve_banded((1,1),bands,weight*v)
    diffusion=float(np.sum(weight*v*g)/3)
    residual=diagonal*g-np.r_[conduct*g[1:],0.]-np.r_[0.,conduct*g[:-1]]-weight*v
    return {'isotropized_green_kubo_diffusion_m2_s':diffusion,
        'diffusion_linear_residual_relative':float(np.linalg.norm(residual)/np.linalg.norm(weight*v)),
        'diffusion_energy_grid_last_eV':float(E[-1])}


def local_state(c,r):
    x=r/c.scale;eq=c.eq;z=c.src['ZA'][0]
    h=float(np.interp(x,eq['x'],eq['h']))*c.ns
    cat=float(np.interp(x,eq['x'],eq['c']))*c.ns
    mass=MP*eq['M']
    omega=c.src['xi']*z*E_CHARGE*c.B/mass
    Er=-1e6*c.T/c.scale*float(c.spline(x,1))
    drift=-Er/c.B
    ionflow=-omega*r
    return {'radius_m':r,'proton_density_m-3':h,'catalyst_density_m-3':cat,
        'ion_temperature_eV':c.T*1e6,'radial_E_V_m':Er,'ExB_azimuthal_m_s':drift,
        'ion_azimuthal_m_s':ionflow,'relative_azimuthal_m_s':ionflow-drift,
        'ion_static_screening_m':sqrt(EPS0*c.T*MEV/(E_CHARGE**2*(h+z*z*cat))),
        'species':[(MP,h,1,c.T*1e6),(mass,cat,z,c.T*1e6)]}


def case(c,r,E,log_lambda=15.,nangle=64,flow=True):
    local=local_state(c,r);U=local['relative_azimuthal_m_s'] if flow else 0.
    coeff=coefficients(E,local['species'],U,log_lambda,nangle)
    tau=3*pi*EPS0*ME**3*C**3/(E_CHARGE**4*c.B**2)
    result=stationary(E,coeff,tau)
    result.update(green_kubo_diffusion(E,coeff,result['density_per_eV']))
    force=result['force_N']
    # In a local Cartesian system x=radial, y=azimuthal, B=B z-hat,
    # the guiding-center force drift is F x B/(q B^2). Opposite charges
    # therefore drift in opposite radial directions under the same drag.
    result['electron_radial_force_drift_m_s']=force/(-E_CHARGE*c.B)
    result['positron_radial_force_drift_m_s']=force/(E_CHARGE*c.B)
    result['isotropic_radiation_time_s']=1.5*tau
    meanE=result['mean_energy_eV'];ve=sqrt(4*meanE*EV/(3*ME))
    rho=ve*ME/(E_CHARGE*c.B)
    result['rms_electron_gyroradius_m']=rho
    result['gyroradius_over_ion_screening']=rho/local['ion_static_screening_m']
    result['direction_rate_times_radiation_time']=result['mean_direction_rate_s-1']*1.5*tau
    result['gyroframe_mean_speed_over_c']=sqrt(2*meanE*EV/ME)/C
    result['log_lambda']=log_lambda;result['include_relative_flow']=flow
    result['local_state']={k:v for k,v in local.items() if k!='species'}
    for k in ['density_per_eV','cdf']:
        result[k]=result[k].tolist()
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);ap.add_argument('--grid',type=int,default=4097)
    args=ap.parse_args()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049'])
    E=np.geomspace(.001,5e4,args.grid)
    cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])
    rows=[]
    for fraction in [.01,.1,.5,.9,.99]:
        r=float(np.interp(fraction,cdf,c.eq['x'])*c.scale)
        for log_lambda in [5.,15.,30.]:
            for flow in [False,True]:
                row=case(c,r,E,log_lambda,flow=flow);row['reaction_fraction_inside']=fraction
                rows.append(row)
                print(fraction,log_lambda,flow,'mean eV',row['mean_energy_eV'],
                    'radiation MeV/s',row['radiation_eV_s']/1e6,
                    'electron radial m/s',row['electron_radial_force_drift_m_s'],
                    'angular/radiation',row['direction_rate_times_radiation_time'],flush=True)
    out={'scope':__doc__,'energy_eV':E.tolist(),'cases':rows,
         'limits':['Prescribed unmagnetized collision operator: log-lambda variation is a sensitivity, not a valid magnetization correction or confidence interval.',
           'Isotropic closure neglects velocity-space anisotropy, radial and axial transport, and finite source/sink terms.',
           'Free-space cyclotron emission omits cavity/absorption effects; relativistic corrections are absent.',
           'Force drift and fixed-speed spatial diffusion are local diagnostics, not measured or calculated residence times.',
           'No electron or positron abundance, atomic rate, instability growth or net plant cost follows from this screen.']}
    payload=json.dumps(out,separators=(',',':'),allow_nan=False).encode()
    (args.outdir/'lepton-collision-screen.json.gz').write_bytes(gzip.compress(payload,mtime=0))


if __name__=='__main__':main()
