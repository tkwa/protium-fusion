"""Weak distant-encounter guiding-center tensor with prescribed straight ions.

The electron magnetic moment is fixed in this limiting model. Unmagnetized
ion force correlations provide the impulse covariance; only parallel impulse
changes kinetic energy, while transverse impulses move the guiding center.
The scalar logarithm specifies an impact interval. Matching to finite-gyroradius
encounters and the source's inhomogeneous dielectric remains separate.
"""
from math import pi,sqrt,erf,exp
import numpy as np
from numba import njit
from electron_wall_orbits import ME,E_CHARGE as EV,EPS0


def species_table(species,logarithm):
    # G, me/mi, vth²; temperature follows from vth² and mass ratio.
    return np.array([[4*pi*(z*EV*EV/(4*pi*EPS0))**2*n*logarithm/ME**2,
        ME/mass,2*T*EV/mass] for mass,n,z,T in species])


@njit(cache=True)
def force_tensor_derivatives(w,G,vth2):
    """Relative-velocity parallel/per-component transverse impulse variance."""
    x=w*w/vth2
    if x<.1:
        D0=G*4/(3*sqrt(pi*vth2))
        # A common series also stabilizes derivatives near zero. Computing
        # erf(sqrt(x))-2*sqrt(x)*exp(-x)/sqrt(pi) loses precision there.
        sp=0.;st=0.;spp=0.;stp=0.;term=1.;dterm=0.
        for n in range(9):
            cp=3./(2*n+3);ct=cp/(2*n+1)
            sp+=cp*term;st+=ct*term;spp+=cp*dterm;stp+=ct*dterm
            dterm=(-term-x*dterm)/(n+1)
            term*=-x/(n+1)
        dp=D0*sp;dt=D0*st
        dpp=D0*spp*2*w/vth2;dtp=D0*stp*2*w/vth2
    else:
        psi=erf(sqrt(x))-2/sqrt(pi)*sqrt(x)*exp(-x)
        derivative=2/sqrt(pi)*sqrt(x)*exp(-x)
        dp=G*psi/(x*w)
        dt=G*((1-1/(2*x))*psi+derivative)/w
        dpp=G/(w*w)*(2*derivative-3*psi/x)
        dtp=G/(w*w)*((3/(2*x)-1)*psi-derivative)
    return dp,dt,dpp,dtp


@njit(cache=True)
def coefficients(vz,U,table):
    """Return physical impulse moments in velocity units, without free gyroturn.

    U is the ion mean y-velocity in the electric-drift frame. State variables
    are X, Y, vz, where X=x+vy/omega and Y=y-vx/omega. The returned diffusion
    tensor is for impulse/m_e (x,y,z); apply X=I_y/omega and Y=-I_x/omega.
    A_y and A_z include the Ito divergence and ion-polarization terms.
    """
    w=sqrt(U*U+vz*vz)
    dxx=0.;dyy=0.;dzz=0.;dyz=0.;dzzp=0.;dyzp=0.;ay=0.;az=0.
    for j in range(table.shape[0]):
        G,ratio,vth2=table[j]
        dp,dt,dpp,dtp=force_tensor_derivatives(w,G,vth2)
        if w>0.:
            s=vz/w;c=U/w;difference=dp-dt
            zz=dt+difference*s*s
            yz=-difference*c*s
            zzprime=(dtp+(dpp-dtp)*s*s)*s+difference*2*s*c*c/w
            yzprime=-U*((dpp-dtp)*vz*vz/(w*w*w)+difference*(U*U-vz*vz)/(w**4))
            yy=dt+difference*c*c
        else:
            zz=dt;yy=dt;yz=0.;zzprime=0.;yzprime=0.
        # me/(2*kTi) = (me/mi)/vth². D dot (v-U) = dp*(v-U).
        polarization=ratio/vth2*dp
        ay+=.5*yzprime+polarization*U
        az+=.5*zzprime-polarization*vz
        dxx+=dt;dyy+=yy;dzz+=zz;dyz+=yz;dzzp+=zzprime;dyzp+=yzprime
    return np.array([ay,az,dxx,dyy,dzz,dyz,dzzp,dyzp])


def state_moments(vz,U,table,omega):
    """Ito drift/diffusion in (X,Y,vz), including radial/axial correlation."""
    ay,az,dxx,dyy,dzz,dyz,dzzp,dyzp=coefficients(vz,U,table)
    A=np.array([ay/omega,0.,az])
    D=np.array([[dyy/omega**2,0.,dyz/omega],[0.,dxx/omega**2,0.],[dyz/omega,0.,dzz]])
    divergence=np.array([dyzp/omega,0.,dzzp])
    return A,D,divergence


def main():
    from pathlib import Path
    import argparse,json
    from scipy.integrate import quad
    from electron_wall_orbits import Cylinder
    from lepton_collision_screen import local_state
    root=Path(__file__).parent
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=root)
    ap.add_argument('--outdir',type=Path,default=root);ap.add_argument('--log',type=float,default=4.)
    args=ap.parse_args()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049'])
    cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])
    radius=float(np.interp(.5,cdf,c.eq['x'])*c.scale)
    local=local_state(c,radius);omega=-EV*c.B/ME
    table=species_table(local['species'],args.log)
    energy=np.r_[0.,np.geomspace(.01,2e4,257)]
    rows=[];averages=[]
    for flow in [False,True]:
        U=local['relative_azimuthal_m_s'] if flow else 0.
        out=[]
        for E in energy:
            vz=sqrt(2*E*EV/ME);q=coefficients(vz,U,table)
            out.append({'axial_energy_eV':float(E),'parallel_velocity_m_s':vz,
                'Ay_m_s2':float(q[0]),'Az_m_s2':float(q[1]),
                'impulse_covariance_m2_s3':[[float(q[2]),0.,0.],[0.,float(q[3]),float(q[5])],[0.,float(q[5]),float(q[4])]],
                'axial_heating_eV_s':float(ME/EV*(vz*q[1]+.5*q[4])),
                'diffusion_only_heating_eV_s':float(.5*ME/EV*q[4]),
                'electron_radial_drift_m_s':float(q[0]/omega)})
        rows.append({'include_relative_flow':flow,'states':out})
        for temperature in [10.,100.,1000.,2311.7367,10000.]:
            sigma=sqrt(temperature*EV/ME)
            def moment(z,which):
                vz=z*sigma;q=coefficients(vz,U,table)
                value=ME/EV*(vz*q[1]+.5*q[4]) if which==0 else q[0]/omega
                return sqrt(2/pi)*exp(-.5*z*z)*value
            p=quad(lambda z:moment(z,0),0,12,epsabs=1e-5,epsrel=3e-9)[0]
            drift=quad(lambda z:moment(z,1),0,12,epsabs=1e-13,epsrel=3e-9)[0]
            averages.append({'include_relative_flow':flow,'prescribed_parallel_temperature_eV':temperature,
                'mean_axial_heating_eV_s':p,'mean_electron_radial_drift_m_s':drift})
    output={'scope':__doc__,'logarithmic_interval':args.log,'B_T':c.B,
        'local_state':local,'cases':rows,'prescribed_Gaussian_averages':averages,
        'limits':['Magnetic moment fixed; no finite-gyroradius or close-encounter contribution.',
            'Straight ion orbits at zeroth order with perturbative ion recoil/polarization.',
            'The fixed logarithm is a scale parameter, not a qualified screening prescription.',
            'Gaussian axial distributions are prescribed probes, not stationary reactor solutions.',
            'No physical lepton birth rate, escape time or plant load is inferred.']}
    (args.outdir/'guiding-center-response.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'cold_with_flow':rows[1]['states'][0],'averages':averages},indent=2))


if __name__=='__main__':main()
