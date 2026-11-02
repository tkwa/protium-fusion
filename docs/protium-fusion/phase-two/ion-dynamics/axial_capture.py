"""Relativistic single-wave capture bound and uniform-axial ion damping screen.

The electrostatic potential is prescribed, phi=Phi*cos(k*z-omega*t).
This is a separatrix availability bound, not an adiabatic capture probability.
The ion damping formula is exact only for a radially uniform total potential,
a Maxwellian parallel velocity distribution, and the linear collisionless model.
No collector, power circuit, or relativistic Maxwell structure is certified.
"""
from pathlib import Path
import json, argparse
import numpy as np
from scipy.special import wofz
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[2]
E=1.602176634e-19
C=299792458.
ME=9.1093837015e-31
MP=1.672621899e-27
MC2=ME*C*C/E
YEAR=365.25*86400

def source(root):
    p=next(x for x in json.loads((root/'network-plant.json').read_text())['cases'] if x['assumptions']['B']==140 and x['assumptions']['heat']==.7 and x['assumptions']['end_in']==.999)
    return p,p['selected_source_profiles']['22049']

def minimum_amplitude(beta, uz=0., uperp=0.):
    """Any-phase trapping requires 2 e Phi >= H0-Hmin; momenta are p/(mc)."""
    gp=np.sqrt(1+uperp*uperp)
    h0=np.sqrt(gp*gp+uz*uz)-beta*uz
    hmin=gp*np.sqrt(1-beta*beta)
    return MC2*(h0-hmin)/2

def initial_phase_fraction(beta,phi,uz=0.,uperp=0.):
    """Fraction of phases already inside a fixed sinusoidal bucket at switch-on.
    A particle's initial energy in the wave frame minus bucket minimum kinetic
    energy must fit Phi*(1-cos phase); charge sign translates the phase by pi.
    This excludes adiabatic passage, collisions and self-consistent wave growth.
    """
    d=2*minimum_amplitude(beta,uz,uperp)
    if phi<=0:return 0.
    a=d/phi-1
    if a>=1:return 0.
    if a<=-1:return 1.
    return float(np.arccos(a)/np.pi)

def uniform_chi(omega,k,T_eV,mass,Z,N):
    vt=np.sqrt(T_eV*E/mass)
    zeta=omega/(np.sqrt(2)*k*vt)
    return -N*(Z*E)**2/(T_eV*E)*(1+1j*np.sqrt(np.pi)*zeta*wofz(zeta))

def uniform_heat(beta,k,phi,s):
    omega=beta*C*k
    rows=[]
    for Z,M,N in [(1,1.,s['proton_inventory_per_m']),(22,48.58228894248968,s['catalyst_inventory_per_m'])]:
        chi=uniform_chi(omega,k,1.6e6,M*MP,Z,N)
        P=-omega/2*chi.imag*phi*phi
        rows.append(dict(Z=Z,heat_W_per_m=float(P),heat_MeV_per_primary=float(P/(s['reaction_rate_per_m_s']*1e6*E))))
    return rows

def main(root,out):
    plant,s=source(root);R=s['radius_m'];rows=[]
    for k in [10.,20.,40.,80.,160.]:
        for beta in [.01,.02,.03,.05,.075,.1,.15,.2,.225,.25,.275,.3,.35,.4,.5,.7,.9,.99]:
            phi=minimum_amplitude(beta)
            heat=uniform_heat(beta,k,phi,s)
            rows.append(dict(k_per_m=k,beta_phase=beta,minimum_total_potential_amplitude_V=phi,voltage_for_half_initial_phases_V=2*phi,
                titanium_qPhi_over_T=22*phi/1.6e6,proton_qPhi_over_T=phi/1.6e6,potential_over_axis_wall=phi/13.9433e6,
                omega_R_over_c=k*beta*R,frequency_Hz=beta*C*k/(2*np.pi),
                uniform_ion_heat_at_threshold_MeV_per_primary=sum(x['heat_MeV_per_primary'] for x in heat),species=heat,
                electrostatic_screen_note='Poisson-only propagation loses accuracy as beta_phase or omega*R/c approaches one.'))
    crossings=[]
    for k in [10.,20.,40.,80.,160.]:
        for duty in [1.,.1,.01]:
            for phase_fraction in [0.,.5,.9]:
                # phi=deltaH/[1+cos(pi f)]; zero is the limiting threshold.
                multiplier=2/(1+np.cos(np.pi*phase_fraction))
                def fun(b):
                    phi=multiplier*minimum_amplitude(b)
                    return duty*sum(x['heat_MeV_per_primary'] for x in uniform_heat(b,k,phi,s))-1.117140444
                b=brentq(fun,.15,.6)
                crossings.append(dict(k_per_m=k,duty=duty,initial_phase_fraction=phase_fraction,beta_equal_extra_margin=b,
                    phi_V=multiplier*minimum_amplitude(b),qualified_geometry=False))
    gamma_end=1+40e6/MC2;beta_end=np.sqrt(1-1/gamma_end**2)
    # Prescribed potential with synchronous phase sin(psi_s)=1/2 gives half peak force.
    stages=[]
    for phi in [100.,1000.,10000.]:
        k=40.;grad=.5*k*phi
        gamma0=1/np.sqrt(1-.1**2)
        length=(gamma_end-gamma0)*MC2/grad
        time=MC2*(np.sqrt(gamma_end**2-1)-np.sqrt(gamma0**2-1))/(C*grad)
        stages.append(dict(phi_V=phi,k_per_m=k,synchronous_gradient_V_m=grad,acceleration_distance_m=length,acceleration_time_s=time,
                          final_phase_beta=beta_end,final_bucket_bounce_omega_rad_s=np.sqrt(E*k*k*phi*np.sqrt(.75)/(ME*gamma_end**3)),
                          claim='Kinematics only; capture, radial Maxwell fields, loss and circuit power unresolved.'))
    data=dict(scope=__doc__,source=s,rows=rows,margin_crossings=crossings,kinematic_acceleration=stages,
      equations={'Hamiltonian':'H/(mc2)=sqrt(1+u_perp2+u_z2)-beta_phase*u_z+q*Phi*cos(phase)/(mc2)',
                 'minimum_kinetic_H':'sqrt(1+u_perp2)/gamma_phase',
                 'any_phase_bound':'2*e*Phi >= mc2*[sqrt(1+u_perp2+u_z2)-beta_phase*u_z-sqrt(1+u_perp2)/gamma_phase]',
                 'linear_uniform_ion_chi':'-N_line*q2/T * [1+zeta*Z(zeta)]',
                 'ion_heat':'-omega/2 * Im(chi) * Phi2'},
      limitations=['A fixed wave has zero initial-phase measure at the minimum trapping amplitude itself.',
                   'Half the initially slow particles lie inside a fixed bucket only at twice the minimum amplitude; 90% needs ~40.9 times minimum.',
                   'An accelerating/chirped wave can have different capture; the fixed-wave result does not prohibit staged acceleration.',
                   'No electrode wall voltage has been equated with the interior total potential.',
                   'Uniform total potential is a diagnostic geometry, not a proven lower bound over every radial profile.',
                   'Axisymmetric m=0 drive carries no canonical azimuthal torque; its ion heating is laboratory heating.',
                   'Heating could replace a compatible existing heater load, limited by the 8.668994 MeV main-stage thermal preparation budget.'],
      source_references=[{'title':'Kubo (1957), Statistical-Mechanical Theory of Irreversible Processes I','url':'https://doi.org/10.1143/JPSJ.12.570'},
                         {'title':'Esarey et al. (2007), Thermal effects in plasma-based accelerators','url':'https://doi.org/10.1063/1.2714022'}])
    out.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    for b in [.1,.2,.25,.3]:
        x=next(x for x in rows if x['k_per_m']==40 and x['beta_phase']==b)
        print(b,x['minimum_total_potential_amplitude_V'],x['uniform_ion_heat_at_threshold_MeV_per_primary'])
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=Path(__file__).with_suffix('.json'));a=ap.parse_args();main(a.root,a.out)
