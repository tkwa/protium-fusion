"""Independent limits and convergence of the transverse clearing diagnostics.

These checks qualify frozen orbits and a rigid-centroid closure. They do not
qualify a self-consistent driven plasma, lepton capture, or electrical gain.
"""
from pathlib import Path
import argparse, json, gzip
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.special import roots_legendre
from ion_dipole_spectrum import trajectories
from rigid_ion_dipole import modes
from electron_wall_orbits import Cylinder
from nonneutral_screen import EPS0, E_CHARGE, MEV
from catalyst_penning import MP


def main(root, data):
    checks=[]; details={}
    def check(name, value, limit):
        checks.append(dict(name=name,metric=float(value),upper_limit=limit,passed=bool(value<limit)))
    def orbit(k,h,total,initial):
        co=np.array([[0.],[.5*k],[0.],[0.]])
        n=int(round(total/h))
        return trajectories(np.asarray(initial,float),1.,1.,0.,h,n,1,0.,1.,co,1e6)
    # Independent closed-form Lorentz trajectory with no radial potential.
    init=np.array([[.3,-.2,.4,.7]])
    pos,kin,err,_,_=orbit(0.,.13,130.,init)
    t=np.arange(pos.shape[1])*.13; u=complex(*init[0,2:]);z0=complex(*init[0,:2])
    exact=z0+u*(1-np.exp(-1j*t))/1j
    check('Exact free magnetic helix',np.max(np.abs(pos[0]-exact)),2e-11)
    check('Exact free magnetic speed',np.max(np.abs(kin[0]-abs(u)**2)),2e-12)
    k=.3;D=np.sqrt(1+4*k);lp=(-1+D)/2;lm=(-1-D)/2
    harmonic=[]
    for lam in [lp,lm]:
        errors=[]
        for h in [.01,.005]:
            p,_,e,_,_=orbit(k,h,40.,[[1.,0.,0.,lam]])
            t=np.arange(p.shape[1])*h;er=float(np.max(np.abs(p[0]-np.exp(1j*lam*t))))
            errors.append(er)
            check(f'Harmonic circular orbit lambda={lam}, h={h}',er,.002)
            check(f'Harmonic energy lambda={lam}, h={h}',e[0],1e-4)
        check(f'Second-order harmonic convergence lambda={lam}',abs(errors[0]/errors[1]-4),.02)
        harmonic.append(dict(frequency=lam,errors=errors))
    details['harmonic_orbits']=harmonic
    # Force an oscillator from rest with a Hann-enveloped resonant circular
    # field. Its final mechanical energy independently checks the factor 1/4
    # in the circular dipole fluctuation-response convention.
    T=100.;E=.01
    def rhs(t,y):
        z=y[0]+1j*y[1];v=y[2]+1j*y[3]
        f=E*np.sin(np.pi*t/T)**2*np.exp(1j*lp*t)
        a=-1j*v-k*z+f
        return [v.real,v.imag,a.real,a.imag]
    sol=solve_ivp(rhs,[0,T],[0,0,0,0],rtol=2e-11,atol=2e-13,max_step=.05)
    y=sol.y[:,-1];actual=.5*np.dot(y[2:],y[2:])+.5*k*np.dot(y[:2],y[:2])
    variance_plus=2/(lp*D);window_integral=T/2
    spectral_prediction=E**2*lp**2*variance_plus*window_integral**2/4
    check('Circular-response normalization from independently forced oscillator',abs(actual/spectral_prediction-1),2e-7)
    details['forced_harmonic']={'energy':actual,'spectral_prediction':spectral_prediction}
    p=next(p for p in json.loads((root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(p['selected_source_profiles']['22049']);Omega=E_CHARGE*c.B/MP;w=-c.eq['kappa']*c.eq['xi']
    # Compare the splitting algorithm with an adaptive ODE in the nonlinear
    # source. Both species use the same physical initial transverse velocity.
    for Z,M in [(1.,1.),(22.,c.eq['M'])]:
        initial=np.array([[4.,.7,.23/np.sqrt(M),-.41/np.sqrt(M)]])
        h=.005;n=4000
        pos,_,err,_,_=trajectories(initial,Z,M,w,h,n,10,c.eq['x'][0],c.eq['x'][1]-c.eq['x'][0],c.spline.c,c.src['xwall'])
        def fun(t,y):
            r=np.hypot(*y[:2]);fac=-Z/M*float(c.spline(r,1))/r+w*(Z/M+w);om=Z/M+2*w
            return [y[2],y[3],fac*y[0]+om*y[3],fac*y[1]-om*y[2]]
        ts=np.arange(pos.shape[1])*h*10
        ref=solve_ivp(fun,[0,n*h],initial[0],t_eval=ts,rtol=2e-12,atol=2e-13,max_step=.01)
        check(f'Nonlinear adaptive-ODE orbit Z={Z}',np.max(np.abs(pos[0]-ref.y[0]-1j*ref.y[1])),3e-5)
    # Rigid mutual-force coefficient via an independent finite-displacement
    # two-dimensional integration, rather than its analytic density overlap.
    rigid=json.loads((data/'rigid-ion-dipole.json').read_text())
    x=c.eq['x'];rhoC=c.eq['c'];cdf=CubicHermiteSpline(x,c.eq['ic'],x*rhoC)
    gl,gw=roots_legendre(512);rr=(gl+1)*c.R/2;wr=gw*c.R/2
    theta=2*np.pi*(np.arange(256)+.5)/256;ct=np.cos(theta);st=np.sin(theta)
    nH=c.ns*c.src['h0']*np.exp(-c.spline(rr/c.scale)-c.eq['ah']*(rr/c.scale)**2)
    def force(displacement):
        xx=rr[:,None]*ct+displacement;yy=rr[:,None]*st;rad=np.hypot(xx,yy)
        charge=2*np.pi*c.scale**2*c.ns*22*E_CHARGE*cdf(rad/c.scale)
        Ex=charge/(2*np.pi*EPS0*rad**2)*xx
        return float(E_CHARGE*np.sum(wr*rr*nH*np.mean(Ex,axis=1))*2*np.pi)
    finite=[]
    for d in [c.scale*.01,c.scale*.005]:
        val=(force(d)-force(-d))/(2*d)
        finite.append(val)
        check(f'Mutual force from displaced density integral, d={d}',abs(val/rigid['cross_force_derivative_N_m2']-1),2e-5)
    details['finite_offset_force_coefficients']=finite
    A=np.array(rigid['total_force_matrix_N_m2']);ml=np.array(rigid['line_mass_kg_m']);Om=np.array(rigid['cyclotron_rad_s'])
    for row in rigid['modes']:
        frequency=complex(row['omega_real_rad_s'],row['omega_imag_rad_s'])
        v=np.array(row['relative_centroid_amplitude_real'])+1j*np.array(row['relative_centroid_amplitude_imag'])
        matrix=np.diag(ml*frequency*(frequency+Om))+A
        check(f'Quadratic eigenproblem residual {frequency}',np.linalg.norm(matrix@v)/(np.linalg.norm(A)+np.linalg.norm(ml*frequency*(frequency+Om))),2e-13)
    # Kohn limit: equal q/m and a uniform circular image field permit rigid
    # common translation. Internal forces cancel exactly; the other mode is dark.
    charges=np.array([1.,2.])*1e-6;masses=charges/2e-6;cycl=np.array([2.,2.])
    roots,Ak,im,cr=modes(masses,charges,cycl,.1,1.)
    coeff=np.sum(charges)*2e-6/(2*np.pi*EPS0)
    for exact in np.roots([1.,2.,coeff]):
        check('Equal-q/m common mode analytic frequency',min(abs(r['omega_real_rad_s']-exact) for r in roots),2e-13)
    check('Internal forces cancel in total centroid equation',np.linalg.norm(np.sum(cr,axis=0)),1e-14)
    for frequency in [-.4,-1.3]:
        response=np.linalg.solve(-np.diag(masses*frequency*(frequency+cycl))-Ak,charges)
        check('Uniform drive excites common translation at equal q/m',abs(response[0]/response[1]-1),2e-13)
    spectra={name:json.load(gzip.open(data/filename,'rt')) for name,filename in [
        ('reference','ion-dipole-spectrum.json.gz'),('half_step','ion-dipole-spectrum-half-step.json.gz'),('independent','ion-dipole-spectrum-independent.json.gz'),
        ('refined','ion-dipole-spectrum-refined.json.gz'),('refined_independent','ion-dipole-spectrum-refined-independent.json.gz')]}
    for name,s in spectra.items():
        for row in s['rows']:
            prefix=f'{name}, Z={row["Z"]}'
            check('Spectral Parseval '+prefix,abs(row['full_PSD_variance_m2']/row['weighted_time_mean_r2_m2']-1),2e-13)
            check('Independent velocity second spectral moment '+prefix,abs(row['full_PSD_velocity_m2_s2']/row['window_derivative_identity_m2_s2']-1),3e-5)
            check('Canonical velocity variance '+prefix,abs(row['time_mean_u2_scaled']*row['mass_over_proton']/2-1),.01)
            check('Stationary spatial variance '+prefix,abs(row['time_mean_r2_m2']/row['initial_mean_r2_m2']-1),.006)
            check('Rotating energy error '+prefix,row['max_rotating_energy_error_over_T'],.0005)
            check('No sampled wall crossings '+prefix,row['sampled_wall_crossings'],.5)
    bands=[]
    for i in range(2):
        ref=spectra['reference']['rows'][i];freq=np.array(ref['rotating_frequency_rad_s'])
        for lo,hi in [(-1.57e9,-1.4e9),(-1.4e9,-.8e9),(-.8e9,-.14e9)]:
            select=(freq+spectra['reference']['ion_rotation_rad_s']>=lo)&(freq+spectra['reference']['ion_rotation_rad_s']<hi)
            values={name:float(np.sum(np.array(s['rows'][i]['bare_rotating_heat_MeV_per_primary_per_E2'])[select])*s['frequency_resolution_rad_s']) for name,s in spectra.items()}
            # These are broad-band numerical comparisons; narrow-bin spectra
            # have larger ensemble fluctuations and no damping is inferred.
            ratio_step=abs(values['half_step']/values['reference']-1);ratio_seed=abs(values['refined_independent']/values['refined']-1)
            check(f'Broad-band orbit step refinement Z={ref["Z"]}, {lo,hi}',ratio_step,.003)
            check(f'Broad-band independent 16384-particle ensembles Z={ref["Z"]}, {lo,hi}',ratio_seed,.12)
            bands.append(dict(Z=ref['Z'],lab_frequency_band=[lo,hi],integrals=values,step_relative_change=ratio_step,
                seed_relative_change=ratio_seed,initial_4096_seed_relative_change=abs(values['independent']/values['reference']-1),
                reference_population_refinement_relative_change=abs(values['refined']/values['reference']-1)))
        e0=ref['max_rotating_energy_error_over_T'];eh=spectra['half_step']['rows'][i]['max_rotating_energy_error_over_T']
        check(f'Nonlinear source energy second-order convergence Z={ref["Z"]}',abs(e0/eh-4),.05)
    details['broad_bands']=bands
    out=dict(scope=__doc__,passed=sum(c['passed'] for c in checks),total=len(checks),checks=checks,details=details)
    (data/'ion-dipole-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['passed','total']}))
    for c in checks:
        if not c['passed']: print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--data',type=Path,default=Path(__file__).parent);a=ap.parse_args();main(a.root,a.data)
