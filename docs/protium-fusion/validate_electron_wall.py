"""Check invariant wall orbits against full relativistic Lorentz trajectories."""
from pathlib import Path
from math import pi,sqrt,acos
import argparse,json
import numpy as np
from scipy.integrate import solve_ivp
from electron_wall_orbits import Cylinder,ME,MEC2,E_CHARGE,EPS0,C
ROOT=Path(__file__).parent


def lorentz(c,K,mu,angle,sign,rtol=2e-12):
    length=ME*C/(E_CHARGE*c.B);time=ME/(E_CHARGE*c.B)
    p=sqrt(K*(K+2*MEC2))/MEC2;g0=1+K/MEC2
    init=np.array([c.R/length,0.,0.,-p*mu,p*sqrt(1-mu*mu)*np.cos(angle),p*sqrt(1-mu*mu)*np.sin(angle)])
    def rhs(t,y):
        r=np.hypot(y[0],y[1]);x=r*length/c.scale
        E=-1e6*c.T/c.scale*float(c.spline(x,1))
        v=y[3:]/sqrt(1+np.dot(y[3:],y[3:]));er=E/(c.B*C)
        force=sign*np.array([er*y[0]/r+v[1],er*y[1]/r-v[0],0.])
        return np.r_[v,force]
    def wall(t,y):return np.hypot(y[0],y[1])-c.R/length
    wall.direction=1;wall.terminal=True
    def turning(t,y):return y[0]*y[3]+y[1]*y[4]
    turning.direction=1;turning.terminal=False
    sol=solve_ivp(rhs,[0,8*pi*g0],init,method='DOP853',events=[wall,turning],
                  rtol=rtol,atol=1e-13,max_step=.04*g0,dense_output=True)
    if not sol.t_events[0].size or not sol.t_events[1].size:raise RuntimeError('Missing orbit event')
    final=sol.y_events[0][0];low=sol.y_events[1][0]
    sample=sol.sol(np.linspace(0,sol.t_events[0][0],301));r=np.hypot(sample[0],sample[1])*length
    gamma=np.sqrt(1+np.sum(sample[3:]**2,axis=0))
    energy=(gamma-1)*MEC2+sign*c.potential_MV(r)
    angular=sample[0]*sample[4]-sample[1]*sample[3]+sign*(sample[0]**2+sample[1]**2)/2
    # Perturbative Lienard radiation integral along the conservative orbit.
    # It measures the size of the omitted radiation reaction, not whether a
    # nearly grazing emitted particle is recaptured by a real surface.
    er=-1e6*c.T/c.scale*c.spline(r/c.scale,1)/(c.B*C)
    beta=sample[3:]/gamma
    ex=er*sample[0]*length/r;ey=er*sample[1]*length/r
    force2=(ex+beta[1])**2+(ey-beta[0])**2
    force2-=(beta[0]*ex+beta[1]*ey)**2
    tau_perp=3*pi*EPS0*ME**3*C**3/(E_CHARGE**4*c.B**2)
    times=np.linspace(0,sol.t_events[0][0],301)*time
    radiated=float(np.trapezoid(MEC2/(2*tau_perp)*gamma**2*force2,times))
    return {'return_time_s':sol.t_events[0][0]*time,'penetration_m':c.R-np.hypot(low[0],low[1])*length,
      'axial_displacement_m':final[2]*length,'max_energy_error_MeV':float(np.max(np.abs(energy-K))),
      'max_relative_angular_error':float(np.max(np.abs(angular-angular[0]))/abs(angular[0])),
      'perturbative_radiation_MeV':radiated,'perturbative_radiation_over_initial_energy':radiated/K}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);checks=[];comparisons=[]
    def check(name,passed,metric=None):checks.append({'name':name,'passed':bool(passed),'metric':metric})
    for sign in [-1,1]:
        for K in [1e-6,1e-3,.01,1.,30.,100.]:
            for mu,angle in [(.5,.7),(.9,3.4)]:
                ref=c.orbit(K,mu,angle,sign);ode=lorentz(c,K,mu,angle,sign)
                depth=abs(ref['penetration_m']-ode['penetration_m'])/ref['penetration_m']
                duration=abs(ref['return_time_s']-ode['return_time_s'])/ref['return_time_s']
                axial=abs(ref['axial_displacement_m']-ode['axial_displacement_m'])/max(abs(ref['axial_displacement_m']),1e-20)
                check(f'Lorentz/invariant orbit sign={sign} K={K:g} mu={mu}',depth<2e-4 and duration<2e-4 and axial<2e-4,{'relative_depth':depth,'relative_time':duration,'relative_axial':axial})
                check(f'Lorentz invariants sign={sign} K={K:g} mu={mu}',ode['max_energy_error_MeV']<max(2e-9,K*2e-9) and ode['max_relative_angular_error']<2e-10,ode)
                comparisons.append({'charge_sign':sign,'initial_MeV':K,'mu':mu,'azimuth':angle,'invariant':ref,'lorentz':ode})
    # The cold planar crossed-field solution is independent of the cylindrical
    # Hamiltonian. Curvature corrections should be small at sub-micron depth.
    cold=c.orbit(0.,.5,0.,-1)
    exact_depth=2*ME*c.src['wall_E_V_m']/(E_CHARGE*c.B**2)
    exact_time=2*pi*ME/(E_CHARGE*c.B)
    check('Cold planar crossed-field depth',abs(cold['penetration_m']/exact_depth-1)<5e-4,cold['penetration_m']/exact_depth-1)
    check('Cold planar gyroperiod',abs(cold['return_time_s']/exact_time-1)<5e-4,cold['return_time_s']/exact_time-1)
    # Turn off E for an exact circular orbit with constant relativistic mass.
    old=c.potential_by_depth_MV;c.potential_by_depth_MV=lambda d:np.asarray(d)*0.
    for sign in [-1,1]:
        for K in [.001,1.,30.]:
            mu=.65;angle=.8;p=sqrt(K*(K+2*MEC2));pp=p*sqrt(1-mu*mu)*np.cos(angle);pz=p*sqrt(1-mu*mu)*np.sin(angle)
            g=sqrt((c.R+pp/(sign*c.kB))**2+(p*mu/c.kB)**2)
            rho=sqrt((p*mu)**2+pp*pp)/c.kB
            cosang=(c.R*c.R-g*g-rho*rho)/(2*g*rho)
            duration=(2*pi-2*acos(np.clip(cosang,-1,1)))*(1+K/MEC2)*ME/(E_CHARGE*c.B)
            ref=c.orbit(K,mu,angle,sign)
            check(f'Exact zero-E circular orbit sign={sign} K={K:g}',
              abs(ref['penetration_m']/(c.R-g+rho)-1)<2e-7 and abs(ref['return_time_s']/duration-1)<2e-7,
              {'depth_relative':ref['penetration_m']/(c.R-g+rho)-1,'time_relative':ref['return_time_s']/duration-1})
    c.potential_by_depth_MV=old
    # Refine the underlying plasma potential, independently of orbit tolerance.
    fine=Cylinder(c.src,nout=24001,rtol=2e-12)
    for r in [0.004,0.0066,0.01,0.019]:
        for sign in [-1,1]:
            k=c.minimum_wall_energy(r,sign);kf=fine.minimum_wall_energy(r,sign)
            check(f'Potential refinement at radius {r:g}, sign={sign}',abs(k/kf-1)<2e-7,abs(k/kf-1))
            # Favorable nearly tangent emission at 1% above the necessary bound
            # must pass the target radius for this monotone reference profile.
            orbit=c.orbit(k*1.01,1e-4,0. if sign<0 else pi,sign,with_time=False)
            check(f'Accessibility just above energy bound r={r:g}, sign={sign}',orbit['minimum_radius_m']<r,orbit['minimum_radius_m'])
    out={'scope':__doc__,'passed':sum(v['passed'] for v in checks),'total':len(checks),'checks':checks,'lorentz_comparisons':comparisons}
    (args.outdir/'electron-wall-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'checks pass')
    for item in checks:
        if not item['passed']:print('FAILED',item)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
