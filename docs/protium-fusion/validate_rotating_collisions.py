"""Check the rotating-frame force and measure finite ion-path curvature.

The exact laboratory Lorentz equation and the independently transformed
rotating equation use the saved nonlinear radial electric field. These tests
check coordinates and ordering; no collisional force is integrated here.
"""
from pathlib import Path
from math import sin,cos,sqrt
import argparse,json,time
import numpy as np
from scipy.integrate import solve_ivp
from electron_wall_orbits import Cylinder,E_CHARGE as EV


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--data',type=Path,default=Path(__file__).parent);args=ap.parse_args();start=time.monotonic()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);scales=json.loads((args.data/'rotating-collision-scales.json').read_text());checks=[];results=[]
    def check(name,value,tol):checks.append({'name':name,'passed':bool(value<tol),'metric':float(value),'upper_limit':tol})
    def turn(x):return np.array([-x[1],x[0]])
    def rotate(x,angle):return np.array([[cos(angle),-sin(angle)],[sin(angle),cos(angle)]])@x
    def field(r):return -c.T*1e6/c.scale*float(c.spline(max(r/c.scale,c.eq['x'][0]),1))
    directions=[np.zeros(2),np.array([1.,0.]),np.array([-1.,0.]),np.array([0.,1.]),np.array([0.,-1.])]
    directions += [np.array([a,b])/sqrt(2) for a in [-1.,1.] for b in [-1.,1.]]
    for row in scales['rows']:
        r0=row['local_state']['radius_m'];w=row['ion_rotation_rad_s'];x=r0/c.scale
        for sp in row['species']:
            z=sp['charge'];mass=sp['mass_kg'];Omega=sp['laboratory_cyclotron_rad_s'];T=mass*sp['thermal_speed_m_s']**2/2
            # Independent density derivative from the saved Gibbs profile.
            coef=c.eq['ah'] if z==1 else c.eq['ac']
            exact=(-z*float(c.spline(x,1))-2*coef*x)/c.scale
            check(f'Rotating force equals thermal density gradient, q={row["reaction_fraction_inside"]}, Z={z}',abs(exact-sp['log_density_first_derivative_m-1'])/max(1,abs(exact)),1e-10)
            if row['reaction_fraction_inside'] not in [.5,.99]:continue
            for encounter in sp['encounters']:
                tau=encounter['time_s'];b=encounter['distance_m'];label=f'q={row["reaction_fraction_inside"]}, Z={z}, {encounter["electron_probe"]}, b/lambda={encounter["screening_length_fraction"]}'
                worst_position=worst_velocity=worst_invariant=0.;errors=[]
                for direction in directions:
                    u0=sp['thermal_speed_m_s']*direction;vrot=u0*tau/r0;vlab=vrot+w*tau*np.array([0.,1.])
                    def derivative(s,y,rotating):
                        xx=y[:2];vv=y[2:];radius=np.linalg.norm(xx)*r0
                        accel=z*EV/mass*field(radius)*xx/(radius/r0)*tau*tau/r0
                        if rotating:accel+=w*(Omega+w)*tau*tau*xx
                        accel-=(Omega+(2*w if rotating else 0))*tau*turn(vv)
                        return np.r_[vv,accel]
                    def invariant(y):
                        xx=y[:2];vv=y[2:]*r0/tau;radius=np.linalg.norm(xx)*r0
                        return .5*mass*np.dot(vv,vv)+z*EV*1e6*float(c.potential_MV(radius))-.5*mass*w*(Omega+w)*radius*radius
                    reference=invariant(np.r_[1.,0.,vrot])
                    for sign in [-1,1]:
                        ts=np.linspace(0,sign,17)
                        a=solve_ivp(lambda s,y:derivative(s,y,True),[0,sign],np.r_[1.,0.,vrot],t_eval=ts,rtol=2e-10,atol=2e-12,max_step=.05)
                        lab=solve_ivp(lambda s,y:derivative(s,y,False),[0,sign],np.r_[1.,0.,vlab],t_eval=ts,rtol=2e-10,atol=2e-12,max_step=.05)
                        if not a.success or not lab.success:raise RuntimeError('Orbit integration failed')
                        for si,yr,yl in zip(ts,a.y.T,lab.y.T):
                            transformed_position=rotate(yl[:2],-w*tau*si)
                            transformed_velocity=rotate(yl[2:]-w*tau*turn(yl[:2]),-w*tau*si)
                            worst_position=max(worst_position,np.linalg.norm(transformed_position-yr[:2])*r0/b)
                            worst_velocity=max(worst_velocity,np.linalg.norm(transformed_velocity-yr[2:])*r0/(tau*sp['thermal_speed_m_s']))
                            worst_invariant=max(worst_invariant,abs(invariant(yr)-reference)/T)
                        # Relative ion/electron-GC displacement; z motion is exactly linear.
                        actual_ion=a.y[:2,-1]*r0
                        actual_electron=r0*np.array([cos(row['relative_GC_rotation_rad_s']*tau*sign),sin(row['relative_GC_rotation_rad_s']*tau*sign)])
                        straight=(u0-np.array([0.,row['relative_GC_rotation_rad_s']*r0]))*tau*sign
                        errors.append(float(np.linalg.norm(actual_ion-actual_electron-straight)/b))
                check('Laboratory versus rotating position '+label,worst_position,1e-7)
                check('Laboratory versus rotating velocity '+label,worst_velocity,1e-7)
                check('Rotating energy invariant '+label,worst_invariant,1e-7)
                results.append({'label':label,'reaction_fraction_inside':row['reaction_fraction_inside'],'charge':z,'electron_probe':encounter['electron_probe'],'screening_length_fraction':encounter['screening_length_fraction'],'position_error_normalized_by_distance':worst_position,'velocity_error_normalized_by_thermal_speed':worst_velocity,'energy_invariant_error_normalized_by_ion_temperature':worst_invariant,'relative_path_departures_over_distance':errors,'maximum_sampled_relative_path_departure_over_distance':max(errors)})
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'orbit_comparisons':results,'sampling':'Initial ion position at the electron GC radius; transverse thermal velocities zero and eight directions at vth; endpoints at plus/minus one characteristic encounter time. No impact-offset ensemble or small-relative-speed tail is integrated.','seconds':time.monotonic()-start}
    (args.data/'rotating-collision-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['passed','total','seconds']}),flush=True)
    for c in checks:
        if not c['passed']:print('FAILED',c)
    for row in results:print(row['label'],'maximum relative path departure',row['maximum_sampled_relative_path_departure_over_distance'])
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
