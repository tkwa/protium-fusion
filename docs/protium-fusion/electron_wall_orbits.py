"""Relativistic wall-emission orbits in the saved cylindrical source fields.

Energy, canonical angular momentum and axial momentum are exact invariants of
the prescribed static, axisymmetric fields. Collisions, emission-driven fields,
roughness, ends, RF perturbations and radiation reaction are separate effects.
"""
from pathlib import Path
from math import pi,sqrt
import json,argparse
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.integrate import quad
from catalyst_penning import equilibrium,MP
from nonneutral_screen import MEV,E_CHARGE,EPS0,C
ME=9.1093837139e-31
MEC2=ME*C*C/MEV
ROOT=Path(__file__).parent


class Cylinder:
    def __init__(self,src,nout=12001,rtol=1e-11):
        self.src=src;self.R=src['radius_m'];self.B=src['B_T'];self.T=src['T_MeV']
        self.eq=eq=equilibrium(*src['ZA'],src['h0'],src['c0'],src['xi'],nout=nout,rtol=rtol)
        self.scale=sqrt(self.T*MEV*MP)/(E_CHARGE*self.B)
        self.ns=EPS0*self.B**2/MP
        self.spline=CubicHermiteSpline(eq['x'],eq['psi'],eq['dpsi'])
        self.wallpsi=float(self.spline(src['xwall']))
        self.kB=E_CHARGE*self.B*C/MEV  # MeV/c momentum per metre

    def potential_MV(self,r):
        return self.T*(self.spline(np.maximum(np.asarray(r)/self.scale,self.eq['x'][0]))-self.wallpsi)

    def potential_by_depth_MV(self,d):
        # Direct local polynomial subtraction avoids loss of precision for the
        # sub-micron magnetic-insulation orbit near the wall.
        xw=self.src['xwall'];dx=np.asarray(d)/self.scale
        idx=np.searchsorted(self.spline.x,xw,side='right')-1
        idx=min(idx,len(self.spline.x)-2)
        near=dx<=xw-self.spline.x[idx]
        y=xw-self.spline.x[idx];a,b,c,_=self.spline.c[:,idx]
        delta=-dx*(3*a*y*y+2*b*y+c)+dx*dx*(3*a*y+b)-a*dx**3
        far=self.potential_MV(self.R-np.asarray(d))
        return np.where(near,self.T*delta,far)

    def minimum_wall_energy(self,r,sign=-1):
        """Least initial kinetic energy over all emission directions to reach r.

        R*p_initial + r*p_inside must at least equal the magnetic-flux term.
        Axial momentum is zero for this minimum, and the wall direction is
        tangent in the favorable azimuthal sense. Reaching smaller radii must
        also pass any intervening barrier; the equation is a necessary bound.
        """
        phi=float(self.potential_MV(r));mag=self.kB*(self.R*self.R-r*r)/2
        def f(k):
            inside=k-sign*phi
            if inside<0:return -mag
            return self.R*sqrt(k*(k+2*MEC2))+r*sqrt(inside*(inside+2*MEC2))-mag
        lower=max(0.,sign*phi)
        if f(lower)>=0:return lower
        upper=max(1.,lower+1.)
        while f(upper)<0:upper*=2
        return brentq(f,lower,upper,xtol=1e-11)

    def maximum_access_depth(self,K,sign=-1):
        # Invert the necessary all-direction energy bound. This can overestimate
        # access if an intervening barrier exists, which is conservative here.
        inner=self.R*1e-7
        if self.minimum_wall_energy(inner,sign)<=K:return self.R
        r=brentq(lambda r:self.minimum_wall_energy(r,sign)-K,inner,self.R,xtol=1e-15,rtol=1e-13)
        return self.R-r

    def orbit(self,K,mu,azimuth,sign=-1,with_time=True):
        p0=sqrt(K*(K+2*MEC2));pp=p0*sqrt(1-mu*mu)*np.cos(azimuth)
        pz=p0*sqrt(1-mu*mu)*np.sin(azimuth);R=self.R
        def state(d):
            r=R-d;local=K-sign*float(self.potential_by_depth_MV(d))
            pphi=(R*pp+sign*self.kB*d*(2*R-d)/2)/r
            f=local*(local+2*MEC2)-pphi*pphi-pz*pz
            if local<0:f=-abs(f)-abs(local)
            return local,pphi,f
        probes=np.r_[0.,np.geomspace(1e-13,R*(1-1e-6),220)]
        last=0.;turn=None
        for d in probes[1:]:
            if state(d)[2]<0:
                turn=brentq(lambda d:state(d)[2],last,d,xtol=2e-17,rtol=1e-13);break
            last=d
        if turn is None:raise RuntimeError('No wall-return turning point in tested interval')
        out={'initial_MeV':K,'inward_direction_cosine':mu,'tangent_azimuth_rad':azimuth,'charge_sign':sign,
             'penetration_m':turn,'minimum_radius_m':R-turn,'turn_radial_momentum2_MeV2':state(turn)[2]}
        if with_time:
            # d=turn*sin(theta)^2 removes the square-root endpoint singularity.
            def integrand(theta,time):
                d=turn*np.sin(theta)**2;local,pphi,f=state(d)
                if f<=0:
                    # Only roundoff at the upper endpoint is allowed. Interior
                    # quadrature nodes are checked by the independent validator.
                    h=max(turn*1e-5,1e-15)
                    deriv=(state(turn)[2]-state(turn-h)[2])/h
                    f=max(-deriv*(turn-d),1e-300)
                factor=2*turn*np.sin(theta)*np.cos(theta)/sqrt(f)
                return factor*((local+MEC2)/C if time else pz)
            out['return_time_s']=2*quad(lambda t:integrand(t,True),0,pi/2,epsabs=1e-19,epsrel=2e-7,limit=150)[0]
            out['axial_displacement_m']=2*quad(lambda t:integrand(t,False),0,pi/2,epsabs=1e-10,epsrel=2e-7,limit=150)[0]
        return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);ap.add_argument('--quick',action='store_true');args=ap.parse_args()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    cylinders={};profiles=[]
    for key in ['22049','23051','24052','25053','25055']:
        c=Cylinder(plant['selected_source_profiles'][key]);cylinders[key]=c;eq=c.eq
        denom=np.interp(c.src['xwall'],eq['x'],eq['ihc']);cdf=eq['ihc']/denom
        quantiles=[]
        for f in [.1,.5,.9,.99,.999]:
            r=float(np.interp(f,cdf,eq['x'])*c.scale)
            quantiles.append({'reaction_fraction_inside':f,'radius_m':r,'potential_MV':float(c.potential_MV(r)),
               'minimum_wall_electron_MeV':c.minimum_wall_energy(r,-1),'minimum_wall_positron_MeV':c.minimum_wall_energy(r,1)})
        profiles.append({'processor':key,'wall_radius_m':c.R,'wall_field_V_m':c.src['wall_E_V_m'],
          'cold_electron_planar_penetration_m':2*ME*c.src['wall_E_V_m']/(E_CHARGE*c.B*c.B),
          'nonrel_electron_gyroperiod_s':2*pi*ME/(E_CHARGE*c.B),
          'perpendicular_radiation_damping_s':3*pi*EPS0*ME**3*C**3/(E_CHARGE**4*c.B*c.B),
          'reaction_quantiles':quantiles})
    c=cylinders['22049'];orbits=[]
    mus=[.1,.5,.9] if args.quick else np.sqrt((np.arange(8)+.5)/8)
    angles=[0.,pi/2,pi,3*pi/2] if args.quick else 2*pi*(np.arange(16)+.5)/16
    energies=[1e-6,1e-3,1.,30.] if args.quick else [1e-6,1e-5,1e-4,1e-3,.01,.1,1.,10.,30.,100.]
    for sign in [-1,1]:
        for K in energies:
            group=[c.orbit(K,float(mu),float(a),sign) for mu in mus for a in angles]
            orbits+=group
            print('sign',sign,'MeV',K,'max depth mm',max(r['penetration_m'] for r in group)*1e3,
                  'max return ps',max(r['return_time_s'] for r in group)*1e12,flush=True)
    bounds=[{'charge_sign':sign,'initial_MeV':K,'maximum_access_depth_m':c.maximum_access_depth(K,sign)} for sign in [-1,1] for K in energies]
    out={'scope':__doc__,'all_direction_access_bounds':bounds,'sampling':'Uniform midpoint samples in mu^2 (Lambert emission) and azimuth; extrema are sampled, not global bounds.',
         'profiles':profiles,'titanium_wall_orbits':orbits,'quick':args.quick,
         'limits':['Prescribed static cylindrical fields; the interior-access energy bound is necessary over all directions.',
          'Finite end regions, collisions, roughness, surface reflection, collective charge and RF fields can change return behavior.',
          'This is a wall-contamination screen; electrons and positrons born inside the plasma need separate axial and radial clearing.']}
    (args.outdir/'electron-wall-orbits.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('Ti quantiles',profiles[0]['reaction_quantiles'])


if __name__=='__main__':main()
