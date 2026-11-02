"""Gyro-averaged trace-electron trajectories with diagnostic Landau collisions.

The collision tensor remains unmagnetized at the encounter scale. Gyrophase
is averaged by sampling, while pitch angle, energy and axial position evolve.
Local ion densities and E x B drift are frozen at one selected source radius.
Absorbing ends model electrons only; this is not a positron-removal design.
"""
from pathlib import Path
from math import pi,sqrt,erf,exp,log,cos,sin
import argparse,json,gzip,time
import numpy as np
from numba import njit,prange
from electron_wall_orbits import Cylinder,ME,MP,E_CHARGE,EPS0,C
from lepton_collision_screen import local_state

ROOT=Path(__file__).parent
EV=E_CHARGE


def prepare_species(species,log_lambda=15.):
    return np.array([[4*pi*(z*EV**2/(4*pi*EPS0))**2*n*log_lambda/ME**2,
                      ME/mass,2*T*EV/mass] for mass,n,z,T in species])


@njit(cache=True)
def collision_coeff(speed,table):
    nu=0.;parallel=0.;perp=0.;energy=0.
    a=4/(3*sqrt(pi))
    for j in range(table.shape[0]):
        G,ratio,vth2=table[j]
        x=speed*speed/vth2
        if x<1e-5:
            series=1-3*x/5+3*x*x/14-x*x*x/18
            nu+=G*(1+ratio)*a/vth2**1.5*series
            parallel+=G*a/sqrt(vth2)*series
            perp+=G*a/sqrt(vth2)*(1-x/5+3*x*x/70-x*x*x/126)
            energy+=4*G/sqrt(pi*vth2)*(exp(-x)-2*ratio*x*series/3)
        else:
            psi=erf(sqrt(x))-2/sqrt(pi)*sqrt(x)*exp(-x)
            dpsi=2/sqrt(pi)*sqrt(x)*exp(-x)
            nu+=G*(1+ratio)*psi/speed**3
            parallel+=G*psi/(x*speed)
            perp+=G*((1-1/(2*x))*psi+dpsi)/speed
            energy+=2*G/speed*(dpsi-ratio*psi)
    # energy is d<w^2>/dt, evaluated without cancellation of large angular terms.
    return nu,parallel,perp,energy


@njit(cache=True)
def rotate_direction(nx,ny,nz,q,step_fraction):
    # The first spherical-harmonic moment is exp(-q). Small fixed-angle steps
    # with random azimuth reproduce that moment exactly and converge to the
    # angular diffusion operator as the maximum q per step decreases.
    if q>20.:
        mu=2*np.random.random()-1;theta=2*pi*np.random.random()
        return sqrt(1-mu*mu)*cos(theta),sqrt(1-mu*mu)*sin(theta),mu
    nsub=max(1,int(np.ceil(q/step_fraction)))
    ca=exp(-q/nsub);sa=sqrt(max(0.,1-ca*ca))
    for _ in range(nsub):
        p=sqrt(nx*nx+ny*ny)
        if p>1e-14:tx=-ny/p;ty=nx/p;tz=0.
        else:tx=1.;ty=0.;tz=0.
        sx=ny*tz-nz*ty;sy=nz*tx-nx*tz;sz=nx*ty-ny*tx
        theta=2*pi*np.random.random();ct=cos(theta);st=sin(theta)
        nx=ca*nx+sa*(ct*tx+st*sx)
        ny=ca*ny+sa*(ct*ty+st*sy)
        nz=ca*nz+sa*(ct*tz+st*sz)
    norm=sqrt(nx*nx+ny*ny+nz*nz)
    return nx/norm,ny/norm,nz/norm


@njit(cache=True)
def one_particle(seed,table,U,tau,L,tmax,step_fraction,birth_eV,burn,absorb,B):
    np.random.seed(seed)
    mu=2*np.random.random()-1
    speed=sqrt(2*birth_eV*EV/ME)
    vp=speed*sqrt(1-mu*mu);vz=speed*mu
    initial_vz=vz
    z=L*np.random.random() if absorb else 0.
    t=0.;rad=0.;coll=0.;radial=0.;steps=0;escaped=0
    weighted_E=0.;weighted_perp=0.;weighted_parallel=0.;weight=0.
    maxE=birth_eV;max_step=0.;negative_energy=0;hist=np.zeros(120);pitch=np.zeros(40)
    # Histogram edges are logarithmic in energy, 0.01 eV to 100 keV.
    logmin=log(.01);dlog=log(1e7)/120
    while t<tmax:
        # This bound is independent of the sampled gyrophase. Letting the
        # time step depend on that random phase would bias the phase average.
        wmin=sqrt((vp-abs(U))**2+vz*vz)
        maxnu,_,_,_=collision_coeff(wmin,table)
        dt=min(step_fraction/max(maxnu,1e-100),.01*tau,tmax-t)
        # Resolve the selected burn-in boundary without a partial histogram bin.
        if t<burn and t+dt>burn:dt=burn-t
        max_step=max(max_step,dt)
        oldvp=vp;oldvz=vz;oldE=.5*ME*(vp*vp+vz*vz)/EV
        # Strang splitting of exact perpendicular radiation damping.
        half_damp=exp(-dt/(4*tau))
        vp1=vp*half_damp
        rad1=.5*ME*(vp*vp-vp1*vp1)/EV
        theta=2*pi*np.random.random()
        vx=vp1*cos(theta);vy=vp1*sin(theta)
        wx=vx;wy=vy-U;wz=vz
        w=sqrt(wx*wx+wy*wy+wz*wz)
        nu,dp,dtan,Adot=collision_coeff(w,table)
        invw=1/max(w,1e-100)
        normal=np.random.normal()
        # Freeze the radial diffusion coefficient during the step. This form
        # is positive in the low-speed Bessel limit (Adot -> 3*dp), preserves
        # the first energy moment, and becomes exact energy conservation in
        # the Lorentz limit. Finite-step errors are checked by refinement.
        w2new=(w+sqrt(dp*dt)*normal)**2+(Adot-dp)*dt
        if w2new<0:
            negative_energy+=1;w2new=0.
        magnitude=sqrt(w2new)
        if w>1e-50:
            nx,ny,nz=rotate_direction(wx*invw,wy*invw,wz*invw,
                dtan*dt/(w*w),step_fraction)
        else:
            mu=2*np.random.random()-1;theta=2*pi*np.random.random()
            nx=sqrt(1-mu*mu)*cos(theta);ny=sqrt(1-mu*mu)*sin(theta);nz=mu
        vx2=magnitude*nx;vy2=U+magnitude*ny;vz2=magnitude*nz
        dvy=vy2-vy
        before=.5*ME*(vx*vx+vy*vy+vz*vz)/EV
        after=.5*ME*(vx2*vx2+vy2*vy2+vz2*vz2)/EV
        vp2=sqrt(vx2*vx2+vy2*vy2)
        vp=vp2*half_damp;vz=vz2
        rad2=.5*ME*(vp2*vp2-vp*vp)/EV
        newE=.5*ME*(vp*vp+vz*vz)/EV
        znew=z+.5*(oldvz+vz)*dt
        if not absorb and t<burn and t+dt>=burn:znew=0.
        # The step is short compared with transport times. Linear interpolation
        # at the first end crossing removes leading geometric overshoot.
        fraction=1.
        if absorb and (znew<0 or znew>L):
            boundary=0. if znew<0 else L
            fraction=(boundary-z)/(znew-z)
            escaped=1
        used=dt*fraction
        rad+=(rad1+rad2)*fraction
        coll+=(after-before)*fraction
        radial+=ME*dvy/(-EV*B)*fraction
        # For the final fractional step, interpolate energy components for a
        # pathwise energy ledger. This also avoids pretending the full-step
        # final velocity was measured at the boundary.
        finalE=oldE+(newE-oldE)*fraction
        perp0=.5*ME*oldvp*oldvp/EV
        perp1=.5*ME*vp*vp/EV
        parallel0=.5*ME*oldvz*oldvz/EV
        parallel1=.5*ME*vz*vz/EV
        midE=oldE+.5*fraction*(newE-oldE)
        midperp=perp0+.5*fraction*(perp1-perp0)
        midparallel=parallel0+.5*fraction*(parallel1-parallel0)
        if t>=burn:
            weighted_E+=used*midE;weighted_perp+=used*midperp
            weighted_parallel+=used*midparallel;weight+=used
            k=int((log(max(midE,.01))-logmin)/dlog)
            hist[min(119,max(0,k))]+=used
            # A midpoint estimate of pitch with the same residence-time weight.
            zm=.5*(oldvz+vz);pm=.5*(oldvp+vp)
            cosine=zm/sqrt(max(zm*zm+pm*pm,1e-100))
            j=int((cosine+1)*20)
            pitch[min(39,max(0,j))]+=used
        t+=used;z=z+fraction*(znew-z);maxE=max(maxE,finalE);steps+=1
        if escaped:break
        if steps>20000000:break
    if not escaped:finalE=.5*ME*(vp*vp+vz*vz)/EV
    # Column definitions are recorded by the Python wrapper.
    result=np.array([t,rad,coll,finalE,weighted_E,weighted_perp,weighted_parallel,
                     weight,radial,float(steps),float(escaped),maxE,max_step,z,float(negative_energy),initial_vz,vz])
    return result,hist,pitch


@njit(cache=True,parallel=True)
def ensemble(n,seed,table,U,tau,L,tmax,step_fraction,birth_eV,burn,absorb,B):
    result=np.zeros((n,17));hist=np.zeros((n,120));pitch=np.zeros((n,40))
    for i in prange(n):
        result[i],hist[i],pitch[i]=one_particle((seed+104729*i)%2147483647,
            table,U,tau,L,tmax,step_fraction,birth_eV,burn,absorb,B)
    return result,hist.sum(axis=0),pitch.sum(axis=0)


def summary(result,hist,pitch,birth,tau):
    n=len(result);times=result[:,0];weight=result[:,7].sum()
    def stat(v):return {'mean':float(v.mean()),'standard_error':float(v.std(ddof=1)/sqrt(n))}
    # Sampling errors use independent trajectories, not individual adaptive steps.
    radiative_rate=result[:,5].sum()/weight/tau
    out={'particles':n,'escaped':int(result[:,10].sum()),'censored':int(n-result[:,10].sum()),
        'residence_s':stat(times),'radiated_eV':stat(result[:,1]),
        'collision_energy_eV':stat(result[:,2]),'exit_or_terminal_energy_eV':stat(result[:,3]),
        'diagnostic_radial_displacement_m':stat(result[:,8]),
        'occupation_mean_energy_eV':float(result[:,4].sum()/weight),
        'occupation_Tperp_over_Tparallel':float(result[:,5].sum()/(2*result[:,6].sum())),
        'radiation_eV_per_occupied_second':float(radiative_rate),
        'pathwise_energy_closure_max_eV':float(np.max(np.abs(result[:,3]-birth-result[:,2]+result[:,1]))),
        'largest_energy_eV':float(result[:,11].max()),'steps_total':int(result[:,9].sum()),
        'maximum_step_s':float(result[:,12].max()),
        'negative_radial_energy_proposals':int(result[:,14].sum()),
        'energy_histogram_edges_eV':np.geomspace(.01,1e5,121).tolist(),
        'occupation_probability_per_energy_bin':(hist/weight).tolist(),
        'pitch_bin_edges':np.linspace(-1,1,41).tolist(),
        'occupation_probability_per_pitch_bin':(pitch/weight).tolist(),
        'residence_quantiles_s':{str(f):float(np.quantile(times,f)) for f in [.1,.5,.9,.99,.999]}}
    return out


def run(root,fraction,n,seed,step_fraction,logL=15.,birth=1.,tmax=.1,closed=False,burn=0.):
    plant=next(p for p in json.loads((root/'network-plant.json').read_text())['cases']
        if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049'])
    cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])
    r=float(np.interp(fraction,cdf,c.eq['x'])*c.scale);local=local_state(c,r)
    table=prepare_species(local['species'],logL)
    tau=3*pi*EPS0*ME**3*C**3/(EV**4*c.B**2)
    t=time.monotonic()
    result,hist,pitch=ensemble(n,seed,table,local['relative_azimuthal_m_s'],tau,
        plant['active_length_m'],tmax,step_fraction,birth,burn,not closed,c.B)
    seconds=time.monotonic()-t
    out=summary(result,hist,pitch,birth,tau)
    out.update({'scope':__doc__,'reaction_fraction_inside':fraction,'radius_m':r,
        'log_lambda':logL,'birth_energy_eV':birth,'step_fraction':step_fraction,'seed':seed,
        'time_limit_s':tmax,'burn_in_s':burn,'absorbing_ends':not closed,'runtime_s':seconds,
        'active_length_m':plant['active_length_m'],'local_state':{k:v for k,v in local.items() if k!='species'}})
    out['scheme']='Frozen-coefficient relative-energy update and energy-conserving spherical rotations; exact perpendicular radiation damping; phase-independent adaptive time step.'
    if closed:
        displacements=result[:,13]
        estimates=displacements**2/(2*(tmax-burn))
        out['finite_time_axial_diffusion_m2_s']={'mean':float(estimates.mean()),'standard_error':float(estimates.std(ddof=1)/sqrt(n))}
    return out,result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);ap.add_argument('--name',default='lepton-escape-mc')
    ap.add_argument('--fraction',type=float,default=.5);ap.add_argument('--particles',type=int,default=4096)
    ap.add_argument('--seed',type=int,default=93811);ap.add_argument('--step',type=float,default=.02)
    ap.add_argument('--time',type=float,default=.1);ap.add_argument('--birth',type=float,default=1.)
    ap.add_argument('--logL',type=float,default=15.);ap.add_argument('--closed',action='store_true')
    ap.add_argument('--burn',type=float,default=0.);args=ap.parse_args()
    out,result=run(args.root,args.fraction,args.particles,args.seed,args.step,args.logL,args.birth,args.time,args.closed,args.burn)
    (args.outdir/(args.name+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    # Full independent-trajectory data allow statistical auditing without
    # retaining every microscopic time step.
    columns=['time_s','radiated_eV','collision_eV','final_eV','energy_time_eV_s',
        'perp_energy_time_eV_s','parallel_energy_time_eV_s','occupation_time_s',
        'radial_force_displacement_m','steps','escaped','maximum_eV','maximum_step_s','final_z_m','negative_energy_proposals','initial_vz_m_s','final_vz_m_s']
    payload=json.dumps({'columns':columns,'rows':result.tolist()},separators=(',',':'),allow_nan=False).encode()
    (args.outdir/(args.name+'-trajectories.json.gz')).write_bytes(gzip.compress(payload,mtime=0))
    print(json.dumps({k:v for k,v in out.items() if 'histogram' not in k and 'probability' not in k and 'bin_edges' not in k},indent=2),flush=True)


if __name__=='__main__':main()
