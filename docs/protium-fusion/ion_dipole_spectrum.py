"""Frozen-equilibrium ion dipole spectra for a candidate rotating clearing field.

Thermal orbits use the actual nonlinear radial field and the exact rotating
frame. The spectral response omits induced plasma fields, collisions and ends.
It is a bare-orbit response kernel, not the absorbed power of a self-consistent
RF mode or a net electrical load. Both heat and torque are retained.
"""
from pathlib import Path
import argparse,json,time,gzip
import numpy as np
from scipy.stats import qmc
from scipy.special import ndtri
from scipy.fft import fft,fftfreq
from numba import njit,prange
from electron_wall_orbits import Cylinder
from nonneutral_screen import MEV,E_CHARGE
from catalyst_penning import MP
ROOT=Path(__file__).parent

@njit
def psider(r,grid0,step,co):
    idx=min(max(int((r-grid0)/step),0),co.shape[1]-1)
    y=r-(grid0+idx*step)
    return ((co[0,idx]*y+co[1,idx])*y+co[2,idx])*y+co[3,idx],(3*co[0,idx]*y+2*co[1,idx])*y+co[2,idx]

@njit
def acceleration(x,y,z,M,w,grid0,step,co):
    r=np.sqrt(x*x+y*y)
    if r<grid0:
        _,d=psider(grid0,grid0,step,co);factor=-z/M*d/grid0+w*(z/M+w)
    else:
        _,d=psider(r,grid0,step,co);factor=-z/M*d/r+w*(z/M+w)
    return factor*x,factor*y

@njit
def energy(x,y,u,v,z,M,w,grid0,step,co):
    r=np.sqrt(x*x+y*y);psi,_=psider(max(r,grid0),grid0,step,co)
    return .5*M*(u*u+v*v)+z*psi-.5*M*w*(z/M+w)*r*r

@njit(parallel=True)
def trajectories(initial,z,M,w,h,nsteps,stride,grid0,step,co,rwall):
    npart=initial.shape[0];nt=nsteps//stride
    pos=np.empty((npart,nt),np.complex128);kin=np.empty((npart,nt),np.float64)
    error=np.zeros(npart);maxr=np.zeros(npart);wall=np.zeros(npart,np.int64)
    om=z/M+2*w;angle=om*h;ca=np.cos(angle);sa=np.sin(angle)
    if abs(om)>1e-14:aa=np.sin(angle)/om;bb=(1-np.cos(angle))/om
    else:aa=h;bb=0.
    for p in prange(npart):
        x,y,u,v=initial[p];H0=energy(x,y,u,v,z,M,w,grid0,step,co)
        for n in range(nsteps):
            if n%stride==0:
                j=n//stride;pos[p,j]=x+1j*y;kin[p,j]=u*u+v*v
                error[p]=max(error[p],abs(energy(x,y,u,v,z,M,w,grid0,step,co)-H0))
                maxr[p]=max(maxr[p],np.sqrt(x*x+y*y))
                if maxr[p]>rwall:wall[p]=1
            fx,fy=acceleration(x,y,z,M,w,grid0,step,co);u+=h*.5*fx;v+=h*.5*fy
            x+=aa*u+bb*v;y+=aa*v-bb*u
            u,v=ca*u+sa*v,ca*v-sa*u
            fx,fy=acceleration(x,y,z,M,w,grid0,step,co);u+=h*.5*fx;v+=h*.5*fy
    return pos,kin,error,maxr,wall

def initial_ensemble(c,z,M,n,seed):
    assert n&(n-1)==0
    samples=qmc.Sobol(4,scramble=True,seed=seed).random_base2(int(np.log2(n)))
    samples=np.clip(samples,np.finfo(float).eps,1-np.finfo(float).eps)
    key='ih' if z==1 else 'ic';cdf=c.eq[key]/np.interp(c.src['xwall'],c.eq['x'],c.eq[key])
    r=np.interp(samples[:,0],cdf,c.eq['x']);theta=2*np.pi*samples[:,1]
    return np.c_[r*np.cos(theta),r*np.sin(theta),ndtri(samples[:,2])/np.sqrt(M),ndtri(samples[:,3])/np.sqrt(M)]

def run(root,out,n=2048,seed=42,dt=2e-12,duration=131.072e-9,stride=8,batch=128):
    start=time.monotonic();plant=next(x for x in json.loads((root/'network-plant.json').read_text())['cases'] if x['assumptions']['B']==140 and x['assumptions']['heat']==.7 and x['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);Omega=E_CHARGE*c.B/MP;w=-c.eq['kappa']*c.eq['xi'];TJ=c.T*MEV
    nsteps=int(round(duration/dt));assert nsteps%stride==0;nt=nsteps//stride;sample_dt=dt*stride;duration=nsteps*dt
    win=np.hanning(nt);wnorm=np.sum(win**2)*sample_dt
    freq=2*np.pi*fftfreq(nt,sample_dt);keep=(np.abs(freq)<1.8e10)
    rows=[]
    for z,M,inventory in [(1,1.,c.src['proton_inventory_per_m']),(c.eq['Z'],c.eq['M'],c.src['catalyst_inventory_per_m'])]:
        init=initial_ensemble(c,z,M,n,seed+z);psd=np.zeros(nt);r2=np.zeros(nt);u2=np.zeros(nt);maxerr=0.;maxr=0.;hits=0
        for lo in range(0,n,batch):
            state=init[lo:lo+batch]
            pos,kin,err,rad,wall=trajectories(state,z,M,w,dt*Omega,nsteps,stride,c.eq['x'][0],c.eq['x'][1]-c.eq['x'][0],c.spline.c,c.src['xwall'])
            pos*=c.scale
            ft=fft(pos*win,axis=1)*sample_dt
            psd+=np.sum(np.abs(ft)**2,axis=0)/wnorm/n
            r2+=np.sum(np.abs(pos)**2,axis=0)/n
            u2+=np.sum(kin,axis=0)/n
            maxerr=max(maxerr,float(np.max(err)));maxr=max(maxr,float(np.max(rad))*c.scale);hits+=int(np.sum(wall))
        # Positive frequency means a field rotating positively in the common
        # ion frame. Its laboratory frequency is frequency + w*Omega.
        power=inventory*(z*E_CHARGE)**2*freq**2*psd/(4*TJ)
        vscale=c.scale*Omega
        dwin=np.pi/((nt-1)*sample_dt)*np.sin(2*np.pi*np.arange(nt)/(nt-1))
        derivative_power=float(np.sum(win**2*u2*vscale**2+dwin**2*r2+win*dwin*np.gradient(r2,sample_dt))/np.sum(win**2))
        rows.append({'Z':z,'mass_over_proton':M,'inventory_per_m':inventory,
            'rotating_frequency_rad_s':freq[keep].tolist(),'dipole_PSD_m2_s':psd[keep].tolist(),
            'bare_rotating_heat_W_per_m_per_E2':power[keep].tolist(),
            'bare_rotating_heat_MeV_per_primary_per_E2':(power[keep]/(c.src['reaction_rate_per_m_s']*MEV)).tolist(),
            'initial_mean_r2_m2':float(r2[0]),'time_mean_r2_m2':float(np.mean(r2)),
            'weighted_time_mean_r2_m2':float(np.sum(r2*win**2)/np.sum(win**2)),
            'full_PSD_variance_m2':float(np.sum(psd)/duration),
            'full_PSD_velocity_m2_s2':float(np.sum(freq**2*psd)/duration),
            'window_derivative_identity_m2_s2':derivative_power,
            'r2_time_cv':float(np.std(r2)/np.mean(r2)),
            'initial_mean_u2_scaled':float(u2[0]),'time_mean_u2_scaled':float(np.mean(u2)),
            'max_rotating_energy_error_over_T':maxerr,'maximum_sampled_radius_m':maxr,'sampled_wall_crossings':hits})
        print('Z',z,'max energy error/T',maxerr,'r2 CV',rows[-1]['r2_time_cv'],'wall crossings',hits,flush=True)
    d={'scope':__doc__,'n':n,'seed':seed,'dt_s':dt,'duration_s':duration,'stride':stride,'window':'Hann; two-sided PSD normalized by integral window squared',
       'source':c.src,'ion_rotation_rad_s':w*Omega,'frequency_resolution_rad_s':2*np.pi/duration,'rows':rows,'seconds':time.monotonic()-start,
       'linear_response_definition':'P_rot/(length E0^2)=N_per_length q^2 delta^2 S_z(delta)/(4 kT), for a circular field E0(cos delta*t,sin delta*t).',
       'work_and_torque':'tau=P_rot/delta; P_lab=(delta+w)*tau. Signed w is negative. Self-consistent induced fields and maintaining the source are outside this bare-orbit kernel.'}
    out.write_bytes(gzip.compress(json.dumps(d,separators=(',',':'),allow_nan=False).encode(),mtime=0));print('seconds',d['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=ROOT/'ion-dipole-spectrum.json.gz');ap.add_argument('--n',type=int,default=2048);ap.add_argument('--seed',type=int,default=42);ap.add_argument('--dt',type=float,default=2e-12);ap.add_argument('--duration',type=float,default=131.072e-9);ap.add_argument('--stride',type=int,default=8);ap.add_argument('--batch',type=int,default=128);a=ap.parse_args();run(a.root,a.out,a.n,a.seed,a.dt,a.duration,a.stride,a.batch)
