#!/usr/bin/env python3
"""Collisionless 2D-position/3D-velocity moving-blade test particles.
Static channel in rotor frame: 0 < y=x+s*z < gap; 0 < z < L.
Thermal reservoirs at1000K are stationary in lab, flux Maxwellian at axial ends.
Specular energy and direct event histories are checked, no fitted surface data.
"""
import numpy as np
import json, math, argparse, hashlib,time
from pathlib import Path
KB=1.380649e-23
MU=1.66053906660e-27
MASS={'He':4.00260325413*MU,'H2':2*1.00782503223*MU,'H':1.00782503223*MU}

def wilson(k,n,z=1.959963984540054):
    p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0.,mid-h),min(1.,mid+h)]

def trace(n=100000,gas='He',u=700.,theta=45.,L_pitch=1.,open_fraction=.65,alpha=0.,direction=1,seed=1,T=1000.,accelerate=True,max_events=10000):
    rng=np.random.default_rng(seed);m=MASS[gas];sig=math.sqrt(KB*T/m);s=math.tan(math.radians(theta));norm=math.hypot(1,s);gap=open_fraction;L=L_pitch
    y=rng.uniform(0,gap,n);z=np.full(n,0. if direction==1 else L)
    vx=rng.normal(-u,sig,n);vz=direction*sig*np.sqrt(-2*np.log(rng.uniform(0,1,n)));vy=rng.normal(0,sig,n)
    vx0=vx.copy();e0=.5*m*(vx*vx+vz*vz+vy*vy);vxl0=vx+u;el0=.5*m*(vxl0*vxl0+vz*vz+vy*vy)
    status=np.zeros(n,dtype=np.int8);bounces=np.zeros(n,dtype=np.int64);flight=np.zeros(n);specrel=0.;roundoff=0.;iterations=0;skipped=0
    for iteration in range(max_events):
        ids=np.flatnonzero(status==0)
        if not len(ids):break
        iterations=iteration+1
        dy=vx[ids]+s*vz[ids]
        tw=np.where(dy>0,(gap-y[ids])/dy,-y[ids]/dy)
        te=np.where(vz[ids]>0,(L-z[ids])/vz[ids],-z[ids]/vz[ids])
        roundoff=max(roundoff,float(max(0.,-np.min(tw))))
        # Roundoff only: a reflected wall starts exactly on its inward boundary.
        tw=np.maximum(tw,0.);te=np.maximum(te,0.)
        exits=te<=tw;ie=ids[exits];iw=ids[~exits]
        status[ie]=np.where(vz[ie]>0,1,-1);flight[ie]+=te[exits]
        z[ie]+=vz[ie]*te[exits];y[ie]+=dy[exits]*te[exits]
        if not len(iw):continue
        dt=tw[~exits];upper=dy[~exits]>0
        z[iw]+=vz[iw]*dt;flight[iw]+=dt;y[iw]=np.where(upper,gap,0.);bounces[iw]+=1
        oldvx=vx[iw].copy();oldvz=vz[iw].copy();oldvy=vy[iw].copy()
        diffuse=rng.random(len(iw))<alpha if alpha else np.zeros(len(iw),dtype=bool)
        isp=iw[~diffuse]
        if len(isp):
            dot=(vx[isp]+s*vz[isp])/(norm*norm)
            before=vx[isp]**2+vz[isp]**2+vy[isp]**2
            vx[isp]-=2*dot;vz[isp]-=2*s*dot
            after=vx[isp]**2+vz[isp]**2+vy[isp]**2
            specrel=max(specrel,float(np.max(np.abs(after-before)/before)))
        idf=iw[diffuse]
        if len(idf):
            inward=np.where(upper[diffuse],-1.,1.)
            vn=sig*np.sqrt(-2*np.log(rng.uniform(0,1,len(idf))));vt=rng.normal(0,sig,len(idf))
            vx[idf]=(inward*vn-s*vt)/norm;vz[idf]=(inward*s*vn+vt)/norm;vy[idf]=rng.normal(0,sig,len(idf))
        # Exact two-reflection skip for a wholly specular corridor. After a full
        # normal zigzag velocity repeats and axial position changes linearly.
        # Only cycles with all vertex z positions strictly inside the ends skip.
        if alpha==0 and accelerate:
            dy2=vx[iw]+s*vz[iw];twall=gap/np.abs(dy2)
            z1=z[iw]+vz[iw]*twall
            dot=dy2/(norm*norm);vz1=vz[iw]-2*s*dot
            dzcycle=(vz[iw]+vz1)*twall
            z2=z[iw]+dzcycle
            hi=np.maximum(np.maximum(z[iw],z1),z2);lo=np.minimum(np.minimum(z[iw],z1),z2)
            safe=(hi<L)&(lo>0)&(dzcycle!=0)
            j=np.flatnonzero(safe)
            if len(j):
                d=dzcycle[j]
                # Conservative floor: leave at least the final potential exit
                # cycle to explicit event detection; exact at all retained points.
                room=np.where(d>0,L-hi[j],lo[j]);cycles=np.floor(room/np.abs(d)).astype(np.int64)
                good=cycles>0;j=j[good];cycles=cycles[good];q=iw[j]
                z[q]+=cycles*dzcycle[j];flight[q]+=cycles*2*twall[j];bounces[q]+=2*cycles;skipped+=int(np.sum(cycles))
    resolved=status!=0;trans=status==direction;k=int(np.sum(trans));unresolved=int(np.sum(~resolved));nret=int(np.sum(status==-direction))
    ef=.5*m*(vx*vx+vz*vz+vy*vy);elf=.5*m*((vx+u)**2+vz*vz+vy*vy);work=m*u*(vx-vx0)
    closure=np.max(np.abs((elf-el0)-(ef-e0)-work)/(KB*T+np.abs(elf)+np.abs(el0)))
    return {'outcome_sha256':hashlib.sha256(status.tobytes()).hexdigest(),'reflection_count_sha256':hashlib.sha256(bounces.tobytes()).hexdigest(),'n':n,'transmitted':k,'returned':nret,'unresolved':unresolved,'T':k/n,'T_wilson95':wilson(k,n),'T_wilson80':wilson(k,n,1.2815515655446004),'mean_reflections':float(np.mean(bounces)),'maximum_reflections':int(np.max(bounces)),'mean_flight_pitch_over_speed_s':float(np.mean(flight)),'explicit_iterations':iterations,'skipped_specular_cycles':skipped,'max_specular_relative_energy_error':specrel,'max_work_heat_energy_closure_relative_error':float(closure),'max_negative_collision_time_roundoff_s':roundoff,'mean_wall_mechanical_work_per_incident_J':float(np.mean(work)),'mean_wall_heat_to_gas_per_incident_J':float(np.mean(ef-e0)),'mean_total_lab_energy_change_per_incident_J':float(np.mean(elf-el0)),'seed':seed}

def case(params,n,seed):
    f=trace(n=n,direction=1,seed=seed,**params);r=trace(n=n,direction=-1,seed=seed+1,**params)
    tf=f['T'];tr=r['T'];o=params['open_fraction'];df=tf-tr
    kval=tf/tr if tr else None
    kinterval=[f['T_wilson95'][0]/r['T_wilson95'][1],f['T_wilson95'][1]/r['T_wilson95'][0]] if r['T_wilson95'][0]>0 else [f['T_wilson95'][0]/r['T_wilson95'][1],None]
    se=math.sqrt(tf*(1-tf)/n+tr*(1-tr)/n)
    return {'parameters':params,'forward':f,'reverse':r,'K0':kval,'K0_conservative_interval_from_marginal95':kinterval,'net_long_cascade_capture':o*df,'net_capture_normal95':[o*(df-1.95996398454*se),o*(df+1.95996398454*se)],'single_row_zero_backpressure_face_capture':o*tf,'conditional_identical_thermalized_stages_for1e8':math.ceil(math.log(1e8)/math.log(kval)) if kval and kval>1 else None,'conditional_identical_thermalized_stages_for1e9':math.ceil(math.log(1e9)/math.log(kval)) if kval and kval>1 else None}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--n',type=int,default=100000);ap.add_argument('--mode',default='grid',choices=['grid','controls','refine']);args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    params=[]
    if args.mode=='grid':
        for gas in ['He','H2']:
            for theta in [20,30,45,60]:
                for L in [1.,2.,3.]:
                    for alpha in [0.,.5,1.]:params.append(dict(gas=gas,u=700.,theta=theta,L_pitch=L,open_fraction=.65,alpha=alpha,T=1000.))
    elif args.mode=='controls':
        for gas in ['He','H2']:
            for theta in [0,45]:
                for alpha in [0.,.5,1.]:params.append(dict(gas=gas,u=0.,theta=theta,L_pitch=2.,open_fraction=.65,alpha=alpha,T=1000.))
        for gas in ['He','H2']:
            params.append(dict(gas=gas,u=700.,theta=0,L_pitch=2.,open_fraction=.65,alpha=0.,T=1000.))
    else:
        for gas in ['He','H2']:
            for u in [600.,750.]:
                for theta in [30,45]:
                    for o in [.5,.8]:
                        for alpha in [0.,1.]:params.append(dict(gas=gas,u=u,theta=theta,L_pitch=1.,open_fraction=o,alpha=alpha,T=1000.))
    results=[];t=time.monotonic()
    for i,p in enumerate(params):
        rr=case(p,args.n,1000000+100000*['grid','controls','refine'].index(args.mode)+2*i);results.append(rr)
        (out/(args.mode+'.json')).write_text(json.dumps({'scope':__doc__,'mode':args.mode,'n_per_direction':args.n,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':results},indent=2)+'\n')
        print(json.dumps({'done':i+1,'of':len(params),'gas':p['gas'],'theta':p['theta'],'u':p['u'],'alpha':p['alpha'],'L_pitch':p['L_pitch'],'open':p['open_fraction'],'K':rr['K0'],'capture':rr['net_long_cascade_capture'],'unresolved':rr['forward']['unresolved']+rr['reverse']['unresolved'],'seconds':round(time.monotonic()-t,1)}),flush=True)
if __name__=='__main__':main()
