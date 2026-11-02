"""Passive orbit-harmonic Vlasov–Poisson response for axisymmetric axial drive.

Each integrable radial orbit is represented by a full radial period. The radial
basis autocorrelation has positive spectral weights, which are broadened by
free Maxwellian parallel motion. This removes single-time-origin phase noise.
"""
from pathlib import Path
import sys,json,argparse,time
import numpy as np
from scipy.special import wofz
from scipy.integrate import simpson
from numba import njit
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).parent))
from axial_basis import matrices,grid,C_LIGHT
from ion_dipole_spectrum import trajectories,initial_ensemble
from electron_wall_orbits import Cylinder
from catalyst_penning import MP
from nonneutral_screen import MEV,E_CHARGE

@njit
def orbit_periods(rad,dt,nphase):
    npart,nt=rad.shape;periods=np.zeros(npart);samples=np.zeros((npart,nphase));errors=np.zeros(npart);npeaks=np.zeros(npart,np.int64)
    for p in range(npart):
        peaks=np.zeros(nt//2);count=0
        for i in range(1,nt-1):
            if rad[p,i]>rad[p,i-1] and rad[p,i]>=rad[p,i+1]:
                # Parabolic maximum timing is second-order in the sample step.
                den=rad[p,i-1]-2*rad[p,i]+rad[p,i+1]
                frac=.5*(rad[p,i-1]-rad[p,i+1])/den if den!=0 else 0.
                peaks[count]=(i+frac)*dt;count+=1
        npeaks[p]=count
        if count<3:continue
        period=(peaks[count-1]-peaks[0])/(count-1);periods[p]=period
        worst=0.
        for j in range(count-1):worst=max(worst,abs(peaks[j+1]-peaks[j]-period)/period)
        errors[p]=worst
        start=peaks[count//2]
        for j in range(nphase):
            at=(start+period*j/nphase)/dt;i=int(at);f=at-i
            samples[p,j]=(1-f)*rad[p,i]+f*rad[p,i+1]
    return periods,samples,errors,npeaks

def build(root,out,n=1024,seed=42,dt=2e-12,duration=32.768e-9,stride=2,ncore=24,nouter=8,nphase=256,nharm=64,batch=64):
    t0=time.monotonic();plant=next(x for x in json.loads((root/'network-plant.json').read_text())['cases'] if x['assumptions']['B']==140 and x['assumptions']['heat']==.7 and x['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);nodes=grid(c.R,ncore,nouter);om=E_CHARGE*c.B/MP;w=-c.eq['kappa']*c.eq['xi'];nsteps=int(round(duration/dt))
    meta=dict(scope=__doc__,n=n,seed=seed,dt_s=dt,stride=stride,duration_s=nsteps*dt,nphase=nphase,nharm=nharm,ncore=ncore,nouter=nouter,source=c.src,rows=[])
    data=dict(nodes=nodes)
    for Z,M,N in [(1,1.,c.src['proton_inventory_per_m']),(c.eq['Z'],c.eq['M'],c.src['catalyst_inventory_per_m'])]:
        initial=initial_ensemble(c,Z,M,n,seed+Z);coeffs=[];omegas=[];per=[];peaks=[];maxerr=0.;hits=0;parseloss=0.;phaseerr=0.;perioderr=0.
        c0=np.zeros((len(nodes),len(nodes)))
        for lo in range(0,n,batch):
            pos,kin,err,rad,wall=trajectories(initial[lo:lo+batch],Z,M,w,dt*om,nsteps,stride,c.eq['x'][0],c.eq['x'][1]-c.eq['x'][0],c.spline.c,c.src['xwall'])
            radii=np.abs(pos)*c.scale
            period,samples,errs,npeaks=orbit_periods(radii,dt*stride,nphase)
            if np.any(period==0):raise ValueError('Radial orbit had fewer than three maxima')
            perioderr=max(perioderr,float(np.max(errs)));per.extend(period);peaks.extend(npeaks)
            # Hat-basis values sampled uniformly in radial orbital phase.
            f=np.zeros((len(period),nphase,len(nodes)))
            inds=np.clip(np.searchsorted(nodes,samples,side='right')-1,0,len(nodes)-2)
            a=(samples-nodes[inds])/(nodes[inds+1]-nodes[inds]);ii=np.arange(len(period))[:,None];jj=np.arange(nphase)[None,:]
            f[ii,jj,inds]=1-a;f[ii,jj,inds+1]=a
            fft=np.fft.rfft(f,axis=1)/nphase
            # Origin at a radial turning point: all radial basis harmonics are
            # real; retaining tiny sine parts would double phase-space features.
            phaseerr=max(phaseerr,float(np.max(np.abs(fft.imag))))
            full_power=np.sum(np.abs(fft[:,0])**2,axis=1)+2*np.sum(np.abs(fft[:,1:-1])**2,axis=(1,2))+np.sum(abs(fft[:,-1])**2,axis=1)
            keep_power=np.sum(fft[:,0].real**2,axis=1)+2*np.sum(fft[:,1:nharm+1].real**2,axis=(1,2))
            parseloss=max(parseloss,float(np.max(np.maximum(0,full_power-keep_power))))
            coeffs.append(fft[:,:nharm+1].real/np.sqrt(n));omegas.append(2*np.pi/period)
            c0+=np.einsum('pqa,pqb->ab',f,f)/(n*nphase)
            maxerr=max(maxerr,float(np.max(err)));hits+=int(np.sum(wall))
        cf=np.concatenate(coeffs);fr=np.concatenate(omegas)
        data['coef_'+str(Z)]=cf;data['radial_omega_'+str(Z)]=fr;data['c0_time_'+str(Z)]=c0
        meta['rows'].append(dict(Z=Z,mass_over_proton=M,inventory_per_m=N,max_energy_error_over_T=maxerr,wall_crossings=hits,
           max_relative_orbit_period_jitter=perioderr,minimum_sampled_peaks=int(min(peaks)),maximum_sine_coefficient=phaseerr,
           maximum_lost_hat_variance=parseloss,minimum_radial_period_s=float(min(per)),maximum_radial_period_s=float(max(per))))
        print('spectral kernel',Z,'seconds',time.monotonic()-t0,'period jitter',perioderr,'lost variance',parseloss,flush=True)
    meta['seconds']=time.monotonic()-t0;data['metadata_json']=np.array(json.dumps(meta));np.savez_compressed(out,**data)
    print('saved',out,'seconds',time.monotonic()-t0,flush=True)

def susceptibility(data,sr,k,omega,nharm=None):
    Z=sr['Z'];co=data['coef_'+str(Z)];fr=data['radial_omega_'+str(Z)];T=1.6*MEV;vt=np.sqrt(T/(sr['mass_over_proton']*MP));scale=np.sqrt(2)*k*vt
    if nharm is not None:co=co[:,:nharm+1]
    harm=np.arange(co.shape[1]);shift=fr[:,None]*harm[None,:]
    # Both signs of each nonzero radial harmonic, each with its own positive
    # dyad. Kubo's numerator contains the drive omega, not shifted omega.
    def dispersion(x):return 1j*np.sqrt(np.pi)*wofz(x)
    wp=1+omega/scale*dispersion((omega+shift)/scale)
    wm=1+omega/scale*dispersion((omega-shift)/scale)
    weight=wp+wm;weight[:,0]*=.5
    F=co.reshape(-1,co.shape[-1]);ww=weight.reshape(-1)
    return -sr['inventory_per_m']*(Z*E_CHARGE)**2/T*(F.T@(ww[:,None]*F))

def solve(data,k,omega,nharm=None):
    meta=json.loads(str(data['metadata_json']));nodes=data['nodes'];K,_,_=matrices(nodes,k)
    cs=[susceptibility(data,sr,k,omega,nharm) for sr in meta['rows']];chi=sum(cs);A=K-chi
    phi=np.r_[np.linalg.solve(A[:-1,:-1],-A[:-1,-1]),1.+0j]
    power=[float(-omega/2*np.vdot(phi,x@phi).imag) for x in cs]
    return phi,power,cs,np.linalg.cond(A[:-1,:-1])

def responses(kernel,out,kvals=(40.,80.,160.),betas=None,nharm=None):
    start=time.monotonic();data=np.load(kernel);meta=json.loads(str(data['metadata_json']));nodes=data['nodes'];src=meta['source'];c=Cylinder(src)
    if betas is None:betas=np.array([0.,.01,.02,.03,.05,.075,.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.4,.5])
    rad=c.eq['x']*c.scale;rw=c.eq['h']*c.eq['c']*rad;cdf=c.eq['ihc']/np.interp(src['xwall'],c.eq['x'],c.eq['ihc']);med=float(np.interp(.5,cdf,rad))
    rows=[]
    for k in kvals:
        for beta in betas:
            omega=k*C_LIGHT*beta;phi,power,cs,cond=solve(data,k,omega,nharm)
            pf=np.interp(rad,nodes,phi.real)+1j*np.interp(rad,nodes,phi.imag)
            mask=rad<=c.R
            rms=np.sqrt(np.trapezoid((rw*abs(pf)**2)[mask],rad[mask])/np.trapezoid(rw[mask],rad[mask]))
            pm=np.interp(med,nodes,phi.real)+1j*np.interp(med,nodes,phi.imag)
            im_eigen=[float(np.linalg.eigvalsh(-x.imag)[0]) for x in cs]
            rows.append(dict(k_per_m=k,beta_phase=float(beta),omega_rad_s=omega,omega_R_over_c=omega*c.R/C_LIGHT,
                axis_transfer_abs=float(abs(phi[0])),median_reaction_transfer_abs=float(abs(pm)),reaction_rms_transfer=float(rms),
                total_ion_heat_MeV_per_primary_per_wall_V2=sum(power)/(src['reaction_rate_per_m_s']*MEV),species_heat_W_per_m_per_wall_V2=power,
                interior_phi_real=phi.real.tolist(),interior_phi_imag=phi.imag.tolist(),condition_number=float(cond),
                lowest_dissipative_eigenvalues=im_eigen))
        print('spectral response',k,'seconds',time.monotonic()-start,flush=True)
    output=dict(scope=__doc__,kernel_metadata=meta,nodes_m=nodes.tolist(),nharm_used=nharm or meta['nharm'],rows=rows,seconds=time.monotonic()-start,
      formulas={'radial_correlation':'C_ab(t)=sum_l <f_a,l f_b,l*> exp(i*l*Omega_r*t)',
                'susceptibility':'chi_ab=−Nq2/T sum_l <f_a,l f_b,l*> [1+omega/(sqrt(2)k vth) Z((omega+l*Omega_r)/(sqrt(2)k vth))]',
                'passivity':'−Im chi is a sum of positive semidefinite dyads with nonnegative weights for omega>0',
                'field':'(K_Poisson−chi) Phi=0 internally; Phi_wall=1',
                'heat':'P/L=−omega/2 Im[Phi† chi Phi]'},
      limits=['Linear collisionless kinetic-Poisson model; omega is laboratory frequency because m=0.',
              'Radial harmonic convergence, orbit period accuracy and independent invariant sampling must be checked.',
              'Finite ends, electrode layout, nonlinear trapping, trace lepton response and Maxwell propagation remain unresolved.',
              'High beta rows are model diagnostics only; they do not validate a physical near-c accelerating field.'])
    out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['kernel','response']);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--kernel',type=Path);ap.add_argument('--n',type=int,default=1024);ap.add_argument('--seed',type=int,default=42);ap.add_argument('--dt',type=float,default=2e-12);ap.add_argument('--duration',type=float,default=32.768e-9);ap.add_argument('--stride',type=int,default=2);ap.add_argument('--ncore',type=int,default=24);ap.add_argument('--nouter',type=int,default=8);ap.add_argument('--nphase',type=int,default=256);ap.add_argument('--nharm',type=int,default=64);a=ap.parse_args()
    if a.mode=='kernel':build(ROOT,a.out,a.n,a.seed,a.dt,a.duration,a.stride,a.ncore,a.nouter,a.nphase,a.nharm)
    else:responses(a.kernel,a.out,nharm=a.nharm)
