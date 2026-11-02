"""Independent controls and bounded convergence report for axial drive screens."""
from pathlib import Path
import sys,json,argparse
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import i0,wofz
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(HERE))
from axial_basis import matrices,grid,C_LIGHT
from axial_spectral import solve,susceptibility
from axial_capture import minimum_amplitude,initial_phase_fraction,uniform_chi,MC2
from electron_wall_orbits import Cylinder
from nonneutral_screen import MEV,E_CHARGE
from catalyst_penning import MP

def run(kernel,main,independent,refined,out):
    data=np.load(kernel);meta=json.loads(str(data['metadata_json']));c=Cylinder(meta['source']);checks=[]
    def add(name,ok,values=None):checks.append(dict(name=name,pass_=bool(ok),values=values))
    # Independently integrated Gaussian response verifies sign and normalization.
    for zeta in [.1,.5,1.,2.,4.]:
        re=quad(lambda u:np.exp(-u*u)*np.cos(2*zeta*u),0,12,epsabs=1e-12)[0]
        im=quad(lambda u:np.exp(-u*u)*np.sin(2*zeta*u),0,12,epsabs=1e-12)[0]
        numerical=1+2j*zeta*(re+1j*im)
        analytic=1+zeta*1j*np.sqrt(np.pi)*wofz(zeta)
        err=abs(numerical-analytic);add('Gaussian Kubo integral zeta='+str(zeta),err<2e-11,err)
    for k in [10.,40.,80.,160.]:
        nodes=grid(c.R,72,24);K,_,_=matrices(nodes,k);f=np.r_[np.linalg.solve(K[:-1,:-1],-K[:-1,-1]),1.]
        exact=i0(k*nodes)/i0(k*c.R);err=np.max(abs(f-exact))/np.max(abs(exact))
        add('Cylindrical vacuum transfer k='+str(k),err<.0003,float(err))
    # Static susceptibility is local. Compare orbital covariance with direct
    # radial-density quadrature and an independent ODE boundary-value result.
    static=[]
    for nc,no in [(24,8),(36,12),(72,24),(144,48)]:
        nodes=grid(c.R,nc,no);K,_,_=matrices(nodes,40.)
        chi=np.zeros_like(K);xx,ww=np.polynomial.legendre.leggauss(12)
        for sr in meta['rows']:
            dens=c.eq['h'] if sr['Z']==1 else c.eq['c'];integ=c.eq['ih'] if sr['Z']==1 else c.eq['ic']
            norm=np.interp(c.src['xwall'],c.eq['x'],integ);local=np.zeros_like(K)
            for i in range(len(nodes)-1):
                lo,hi=nodes[i:i+2];r=(xx+1)*(hi-lo)/2+lo;x=r/c.scale
                prob=x*np.interp(x,c.eq['x'],dens)/norm/c.scale
                fs=np.stack([(hi-r)/(hi-lo),(r-lo)/(hi-lo)]);weight=ww*(hi-lo)/2*prob
                local[i:i+2,i:i+2]+=np.einsum('aq,bq,q->ab',fs,fs,weight)
            chi-=sr['inventory_per_m']*(sr['Z']*E_CHARGE)**2/(1.6*MEV)*local
        A=K-chi;phi=np.r_[np.linalg.solve(A[:-1,:-1],-A[:-1,-1]),1.]
        x0=c.eq['x'][0];kw=40*c.scale
        def rhs(x,y):return [y[1],(np.interp(x,c.eq['x'],c.eq['h']+22**2*c.eq['c'])+kw*kw)*y[0]-y[1]/x]
        sol=solve_ivp(rhs,[x0,c.src['xwall']],[1.,0.],rtol=2e-10,atol=1e-12,dense_output=True,max_step=.01)
        exact=sol.sol(np.maximum(nodes/c.scale,x0))[0]/sol.y[0,-1]
        medr=np.interp(.5,c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc']),c.eq['x'])*c.scale
        at=np.interp(medr,nodes,phi);ex=float(sol.sol(medr/c.scale)[0]/sol.y[0,-1])
        static.append(dict(ncore=nc,axis_FEM=phi[0],axis_ODE=exact[0],axis_relative_error=phi[0]/exact[0]-1,
                           median_FEM=at,median_ODE=ex,median_relative_error=at/ex-1))
    add('Static screening converges with radial refinement',abs(static[-1]['axis_relative_error'])<.006 and abs(static[-1]['median_relative_error'])<.002,static)
    for sr in meta['rows']:
        Z=sr['Z'];co=data['coef_'+str(Z)];part=co.sum(axis=2);e0=np.max(abs(part[:,0]-1/np.sqrt(meta['n'])));eh=np.max(abs(part[:,1:]))
        add('Partition of unity in radial harmonics Z='+str(Z),e0<1e-12 and eh<1e-12,[float(e0),float(eh)])
        add('Hamiltonian orbit accuracy Z='+str(Z),sr['max_energy_error_over_T']<.001 and sr['wall_crossings']==0,sr)
        for omega in [1e7,1e8,1e9,3e9,1e10]:
            chi=susceptibility(data,sr,40.,omega);summed=chi.sum()
            exact=uniform_chi(omega,40.,1.6e6,sr['mass_over_proton']*MP,Z,sr['inventory_per_m'])
            err=abs(summed-exact)/abs(exact)
            add('Uniform axial susceptibility Z='+str(Z)+' omega='+str(omega),err<2e-8,float(err))
    harmonic=[]
    for k,beta in [(40.,.1),(40.,.2),(40.,.3),(80.,.1),(80.,.2),(80.,.3)]:
        ph64,p64,_,_=solve(data,k,k*C_LIGHT*beta,64);ph32,p32,_,_=solve(data,k,k*C_LIGHT*beta,32)
        er=abs(sum(p32)/sum(p64)-1);ep=np.max(abs(ph32-ph64))/np.max(abs(ph64));harmonic.append(dict(k=k,beta=beta,power_relative_change=er,field_relative_change=float(ep)))
    add('Radial harmonic truncation for reported bands',max(x['power_relative_change'] for x in harmonic)<.01 and max(x['field_relative_change'] for x in harmonic)<.003,harmonic)
    reports={n:json.loads(p.read_text()) for n,p in [('main',main),('independent',independent),('refined',refined)]}
    convergence=[];derived=[]
    for row in reports['main']['rows']:
        k=row['k_per_m'];b=row['beta_phase'];p=np.array(row['species_heat_W_per_m_per_wall_V2'])
        add('Passive species absorption k='+str(k)+' beta='+str(b),min(p)>=-1e-15,p.tolist())
        if b not in [.1,.2,.25,.3]:continue
        r2={n:next(r for r in d['rows'] if r['k_per_m']==k and r['beta_phase']==b) for n,d in reports.items()}
        h=row['total_ion_heat_MeV_per_primary_per_wall_V2'];tr=row['median_reaction_transfer_abs'];phi=minimum_amplitude(b);wall=phi/tr
        rel={n:r['total_ion_heat_MeV_per_primary_per_wall_V2']/h-1 for n,r in r2.items() if n!='main'}
        convergence.append(dict(k=k,beta=b,**rel))
        derived.append(dict(k_per_m=k,beta_phase=b,median_trapping_threshold_V=phi,wall_voltage_for_median_threshold_V=wall,
                            continuous_heat_MeV_per_primary=h*wall*wall,heat_for_half_initial_median_phases_MeV_per_primary=4*h*wall*wall,
                            titanium_qPhi_over_T_at_median=22*phi/1.6e6,seed_relative_change=rel['independent'],mesh_relative_change=rel['refined']))
    # Phase-measure bound independently counted from Hamiltonian inequalities.
    phase=(np.arange(200000)+.5)*2*np.pi/200000
    for beta in [.01,.1,.3,.9]:
        for factor in [1.,2.,5.,40.863458]:
            phi=factor*minimum_amplitude(beta);H=np.ones(len(phase))+phi/MC2*np.cos(phase);sep=np.sqrt(1-beta*beta)+phi/MC2
            measured=np.mean(H<sep);pred=initial_phase_fraction(beta,phi)
            add('Relativistic separatrix phase measure beta='+str(beta)+' factor='+str(factor),abs(measured-pred)<1e-5,[float(measured),pred])
    r=next(r for r in reports['main']['rows'] if r['k_per_m']==40 and r['beta_phase']==.1)
    wall=1000.;heat=sum(r['species_heat_W_per_m_per_wall_V2'])*wall**2;Nsum=sum(x['inventory_per_m'] for x in meta['rows']);U=1.5*Nsum*1.6*MEV
    rate=c.src['reaction_rate_per_m_s'];shot=dict(length_m=10.,burn_s=.001,primary_pn_per_shot=rate*10*.001,
        main_thermal_inventory_J_per_m=U,source_primary_rate_per_m_s=rate,
        eventual_full_network_nuclear_W_per_m_upper_comparison=rate*16.17665*MEV,
        beta_point_one_wall_1kV_heat_W_per_m=heat,one_microsecond_burst_heat_J_per_m=heat*1e-6,
        average_heating_MeV_per_primary_at_one_Hz=heat*1e-6/(rate*MEV),
        median_fixed_bucket_initial_phase_fraction=initial_phase_fraction(.1,1000*r['median_reaction_transfer_abs']),
        scope='The0.1c wave only tests injection/capture; its~2.6keV particles cannot clear40MV end barriers. Whole-plasma dump/reload is a separate method.')
    outdata=dict(scope=__doc__,checks=checks,all_pass=all(x['pass_'] for x in checks),count=len(checks),convergence=convergence,derived_capture_costs=derived,shot_diagnostic=shot,
      uncertainty_note='Weak radial-harmonic absorption has~15–25% independent-ensemble variation at selected points. This is reported explicitly; no percent-level tail precision is claimed.')
    out.write_text(json.dumps(outdata,indent=2,allow_nan=False)+'\n');print('checks',len(checks),'pass',outdata['all_pass']);print(json.dumps(shot,indent=2))
    for ck in checks:
        if not ck['pass_']:print('FAILED',ck)
    if not outdata['all_pass']:raise SystemExit(1)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--kernel',type=Path,required=True);ap.add_argument('--main',type=Path,required=True);ap.add_argument('--independent',type=Path,required=True);ap.add_argument('--refined',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();run(a.kernel,a.main,a.independent,a.refined,a.out)
