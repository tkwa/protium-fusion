"""Spherical Lagrangian expansion with a local nuclear reaction network.

Cell masses and isotope abundances are Lagrangian. Radii and velocities are
at shell boundaries; pressure and temperature are cell quantities. The same
pressure (including artificial viscosity) enters force and PdV work. Thus
kinetic/internal energy exchange cancels in the semidiscrete equations.
BDF time integration has finite tolerance error, measured in every result.

This begins with compressed fuel at rest. It is not an implosion model.
Self-gravity is omitted at the small pulse masses examined here.
Optional grey flux-limited radiation diffusion uses a constant opacity.
"""
import argparse
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix,coo_matrix
from scipy.optimize._numdiff import group_columns
import pulse_model as pm
from fast_eos import FastEOS,KB_CGS,NA,ARAD_CGS

NN=pm.NN;NR=pm.NR;NV=NN+7
STOICH=pm.STOICH;COEFF=pm.COEFF;RID=pm.RATE_INDEX
PREF=pm.PREF;DEXP=pm.DEXP;REACT=pm.REACT
Q=pm.Q;QNU=pm.QNU;ZDEL=pm.ZZ@pm.STOICH
MEV_PER_U_CGS=pm.MEV/pm.MU*1e4
TI=NN;RI=NN+1;VI=NN+2;NU_I=NN+3;RAD_I=NN+4;HEAT_I=NN+5;CHARGE_I=NN+6


@njit(cache=True)
def network(temp,rho,y):
    nz=len(temp)
    dy=np.zeros((nz,NN));qtot=np.zeros(nz);qnu=np.zeros(nz);qcharge=np.zeros(nz)
    for i in range(nz):
        t=temp[i]/1e9
        tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
        rates=np.zeros(NR)
        for s in range(len(COEFF)):
            exponent=0.
            for k in range(7): exponent+=COEFF[s,k]*tf[k]
            rates[RID[s]]+=np.exp(exponent)
        for j in range(NR):
            f=rates[j]*PREF[j]*rho[i]**DEXP[j]
            for k in range(3):
                n=REACT[j,k]
                # The mass-action polynomial remains differentiable during
                # Newton iterations. Clipping trace abundances made the stiff
                # solve stall at zero; physical solutions are nonnegative.
                if n>=0: f*=y[i,n]
            for n in range(NN): dy[i,n]+=STOICH[n,j]*f
            qtot[i]+=Q[j]*f*MEV_PER_U_CGS
            qnu[i]+=QNU[j]*f*MEV_PER_U_CGS
            qcharge[i]+=ZDEL[j]*f
    return dy,qtot,qnu,qcharge


class IdealEOS:
    """Verification-only gamma=5/3 ideal ions and electrons."""
    def state(self,rho,T,ye,yi):
        p=rho*KB_CGS*NA*T*(ye+yi);e=1.5*p/rho
        return SimpleNamespace(p=p,e=e,dp_drho=p/rho,dp_dT=p/T,
          de_drho=np.zeros_like(e),de_dT=e/T,c_v=e/T,
          gamma1=np.full_like(e,5/3))


def initial(stage):
    if stage=='He':
        y=np.zeros(NN);y[pm.SPECIES.index('he4')]=1/4
        return y,'c12',1/12
    cat,protons,dest={'A':('n14',1,'o15'),'B':('n15',3,'o14'),
      'B1':('n15',1,'c12'),'B2':('c12',2,'o14'),
      'C13':('c13',1,'n14')}[stage]
    idx=pm.SPECIES.index(cat);denom=pm.AA[idx]+protons
    y=np.zeros(NN);y[idx]=1/denom;y[pm.SPECIES.index('p')]=protons/denom
    return y,dest,1/denom


def simulate(stage='A',temp0=3e8,rho0_si=1e8,radius0=100.,*,zones=32,
             nuclear=True,ideal_test=False,radiation=False,kappa=.001,
             final_time_units=2.,rtol=1e-6,art_linear=.5,art_quad=2.,
             save_trace=None,reduced_a=False,initial_profile=None,atol=1e-11,
             carbon13_fraction=0.,viscosity_mode="volume_gate", trace_samples=240):
    eos=IdealEOS() if ideal_test else FastEOS()
    y0,target,cats0=initial(stage)
    if carbon13_fraction:
        if stage!='B2' or not 0<=carbon13_fraction<=1:
            raise ValueError('Carbon isotope mixture requires B2 and a fraction in [0,1]')
        y0=np.zeros(pm.NN);denom=14+carbon13_fraction
        y0[pm.SPECIES.index('p')]=2/denom
        y0[pm.SPECIES.index('c12')]=(1-carbon13_fraction)/denom
        y0[pm.SPECIES.index('c13')]=carbon13_fraction/denom
        cats0=1/denom
    species=pm.SPECIES;AA=pm.AA;ZZ=pm.ZZ;MASS=pm.MASS
    NN=len(species);NV=NN+7
    if reduced_a:
        if stage!='A': raise ValueError('Reduced network is only for stage A')
        ids=[pm.SPECIES.index(n) for n in ['p','n14','o15']]
        species=['p','n14','o15'];y0=y0[ids]
        AA=AA[ids];ZZ=ZZ[ids];MASS=MASS[ids];NN=3;NV=NN+7
        forward=pm.LIB.get_rate_by_name('n14(p,g)o15')
        reverse=pm.LIB.get_rate_by_name('o15(,p)n14')
        af=np.array([s.a for s in forward.sets]);ar=np.array([s.a for s in reverse.sets])
        q_capture=MASS[0]+MASS[1]-MASS[2]
    TI=NN;RI=NN+1;VI=NN+2;NU_I=NN+3;RAD_I=NN+4;HEAT_I=NN+5;CHARGE_I=NN+6
    rho0=rho0_si/1000;R0=radius0*100
    x0=np.linspace(0,1,zones+1)
    mass_total=4*np.pi/3*rho0*R0**3
    mass=mass_total*np.diff(x0**3)
    node_mass=np.concatenate([.5*(mass[:-1]+mass[1:]),[.5*mass[-1]]])
    ye0=float(ZZ@y0);yi0=float(np.sum(y0))
    center=eos.state(rho0,temp0,ye0,yi0)
    cs0=float(np.sqrt(center.gamma1*center.p/rho0));tunit=R0/cs0
    uunit=float(center.e)
    temp_profile=np.ones(zones)
    if ideal_test:
        r2mean=.6*np.diff(x0**5)/np.diff(x0**3)
        temp_profile=1-r2mean
    z0=np.zeros((zones,NV));z0[:,:NN]=y0
    z0[:,TI]=temp_profile;z0[:,RI]=x0[1:]
    if initial_profile is not None:
        z0[:,TI]=np.asarray(initial_profile['T_K'])/temp0
        z0[:,RI]=np.asarray(initial_profile['outer_radius_m'])*100/R0
        z0[:,VI]=np.asarray(initial_profile['velocity_m_s'])*100/cs0
    init_r=np.concatenate([[0.],z0[:,RI]*R0])
    init_rho=mass/(4*np.pi/3*np.diff(init_r**3))
    initstate=eos.state(init_rho,temp0*z0[:,TI],
                        np.full(zones,ye0),np.full(zones,yi0))
    initial_internal=float(mass@initstate.e)
    initial_kinetic=float(.5*node_mass@((z0[:,VI]*cs0)**2))
    initial_energy=initial_internal+initial_kinetic
    # Cell density depends on two bounding nodes. Radiation's inner face
    # therefore reaches two blocks left; force/work alone are tridiagonal.
    sparsity=lil_matrix((zones*NV,zones*NV),dtype=int)
    for i in range(zones):
        for j in range(max(0,i-(2 if radiation else 1)),min(zones,i+2)):
            sparsity[i*NV:(i+1)*NV,j*NV:(j+1)*NV]=1
    sparsity=sparsity.tocsr()
    jac_rows,jac_cols=sparsity.nonzero()
    column_groups=group_columns(sparsity)

    def physical(z):
        yy=z[:,:NN]
        rr=np.concatenate([[0.],z[:,RI]*R0])
        vv=np.concatenate([[0.],z[:,VI]*cs0])
        widths=np.diff(rr)
        if np.any(widths<=0): raise ValueError('Shell crossing')
        volumes=4*np.pi/3*np.diff(rr**3)
        density=mass/volumes
        temp=z[:,TI]*temp0
        yclip=np.maximum(yy,0)
        total=yclip@AA
        ye=(yclip@ZZ)/total;yi=yclip.sum(axis=1)/total
        state=eos.state(density,temp,ye,yi)
        return yy,rr,vv,density,temp,ye,yi,state

    def rhs(t,zflat):
        z=zflat.reshape((zones,NV))
        yy,rr,vv,rho,temp,ye,yi,state=physical(z)
        if nuclear and reduced_a:
            tt=temp/1e9
            tf=np.array([np.ones(zones),1/tt,tt**(-1/3),tt**(1/3),tt,tt**(5/3),np.log(tt)])
            rf=np.exp(af@tf).sum(axis=0);rrate=np.exp(ar@tf).sum(axis=0)
            yp=np.maximum(yy,0)
            flow=rf*rho*yp[:,0]*yp[:,1]-rrate*yp[:,2]
            dy=np.column_stack([-flow,-flow,flow])
            qtot=q_capture*flow*MEV_PER_U_CGS;qnu=np.zeros(zones);qcharge=np.zeros(zones)
        elif nuclear: dy,qtot,qnu,qcharge=network(temp,rho,yy)
        else: dy=np.zeros_like(yy);qtot=np.zeros(zones);qnu=qtot.copy();qcharge=qtot.copy()
        # Isotropic artificial pressure dissipates only in volume compression.
        # A negative radial velocity difference alone does not imply negative
        # spherical divergence. Applying q in expansion violates entropy.
        dvdt=4*np.pi*np.diff(rr**2*vv)
        compress=np.maximum(-np.diff(vv),0.)
        if viscosity_mode == "volume_gate":
            compress=np.where(dvdt<0,compress,0.)
        elif viscosity_mode == "compression_cap":
            compress=np.minimum(compress,np.maximum(-dvdt*rho/mass*np.diff(rr),0.))
        elif viscosity_mode == "divergence":
            # Calibrate to the old radial gradient for homologous flow.
            compress=np.maximum(-dvdt*rho/mass*np.diff(rr)/3.,0.)
        elif viscosity_mode != "legacy":
            raise ValueError("Unknown viscosity mode")
        cs=np.sqrt(np.maximum(state.gamma1*state.p/rho,0.))
        qvis=rho*(art_quad*compress**2+art_linear*cs*compress)
        peff=state.p+qvis
        force=4*np.pi*rr[1:]**2*(peff-np.concatenate([peff[1:],[0.]]))
        accel=force/node_mass
        dvdt=4*np.pi*np.diff(rr**2*vv)
        rhodot=-rho**2*dvdt/mass
        qrad=np.zeros(zones)
        if radiation:
            # cgs opacity; 1 m2/kg = 10 cm2/g.
            opacity=kappa*10
            er=ARAD_CGS*temp**4
            rc=.5*(rr[1:]+rr[:-1])
            rho_face=.5*(rho[1:]+rho[:-1])
            e_face=.5*(er[1:]+er[:-1])
            grad=np.diff(er)/np.diff(rc)
            ratio=np.abs(grad)/(opacity*rho_face*np.maximum(e_face,1e-300))
            limiter=(2+ratio)/(6+3*ratio+ratio**2)
            flux_inner=-pm.C*100*limiter*grad/(opacity*rho_face)
            flux_outer=pm.C*100*er[-1]/(4+3*opacity*rho[-1]*(rr[-1]-rc[-1]))
            lum=np.concatenate([[0.],4*np.pi*rr[1:-1]**2*flux_inner,
                                [4*np.pi*rr[-1]**2*flux_outer]])
            qrad=np.diff(lum)/mass
        de_dyi=1.5*KB_CGS*NA*temp
        uion=de_dyi*yi
        de_dye=(state.e+rho*state.de_drho-uion)/ye
        tdot=(qtot-qnu-qrad+(peff/rho**2-state.de_drho)*rhodot
              -de_dyi*np.sum(dy,axis=1)-de_dye*(dy@ZZ))/state.c_v
        out=np.zeros_like(z)
        out[:,:NN]=dy*tunit;out[:,TI]=tdot*tunit/temp0
        out[:,RI]=vv[1:]*tunit/R0;out[:,VI]=accel*tunit/cs0
        out[:,NU_I]=qnu*tunit/uunit;out[:,RAD_I]=qrad*tunit/uunit
        out[:,HEAT_I]=qtot*tunit/uunit;out[:,CHARGE_I]=qcharge*tunit
        return out.ravel()

    def jacobian(t,zflat):
        # Fixed-scale grouped differences avoid adaptive-step overflow on
        # exactly absent isotopes, whose Jacobian columns can be zero.
        zz=zflat.reshape((zones,NV))
        hh=1e-7*np.maximum(np.abs(zz),1e-2)
        hh[:,:NN]=1e-7*np.maximum(np.abs(zz[:,:NN]),1e-8)
        hh=hh.ravel()
        f0=rhs(t,zflat)
        entries=np.empty(len(jac_rows))
        for group in np.unique(column_groups):
            active=(column_groups==group)
            ff=rhs(t,zflat+hh*active)-f0
            select=(column_groups[jac_cols]==group)
            entries[select]=ff[jac_rows[select]]/hh[jac_cols[select]]
        return coo_matrix((entries,(jac_rows,jac_cols)),shape=sparsity.shape).tocsc()

    def hot(t,z): return np.max(z.reshape((zones,NV))[:,TI])*temp0-2.8e9
    hot.terminal=True;hot.direction=1
    begin=time.monotonic()
    sol=solve_ivp(rhs,(0,final_time_units),z0.ravel(),method='BDF',
      rtol=rtol,atol=atol,jac=jacobian,max_step=.3/zones,events=hot)
    z=sol.y[:,-1].reshape((zones,NV))
    yy,rr,vv,rho,temp,ye,yi,state=physical(z)
    kinetic=float(.5*node_mass@(vv[1:]**2))
    internal=float(mass@state.e)
    nuclear_energy=float(mass@(-(yy-y0)@MASS))*MEV_PER_U_CGS
    neutrinos=float(mass@z[:,NU_I])*uunit
    escaped=float(mass@z[:,RAD_I])*uunit
    residual=(kinetic+internal+neutrinos+escaped-initial_energy-nuclear_energy)/initial_energy
    ymean=mass@yy/mass_total
    record={'model':'radial Lagrangian pulse expansion',
      'stage':stage,'T0_K':temp0,'rho0_kg_m3':rho0_si,'R0_m':radius0,'zones':zones,
      'nuclear_enabled':nuclear,'ideal_gas_verification_profile':ideal_test,
      'initial_carbon13_number_fraction':carbon13_fraction,
      'reduced_three_species_A_network':reduced_a,
      'radiation_diffusion':radiation,'opacity_m2_kg':kappa,'rtol':rtol,
      'atol':atol,'initial_profile_description':None if initial_profile is None
         else initial_profile.get('description','supplied spatial inflow'),
      'artificial_viscosity_linear':art_linear,'artificial_viscosity_quadratic':art_quad,
      'initial_surface_radius_m':float(init_r[-1]/100),
      'initial_internal_energy_J':initial_internal*1e-7,
      'initial_kinetic_energy_J':initial_kinetic*1e-7,
      'elapsed_s':float(sol.t[-1]*tunit),'elapsed_sound_times':float(sol.t[-1]),
      'mass_kg':mass_total/1000,'initial_energy_J':initial_energy*1e-7,
      'nuclear_energy_J':nuclear_energy*1e-7,'final_internal_J':internal*1e-7,
      'final_kinetic_J':kinetic*1e-7,'weak_energy_escaped_upper_J':neutrinos*1e-7,
      'radiation_escaped_J':escaped*1e-7,
      'deposited_nuclear_over_initial_energy':(nuclear_energy-neutrinos)/initial_energy,
      'energy_residual_fraction_initial':residual,
      'max_cell_baryon_residual':float(np.max(np.abs(yy@AA-1))),
      'minimum_mass_fraction':float(np.min(yy*AA)),
      'desired_product_per_initial_catalyst':float(ymean[species.index(target)]/cats0),
      'peak_temperature_K':float(np.max(sol.y[TI::NV])*temp0),
      'surface_expansion_factor':float(rr[-1]/R0),
      'output_mass_fractions':dict(zip(species,(AA*ymean).tolist())),
      'profile':{'outer_radius_m':(rr[1:]/100).tolist(),'velocity_m_s':(vv[1:]/100).tolist(),
         'rho_kg_m3':(rho*1000).tolist(),'T_K':temp.tolist(),
         'desired_product_per_initial_catalyst':(yy[:,species.index(target)]/cats0).tolist(),
         'mass_fraction_of_total':(mass/mass_total).tolist()},
      'success':bool(sol.success),'message':sol.message,'nfev':sol.nfev,'njev':sol.njev,
      'termination_reason':('solver_failure' if not sol.success else 'temperature_guard' if len(sol.t_events[0]) else 'final_time'),
      'viscosity_mode':viscosity_mode,
      'integrated_nuclear_energy_J':float(mass@z[:,HEAT_I])*uunit*1e-7,
      'integrated_charge_change_per_baryon':float(mass@z[:,CHARGE_I]/mass_total),
      'initial_charge_per_baryon':float(y0@ZZ),
      'final_charge_per_baryon':float(ymean@ZZ),
      'state_variables':species+['T/T0','r/R0','v/cs0','weak_energy/uunit','radiated_energy/uunit','nuclear_energy/uunit','weak_charge_change'],
      'velocity_unit_cm_s':cs0,'energy_unit_erg_g':uunit,'time_unit_s':tunit,
      'wall_seconds':time.monotonic()-begin,
      'limitations':[('initially compressed fuel; no implosion or driver'
         if initial_profile is None else 'specified convergent flow; its preparation and driver are not modeled'),
        'spherical symmetry; no multidimensional instabilities or mixing',
        'no screening or Coulomb EOS correction',
        ('three-species N14+p capture and reverse only; no weak decay' if reduced_a
         else '27-isotope network truncated at Mg'),
        'constant opacity if radiation enabled; no electron conduction',
        'hot beta Q entirely assigned to neutrinos, conservatively',
        'no catalyst processing, cold-decay, collection or plant-mass model']}
    if ideal_test:
        s=sol.t[-1];scale=np.sqrt(1+1.2*s*s)
        expectedr=x0[1:]*scale
        expectedv=x0[1:]*(1.2*s/scale)
        expectedT=temp_profile/scale**2
        record['analytic_error']={
          'radius_max_fraction_R0':float(np.max(np.abs(z[:,RI]-expectedr))),
          'velocity_max_fraction_cs0':float(np.max(np.abs(z[:,VI]-expectedv))),
          'temperature_mass_weighted_relative_L1':float(mass@np.abs(z[:,TI]-expectedT)/(mass@expectedT))}
    if save_trace:
        Path(save_trace).with_suffix('.pretrace.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        selected=np.unique(np.rint(np.linspace(0,len(sol.t)-1,min(trace_samples,len(sol.t)))).astype(int))
        np.savez_compressed(save_trace,time_s=sol.t[selected]*tunit,state=sol.y[:,selected],
                            species=species,variables_per_zone=NV,zones=zones,cell_mass_kg=mass/1000,
                            velocity_unit_cm_s=cs0,energy_unit_erg_g=uunit,time_unit_s=tunit)
    return record


