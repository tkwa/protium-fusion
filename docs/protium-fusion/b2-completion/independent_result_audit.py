"""Independent trace arithmetic and semidiscrete conservation audit.

Does not import radial_pulse or radial_completion. Reuses FastEOS only as the
constitutive relation; it is not an independent EOS or nuclear-rate validation.
All trace samples are checked. Sparse-time viscosity integrals are diagnostic.
"""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import pynucastro as pn

BASE=Path(__file__).resolve().parent.parent
if not (BASE/'fast_eos.py').exists():
    BASE=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
sys.path.insert(0,str(BASE))
from fast_eos import FastEOS
MU=1.66053906660e-27
MEV=1.602176634e-13
ERG_PER_G_PER_MEV_BARYON=MEV/MU*1e4

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def audit(result_path,trace_path,expected_sound_times=4.):
    result_path=Path(result_path);trace_path=Path(trace_path)
    d=json.loads(result_path.read_text()); f=np.load(trace_path)
    species=list(f['species']); nuclei=[pn.Nucleus(n) for n in species]
    A=np.array([n.A for n in nuclei]); Z=np.array([n.Z for n in nuclei])
    atomic=np.array([n.mass for n in nuclei]); excess=np.array([n.dm for n in nuclei])
    nn=len(species);nv=int(f['variables_per_zone']);nz=int(f['zones'])
    assert nv in [nn+5,nn+7],(nn,nv)
    times=f['time_s'];states=f['state'].T.reshape((-1,nz,nv))
    rho0=d['rho0_kg_m3']/1000;R0=d['R0_m']*100;T0=d['T0_K']
    M=4*np.pi/3*rho0*R0**3
    mass=M*np.diff(np.linspace(0,1,nz+1)**3)
    node_mass=np.r_[.5*(mass[:-1]+mass[1:]),.5*mass[-1]]
    initial_y=states[0,:,:nn]
    center=FastEOS().state(rho0,T0,initial_y[0]@Z,initial_y[0].sum())
    cs0=float(np.sqrt(center.gamma1*center.p/rho0));uunit=float(center.e)
    tunit=R0/cs0
    units={'mass_kg':M/1000,'velocity_unit_cm_s':cs0,'energy_unit_erg_g':uunit,'time_unit_s':tunit}
    metadata_comparison={k:float(abs(float(f[k])-v)/max(abs(v),1e-300)) for k,v in units.items() if k in f}
    if 'cell_mass_kg' in f:
        metadata_comparison['cell_mass_max_relative']=float(np.max(np.abs(f['cell_mass_kg']-mass/1000)/(mass/1000)))
    eos=FastEOS();rows=[]; all_x=[];all_baryon=[];all_integrated_q=[];all_integrated_z=[]
    modes=['legacy','volume_gate','divergence','compression_cap'];vis={k:[] for k in modes}
    forcework={k:[] for k in modes}; min_width=np.inf; min_rho=np.inf;min_temp=np.inf
    maximum_sampled_T=0.; max_allcell_charge_error=0.
    a=d.get('artificial_viscosity_linear',.5);b=d.get('artificial_viscosity_quadratic',2.)
    for t,z in zip(times,states):
        y=z[:,:nn];r=np.r_[0.,z[:,nn+1]*R0];v=np.r_[0.,z[:,nn+2]*cs0]
        width=np.diff(r);V=4*np.pi/3*np.diff(r**3);rho=mass/V;T=z[:,nn]*T0
        assert np.all(width>0),'shell crossing'
        clipped=np.maximum(y,0);norm=clipped@A
        es=eos.state(rho,T,(clipped@Z)/norm,clipped.sum(axis=1)/norm)
        K=.5*np.sum(node_mass*v[1:]**2);U=np.sum(mass*es.e)
        # Independently use tabulated mass excess, removing the shared A*u term.
        dy=y-initial_y
        qmass=-np.sum(mass*(dy@excess))*ERG_PER_G_PER_MEV_BARYON
        qfull=-np.sum(mass*(dy@atomic))*ERG_PER_G_PER_MEV_BARYON
        weak=np.sum(mass*z[:,nn+3])*uunit;rad=np.sum(mass*z[:,nn+4])*uunit
        charge=np.sum(mass*(dy@Z))/M
        row={'t_s':float(t),'K_J':K*1e-7,'U_J':U*1e-7,'nuclear_mass_excess_J':qmass*1e-7,
             'nuclear_full_atomic_mass_J':qfull*1e-7,'weak_J':weak*1e-7,'radiation_J':rad*1e-7,
             'nuclear_charge_change_per_baryon':charge}
        if nv==nn+7:
            qint=np.sum(mass*z[:,nn+5])*uunit;zint=np.sum(mass*z[:,nn+6])/M
            row.update(integrated_nuclear_J=qint*1e-7,integrated_charge_change_per_baryon=zint)
            all_integrated_q.append(qint-qmass);all_integrated_z.append(zint-charge)
            max_allcell_charge_error=max(max_allcell_charge_error,float(np.max(abs(z[:,nn+6]-dy@Z))))
        rows.append(row);all_x.append(y*A);all_baryon.append(y@A-1)
        min_width=min(min_width,float(width.min()));min_rho=min(min_rho,float(rho.min()));min_temp=min(min_temp,float(T.min()));maximum_sampled_T=max(maximum_sampled_T,float(T.max()))
        # Independent finite-volume reconstruction. Sum_i[-q_i*dV_i/dt]
        # must be nonnegative locally for scalar viscous entropy production.
        Vdot=4*np.pi*np.diff(r*r*v)
        radial_compress=np.maximum(-np.diff(v),0.)
        soundspeed=np.sqrt(np.maximum(es.gamma1*es.p/rho,0.))
        compresses={'legacy':radial_compress,
          'volume_gate':np.where(Vdot<0,radial_compress,0.),
          'divergence':np.maximum(-Vdot/V*width/3,0.),
          'compression_cap':np.minimum(radial_compress,np.maximum(-Vdot/V*width,0.))}
        for mode,c in compresses.items():
            q=rho*(b*c*c+a*soundspeed*c)
            viscous_power=-q*Vdot
            effective_p=es.p+q
            boundary_force=4*np.pi*r[1:]**2*(effective_p-np.r_[effective_p[1:],0.])
            kinetic_work=np.sum(boundary_force*v[1:]);internal_work=-np.sum(effective_p*Vdot)
            scale=max(np.sum(abs(boundary_force*v[1:])),np.sum(abs(effective_p*Vdot)),1.)
            forcework[mode].append(float((kinetic_work+internal_work)/scale))
            vis[mode].append({'t_s':float(t),'negative_power_W':float(np.sum(np.minimum(viscous_power,0.))*1e-7),
              'positive_power_W':float(np.sum(np.maximum(viscous_power,0.))*1e-7),
              'negative_cell_count':int(np.count_nonzero(viscous_power<0)),
              'minimum_specific_viscous_power_W_kg':float(np.min(viscous_power/mass)*1e-4)})
    E0=rows[0]['K_J']+rows[0]['U_J']; final=rows[-1];xmean=A*(mass@states[-1,:,:nn])/M
    residuals=[(r['K_J']+r['U_J']+r['weak_J']+r['radiation_J']-E0-r['nuclear_mass_excess_J'])/E0 for r in rows]
    reconstructed={'mass_kg':M/1000,'elapsed_s':float(times[-1]),'elapsed_sound_times':float(times[-1]/tunit),
      'initial_energy_J':E0,'initial_internal_energy_J':rows[0]['U_J'],'initial_kinetic_energy_J':rows[0]['K_J'],
      'nuclear_energy_J':final['nuclear_full_atomic_mass_J'],'final_internal_J':final['U_J'],
      'final_kinetic_J':final['K_J'],'weak_energy_escaped_upper_J':final['weak_J'],
      'radiation_escaped_J':final['radiation_J'],'output_mass_fractions':dict(zip(species,xmean.tolist()))}
    comparisons={k:float((v-d[k])/max(abs(v),abs(d[k]),1e-300)) for k,v in reconstructed.items() if k!='output_mass_fractions' and k in d}
    comparisons['mass_fraction_max_absolute']=float(max(abs(xmean[i]-d['output_mass_fractions'][s]) for i,s in enumerate(species)))
    mode=d.get('viscosity_mode','legacy')
    vsummary={}
    for m,rs in vis.items():
        vsummary[m]={'negative_cell_samples':sum(r['negative_cell_count'] for r in rs),
          'most_negative_total_viscous_power_W':min(r['negative_power_W'] for r in rs),
          'largest_positive_total_viscous_power_W':max(r['positive_power_W'] for r in rs),
          'sample_trapezoid_negative_viscous_work_J':float(np.trapezoid([r['negative_power_W'] for r in rs],times)),
          'max_force_internal_work_cancellation_relative':max(abs(x) for x in forcework[m])}
    completion=bool(d['success'] and d.get('termination_reason')=='final_time' and abs(times[-1]/tunit-expected_sound_times)<1e-9)
    report={'scope':'Independent raw-trace arithmetic, elemental conservation, and discrete pressure work audit. Shared EOS/mass library and rate convention; no independent rates or driver validation.',
      'inputs':{'result':str(result_path),'trace':str(trace_path),'result_sha256':digest(result_path),'trace_sha256':digest(trace_path),'eos_sha256':digest(BASE/'electron-eos.npz')},
      'samples':len(times),'units_reconstructed':units,'metadata_relative_errors':metadata_comparison,
      'reported_field_relative_errors':comparisons,'reconstructed_result':reconstructed,
      'completion':{'requested_sound_times':expected_sound_times,'complete':completion,'reported_success':d['success'],'reported_termination':d.get('termination_reason')},
      'conservation':{'max_baryon_absolute_all_cells_samples':float(np.max(np.abs(all_baryon))),
        'minimum_mass_fraction_all_cells_samples':float(np.min(all_x)),
        'max_energy_residual_fraction_initial_all_samples':max(abs(x) for x in residuals),
        'final_energy_residual_fraction_initial':residuals[-1],
        'nuclear_mass_excess_minus_full_atomic_J':final['nuclear_mass_excess_J']-final['nuclear_full_atomic_mass_J'],
        'final_nuclear_charge_change_per_baryon':final['nuclear_charge_change_per_baryon'],
        'compensating_leptonic_charge_change_per_baryon':-final['nuclear_charge_change_per_baryon'],
        'independent_integrated_heat_available':nv==nn+7,
        'max_integrated_heat_minus_mass_defect_fraction_initial':max(abs(x)*1e-7/E0 for x in all_integrated_q) if all_integrated_q else None,
        'max_integrated_weak_charge_minus_abundance_charge':max(abs(x) for x in all_integrated_z) if all_integrated_z else None,
        'max_cell_integrated_weak_charge_error':max_allcell_charge_error if nv==nn+7 else None},
      'geometry':{'minimum_shell_width_m':min_width/100,'minimum_density_kg_m3':min_rho*1000,'minimum_T_K':min_temp,'max_sampled_T_K':maximum_sampled_T,
        'final_minimum_velocity_m_s':float(np.min(states[-1,:,nn+2])*cs0/100),
        'final_velocity_decreases_outward':int(np.count_nonzero(np.diff(states[-1,:,nn+2])<0))},
      'viscosity_mode_used':mode,'viscosity_reconstruction':vsummary,'sample_energies':rows,
      'limits':['Nuclear charge is not conserved alone under beta decay; signed weak-charge accumulator closes its change.',
        'Integrated nuclear heating is an independent saved ODE variable, but shares reaction flows and mass/Q source with abundance integration.',
        'Force/work cancellation is an algebraic semidiscrete check; positive viscosity entropy is a distinct requirement.',
        'Sparse trace quadrature is not an accurate total viscous-work integral or a resolved shock study.',
        'Conservation does not establish correct opacity, EOS, nuclear rates, multidimensional dynamics, driver preparation, or collection.']}
    numerical_checks={
      'reported_fields_reproduce':max(abs(x) for x in comparisons.values())<1e-7,
      'metadata_units_reproduce':max(metadata_comparison.values(),default=0)<1e-12,
      'baryon_conservation_all_samples':report['conservation']['max_baryon_absolute_all_cells_samples']<1e-10,
      'roundoff_negativity_only':report['conservation']['minimum_mass_fraction_all_cells_samples']>-1e-7,
      'energy_residual_below_1e_4_input':report['conservation']['max_energy_residual_fraction_initial_all_samples']<1e-4,
      'semidiscrete_force_work_cancels':vsummary[mode]['max_force_internal_work_cancellation_relative']<1e-12,
      'viscous_entropy_sign_valid_at_samples':vsummary[mode]['negative_cell_samples']==0,
      'requested_final_time_completed':completion}
    if nv==nn+7:
        numerical_checks['integrated_nuclear_energy_closes']=report['conservation']['max_integrated_heat_minus_mass_defect_fraction_initial']<1e-8
        numerical_checks['integrated_weak_charge_closes']=report['conservation']['max_cell_integrated_weak_charge_error']<1e-10
    numerical_checks={k:bool(v) for k,v in numerical_checks.items()}
    report['checks']=numerical_checks
    report['all_checks_pass']=all(numerical_checks.values())
    return report

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--result',required=True);ap.add_argument('--trace',required=True);ap.add_argument('--output',required=True);ap.add_argument('--expected-sound-times',type=float,default=4.)
    a=ap.parse_args();r=audit(a.result,a.trace,a.expected_sound_times);Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:r[k] for k in ['completion','conservation','viscosity_mode_used','viscosity_reconstruction']},indent=2))
if __name__=='__main__':main()
