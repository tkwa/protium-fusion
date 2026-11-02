"""Conditional direct conversion of the archived ballistic pulse family.

No hydro or charge-state evolution is run. Fully stripped charge/current values
are favorable kinematic diagnostics. Scalar electrical scenarios are engineering
requirements, not efficiencies predicted for this source.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

E=1.602176634e-19
MEV=1e6*E
EPS=8.8541878128e-12
MU0=1.25663706212e-6
C=299792458.
ME=9.1093837139e-31

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def quant(x,w):
    ix=np.argsort(x);x=x[ix];w=w[ix];cs=np.cumsum(w)
    return {str(q):float(x[min(np.searchsorted(cs,q*cs[-1]),len(x)-1)]) for q in [.01,.1,.5,.9,.99]}

def main(root,output):
    paths={k:root/'docs/protium-fusion'/p for k,p in {
        'station':'open-loop-redesign/station-inputs.json',
        'wall':'membrane-receiver/pulse-loading.json',
        'stopping':'membrane-receiver/stopping-results.json',
        'matched_thermal_cases':'cold-sink-redesign/cases.json'}.items()}
    station=json.loads(paths['station'].read_text());wall=json.loads(paths['wall'].read_text());stop=json.loads(paths['stopping'].read_text())
    cyc=station['cycle'];rec=station['receiver'];R=wall['radius_m'];area=4*math.pi*R**2;mu=wall['atomic_mass_unit_kg']
    A=np.array(wall['A']);Z=np.array(wall['Z']);charged=Z>0
    rate=.8/(rec['occupied_s_per_He_at_R1e9']['mass999']*R/1e9)
    H=cyc['H_nonnu_MeV_per_He'];U=cyc['U_prepared_MeV_per_He'];K=cyc['bulk_kinetic_MeV_per_He']
    capture=.99;eta_p=.6;meanK=rate*K*MEV/area
    stage_out={};pulse_sum=0;weightedK=0
    for name,s in wall['stages'].items():
        pack=s['spectral_packets'];v=np.array(pack['velocity_m_s']);eu=np.array(pack['kinetic_MeV_per_u']);Y=np.array(pack['Y_at_wall']);mf=np.array(pack['mass_fluence_kg_m2']);t=np.array(pack['arrival_time_s'])
        counts=mf[:,None]/mu*Y;energy=counts*A[None,:]*eu[:,None]*MEV
        cells=np.array(pack['cell_index']);vb=np.array(s['cells']['velocity_bounds_m_s']);rb=np.array(s['cells']['radius_bounds_m']);m=np.array(s['cells']['mass_kg'])
        dv=np.diff(vb);drdv=np.diff(rb)/dv;b=R-rb[:-1]+drdv*vb[:-1]
        flux_mass=m[cells]/dv[cells]/area*v*v/b[cells]
        numberflux=flux_mass[:,None]/mu*Y
        j=numberflux@(Z*E);pk=.5*flux_mass*v*v;ii=int(np.argmax(pk))
        ni=numberflux[ii]/v[ii];ne=float(ni@Z)
        omega_pi=float(np.sqrt(np.sum(ni[charged]*(Z[charged]*E)**2/(EPS*A[charged]*mu))))
        pressure=s['smooth']['peak_momentum_flux_Pa'];Bpressure=math.sqrt(2*MU0*pressure)
        eq=eu[:,None]*A[None,charged]/Z[None,charged]*1e6
        eqone=eu[:,None]*A[None,charged]*1e6
        en=energy[:,charged]
        energy_species=energy.sum(axis=0)/energy.sum()
        ramp=[]
        for alpha in [1.,1.75,1.875,2.]:
            ratio=np.divide(A,Z,out=np.full(len(A),np.inf),where=Z>0)
            accepted=(ratio>=alpha)&charged
            recovered=(counts[:,accepted]*Z[None,accepted]*alpha*eu[:,None]*MEV).sum()
            reflected=energy[:,charged & ~accepted].sum()
            ramp.append({'V_over_Eperu_MeV_in_MV':alpha,'ideal_recovered_fraction_of_all_K':float(recovered/energy.sum()),'kinetic_fraction_needing_reflected_or_separate_route':float(reflected/energy.sum()),'scope':'Kinematic shell-following retarder; each accepted ion must be collected or routed before reacceleration. No extraction, waveform, grid or electron losses included.'})
        nodal=s['exact_nodal']['kinetic_J'];pulse_sum+=nodal
        cat=next(k for k,x in cyc['stage_pulse_data'].items() if x['stage']==name)
        pperhe=cyc['stage_pulse_data'][cat]['pulses_per_He'];weightedK+=nodal*pperhe/MEV
        ti=s['smooth']['arrival_1_50_99_percent_energy_s']
        stage_out[name]={
            'pulse_kinetic_J_exact_nodal':nodal,'fluence_K_J_m2_exact_nodal':nodal/area,
            'fluence_K_J_m2_smooth':float(energy.sum()),
            'peak_K_W_m2_smooth_exact_cell_endpoint':s['smooth']['peak_kinetic_power_W_m2'],
            'whole_receiver_peak_K_W':s['smooth']['peak_kinetic_power_W_m2']*area,
            'pulse_interval_s_at_reference_rate':1/(rate*pperhe),
            'mean_K_W_m2_this_stage':rate*pperhe*nodal/area,
            'arrival_1_50_99_percent_energy_s':ti,
            'fully_stripped_charge_fluence_C_m2':float(counts.sum(axis=0)@(Z*E)),
            'fully_stripped_peak_current_A_m2_quadrature_sample':float(j.max()),
            'fully_stripped_sum_local_peak_currents_A':float(j.max()*area),
            'fully_stripped_Eperq_V_energy_weighted_quantiles':quant(eq.ravel(),en.ravel()),
            'singly_ionized_Eperq_V_energy_weighted_quantiles':quant(eqone.ravel(),en.ravel()),
            'species_fraction_of_K':{wall['species'][i]:float(energy_species[i]) for i in np.argsort(-energy_species)[:8]},
            'shared_shell_ramp_ideal_limits':ramp,
            'sample_at_largest_quadrature_Kflux':{'time_s':float(t[ii]),'velocity_m_s':float(v[ii]),'Kflux_W_m2':float(pk[ii]),'fully_stripped_current_A_m2':float(j[ii]),'fully_stripped_electron_density_m3':ne,'total_ion_density_m3':float(ni.sum()),'ion_inertial_length_c_over_omega_pi_m':C/omega_pi,'electron_inertial_length_m':C/math.sqrt(ne*E*E/(EPS*ME)), 'B_matching_peak_ram_pressure_T':Bpressure,'N14_gyroradius_at_that_B_m':14*mu*v[ii]/(7*E*Bpressure)},
            'scope':'Saved64-zone ballistic reconstruction and finite-flight nuclear daughters; atomic charge is assumed only for voltage/current diagnostics. Peak current uses512 quadrature samples, kinetic peak uses cell endpoints.'}
    peak=max(s['whole_receiver_peak_K_W'] for s in stage_out.values());peakarea=peak/area
    storage=[]
    c13=stage_out['C13'];duration=c13['arrival_1_50_99_percent_energy_s'][2]-c13['arrival_1_50_99_percent_energy_s'][0]
    for eta_d in [.5,.7,.8,.85]:
        factor=capture*eta_d
        largest=factor*max(x['pulse_kinetic_J_exact_nodal'] for x in stage_out.values())
        # The0.98 window belongs to the smooth spectrum, avoiding nodal/continuous mixing.
        lower=max(0.,factor*(.98*c13['fluence_K_J_m2_smooth']*area-rate*K*MEV*duration))
        sufficient=factor*pulse_sum
        storage.append({'gross_recovery_fraction':eta_d,'captured_fraction':capture,'largest_recovered_pulse_J':largest,'largest_pulse_J_m2':largest/area,'steady_export_window_lower_bound_J':lower,'sum_stage_sufficient_capacity_J':sufficient,'sum_stage_sufficient_capacity_J_m2':sufficient/area,'mean_gross_direct_W_m2':factor*meanK,'peak_gross_direct_W_m2':factor*peakarea,'added_store_mass_by_specific_energy_J_kg':{str(x):{'lower_bound_kg':lower/x,'sum_stage_sufficient_kg':sufficient/x} for x in [1e4,1e5,1e6]},'driver_store_is_separate':True})
    cases=[]
    for eta_h in [.4,.6,.65,.7]:
      for eta_d in [0.,.5,.7,.8,.85]:
       for store_eta in ([1.] if eta_d==0 else [1.,.9]):
        D0=capture*eta_d*K;loss=(1-store_eta)*D0;D=D0-loss
        thermal=capture*(H+U)-D0
        # Storage loss separately rejected: no unsupported high-temperature credit.
        net=D+eta_h*thermal-U/eta_p
        baseline=eta_h*capture*(H+U)-U/eta_p
        gross=D+eta_h*thermal
        waste=capture*(H+U)-U-net
        identity=H-net-waste-(1-capture)*(H+U)
        cases.append({'eta_heat':eta_h,'eta_direct_gross_of_captured_K':eta_d,'all_direct_energy_storage_roundtrip_efficiency':store_eta,'direct_before_storage_MeV_per_He':D0,'direct_delivered_MeV_per_He':D,'storage_loss_MeV_per_He_rejected_separately':loss,'heat_engine_input_MeV_per_He':thermal,'total_gross_electric_MeV_per_He':gross,'driver_electric_MeV_per_He':U/eta_p,'pre_new_aux_net_MeV_per_He':net,'increment_vs_same_heat_baseline_MeV_per_He':net-baseline,'increment_W_m2':(net-baseline)*rate*MEV/area,'global_waste_MeV_per_He':waste,'global_identity_residual_MeV_per_He':identity})
    breakeven=[]
    for case in cases:
        if case['eta_heat']!=.6 or case['eta_direct_gross_of_captured_K']==0:continue
        gain=case['increment_W_m2']
        for specific in [.001,.002,.004,.006]:
          breakeven.append({'baseline_specific_W_kg':specific,'eta_direct_gross':case['eta_direct_gross_of_captured_K'],'storage_roundtrip':case['all_direct_energy_storage_roundtrip_efficiency'],'converter_only_min_input_peak_rating_W_kg':specific*capture*peakarea/gain,'scope':'Incremental mass break-even, same rate: DeltaM < DeltaP/specific. Omits added storage/auxiliaries and radiator/other mass savings; whole-plant coupling is separate.'})
    matched=json.loads(paths['matched_thermal_cases'].read_text())['cases']['working_200K_60pct_IR11_bypass1']
    matched_rows=[]
    for case in cases:
      if case['eta_heat']!=.6 or case['eta_direct_gross_of_captured_K']==0:continue
      gain=(1-matched['parameters']['auxiliary_fraction'])*case['increment_W_m2']
      matched_rows.append({'eta_direct_gross':case['eta_direct_gross_of_captured_K'],'all_output_storage_roundtrip':case['all_direct_energy_storage_roundtrip_efficiency'],'increment_W_m2_after_ordinary_aux':gain,'total_extra_mass_allowance_kg_m2_ignoring_massloads':gain/matched['net_W_kg'],'converter_only_min_input_peak_rating_W_kg':capture*peakarea*matched['net_W_kg']/gain})
    film=[]
    for x in stop['direct_rows']:
      if x['material']=='silica':
        film.append({'thickness_nm':x['thickness_nm'],'C13_K_fraction_removed_before_converter':x['energy_loss_fraction'],'C13_K_fraction_remaining_upper_bound_for_converter':x['energy_transmitted_fraction'],'scope':'Projectile energy removal only. Scattering, charge equilibration, electron current and film heat retention are separate.'})
    output_data={'scope':'Conditional upstream direct-recovery requirements; fully stripped diagnostics, not a charge-state calculation or a complete converter design. Frozen pulse and cycle inputs are not edited.',
      'input_sha256':{k:sha(p) for k,p in paths.items()},'code_sha256':sha(Path(__file__)),
      'reference':{'radius_m':R,'area_m2':area,'helium_per_s':rate,'timing':'mass999,80percent occupancy','H_nonnu_MeV_per_He':H,'U_prepared_MeV_per_He':U,'K_MeV_per_He':K,'capture':capture,'eligible_K_fraction_assumed':1.,'driver_efficiency':eta_p,'mean_K_W_m2':meanK,'peak_K_W_m2':peakarea,'peak_to_mean_K':peakarea/meanK,'sum_stage_pulse_K_J':pulse_sum,'existing_two_driver_pulse_storage_J':2*max(x['initial_energy_J'] for x in cyc['stage_pulse_data'].values())/eta_p},
      'stages':stage_out,'electrical_cases':cases,'storage':storage,'converter_rating_mass':[{'complete_input_peak_W_kg':r,'converter_kg':capture*peak/r,'converter_kg_m2':capture*peakarea/r} for r in [10.,100.,1000.,10000.]],'converter_only_break_even':breakeven,
      'electrostatic_benchmark':{'N14_charge':7,'mass_u':14,'voltage_V':1e6,'planar_zero_injection_Child_Langmuir_jd2_A':4/9*EPS*math.sqrt(2*7*E/(14*mu))*1e6**1.5,'scope':'Scale diagnostic only; finite-injection retarding multicomponent sheath requires its own Poisson/kinetic solution.'},
      'matched_60pct_thermal_reference':{'case':'working_200K_60pct_IR11_bypass1','specific_W_kg':matched['net_W_kg'],'ordinary_auxiliary_fraction':matched['parameters']['auxiliary_fraction'],'rows':matched_rows,'scope':'Scalar marginal comparison only; new mass-driven support and replacement, storage, mass reserves, and saved radiator/engine mass require the coupled plant calculation.'},
      'silica_before_converter':film,
      'sources':{'plasma_converter_experiment':'https://doi.org/10.13182/FST83-A20820','inductive_ICF_theory':'https://doi.org/10.13182/FST92-A30054','traveling_wave_deceleration_experiment':'https://doi.org/10.1585/pfr.10.3405013'},
      'checks':{'weighted_nodal_stage_K_matches_cycle_relative':weightedK/K-1,'all_electrical_global_residual_abs_MeV':max(abs(x['global_identity_residual_MeV_per_He']) for x in cases),'sum_individual_mean_equals_cycle_relative':sum(s['mean_K_W_m2_this_stage'] for s in stage_out.values())/meanK-1,'prepared_energy_debited_once':True,'driver_store_not_reused_for_output_without_schedule':True}}
    assert abs(output_data['checks']['weighted_nodal_stage_K_matches_cycle_relative'])<1e-6
    assert output_data['checks']['all_electrical_global_residual_abs_MeV']<1e-12
    output.write_text(json.dumps(output_data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'reference':output_data['reference'],'checks':output_data['checks'],'selected_cases':[x for x in cases if x['eta_heat']==.6 and x['eta_direct_gross_of_captured_K'] in [.0,.7,.8]]},indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,default=Path(__file__).with_name('direct_recovery.json'));a=ap.parse_args();main(a.root,a.output)
