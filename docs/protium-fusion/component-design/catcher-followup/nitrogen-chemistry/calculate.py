#!/usr/bin/env python3
"""Auditable inventory and acceptance calculations. No fitted nitrogen-release rates.
Run: python3 calculate.py. Standard library only. Preserved inputs.json is sufficient.
"""
from pathlib import Path
import json, math, hashlib
ROOT=Path(__file__).resolve().parent
KB_EV=8.617333262145e-5
EV_J=1.602176634e-19
NA=6.02214076e23

def release_kernel(times_s, temperatures_K, frequency_s, barrier_eV):
    """First-order trapped cohort; trapezoidal integral on caller's actual T(t).
    Each output is released fraction of THIS specified population. This is an
    interface, not a fit to implanted graphite or to an entire material inventory.
    Resolve rapid temperature excursions and don't interpolate microsecond peaks
    across a seconds grid. Add source cohorts and chemistry/routing separately.
    """
    assert len(times_s)==len(temperatures_K)>0 and frequency_s>=0 and barrier_eV>=0
    hazard=0.0; result=[0.0]
    for i in range(1,len(times_s)):
        dt=times_s[i]-times_s[i-1];assert dt>=0
        a=frequency_s*math.exp(-barrier_eV/(KB_EV*temperatures_K[i-1]))
        b=frequency_s*math.exp(-barrier_eV/(KB_EV*temperatures_K[i]))
        hazard+=dt*(a+b)/2
        result.append(-math.expm1(-hazard))
    return result

def second_order_release_fraction(n0_m3, integrated_k_m3):
    """Illustrative 2N->N2 loss: dn/dt=-2k(t)n², with no trapping/source.
    integrated_k_m3 = integral k(T)dt. Diffusion, defects and competing chemistry
    can invalidate this homogeneous-law reduction; no k is assigned here.
    """
    assert n0_m3>=0 and integrated_k_m3>=0
    x=2*n0_m3*integrated_k_m3
    return x/(1+x)

def main():
    d=json.loads((ROOT/'inputs.json').read_text());checks=[]
    def check(name, condition):
        checks.append({'name':name,'passed':bool(condition)})
        if not condition:raise AssertionError(name)
    area=d['area_m2'];rho=d['graphite_density_kg_m3'];year=d['year_s'];am=d['atomic_mass_kg'];nc=rho/am['c12']
    out={'scope':d['scope'],'source_sha256':d['source_sha256'],'input_sha256':hashlib.sha256((ROOT/'inputs.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'stages':{},'nominal_stopping_depth_scope':'The depth bins are straight CSDA estimates; no straggling, channeling, chemical redistribution or post-decay recoil is included. Concentrations are bin averages, with no release/rework.'}
    annualN=annualO=annualH=0.;annualBins=None;mean_atomic=0.
    stagekeys={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
    for name,x in d['stages'].items():
        spec=x['species'];n=sum(v['nuclei_m2'] for v in spec.values());o=spec['o14']['nuclei_m2']+spec['o15']['nuclei_m2'];ny=n*year/x['period_s']
        atomic_mass=sum((spec[nuc]['nuclei_m2']+spec[oxy]['nuclei_m2'])*am[nuc] for nuc,oxy in [('n14','o14'),('n15','o15')]);mean_atomic+=area*atomic_mass/x['period_s']
        ledger=sum(d['stage_fully_decayed_N_output_atomic_kg'][stagekeys[name]].values())
        relative=(atomic_mass*area-ledger)/ledger
        check(name+' wall O+N equals fully-decayed N ledger',abs(relative)<1e-11)
        bins=[]
        for i,b in enumerate(spec['n14']['CSDA_stopping_depth_bins']):
            lo,hi=b['depth_um'];cnt=sum(v['CSDA_stopping_depth_bins'][i]['nuclei_m2'] for v in spec.values());c=cnt/((hi-lo)*1e-6*nc)
            bins.append({'depth_um':[lo,hi],'N_nuclei_m2':cnt,'N_per_C_after_one_shot_no_release':c})
        check(name+' depth bins conserve N',abs(sum(b['N_nuclei_m2'] for b in bins)/n-1)<1e-12)
        if annualBins is None:annualBins=[{'depth_um':b['depth_um'],'N_nuclei_m2_year':0.} for b in bins]
        for target,b in zip(annualBins,bins):target['N_nuclei_m2_year']+=b['N_nuclei_m2']*year/x['period_s']
        out['stages'][name]={'N_nuclei_cm2_pulse':n/1e4,'N_mass_ug_m2_atomic':atomic_mass*1e9,'O_precursor_fraction_of_eventual_N_nuclei':o/n,'H_to_eventual_N_nuclei':x['H_nuclei_m2']/n,'pulses_year':year/x['period_s'],'period_s':x['period_s'],'ledger_relative_difference':relative,'depth_bins':bins,'species':{k:{'fraction_of_eventual_N':v['nuclei_m2']/n,'mean_incident_MeV':v['kinetic_J_m2']/v['nuclei_m2']/EV_J/1e6 if v['nuclei_m2'] else 0,'CSDA_depth_um_quantiles_number':v['CSDA_depth_um_quantiles_number']} for k,v in spec.items()}}
        annualN+=ny;annualO+=o*year/x['period_s'];annualH+=x['H_nuclei_m2']*year/x['period_s']
    for b in annualBins:
        lo,hi=b['depth_um'];b['N_per_C_after_year_no_release']=b['N_nuclei_m2_year']/((hi-lo)*1e-6*nc)
    totalmass=sum(d['mean_fully_decayed_N_output_atomic_kg_s'].values())
    check('mean atomic flow matches independent cold ledger',abs(mean_atomic/totalmass-1)<1e-11)
    out['annual']={'eventual_N_nuclei_cm2':annualN/1e4,'N_mass_mg_m2_atomic':totalmass*year/area*1e6,'O_precursor_fraction_of_eventual_N_nuclei':annualO/annualN,'H_to_eventual_N_nuclei':annualH/annualN,'no_extraction_depth_bins':annualBins,'uniform_layer_N_per_C_no_release':{str(um):annualN/(rho*um*1e-6/am['c12']) for um in [1,5,10,100,150]}}
    threshold=next(x for x in d['gas_N2_thresholds'] if x['global_capture_g']==.01)
    fmax=threshold['fraction_of_total_returned_N_entering_cavity_as_N2_for10percent_total_column_increase']
    out['gas_acceptance']={'scope':'Same isotope mixture and1000K molecular-flow model as frozen gas/chemistry.json;10% increase in total H2+He steady column, not a universal pellet requirement. Liberation and cavity routing are separate; no value for either is predicted.','threshold':threshold,'routing_scenarios':[{'liberated_fraction_as_N2':r,'maximum_fraction_of_liberated_N2_reaching_cavity':min(1,fmax/r),'N2_atomic_mass_processed_kg_s':r*totalmass,'column_increase_if_front_fraction_0p001_relative_to_baseline':.1*r*.001/fmax} for r in [1e-5,1e-4,.001,.01,.1,1.]]}
    out['process']={'fully_decayed_N_atomic_kg_s':totalmass,'hold_up_scenarios':[{'residence_s':t,'N_atomic_kg':t*totalmass,'N_kg_m2':t*totalmass/area} for t in [1,60,3600,86400,30*86400]],'electric_work_per_eV_per_N_atom_W_m2_at50pct_efficiency':annualN/year*EV_J/.5,'scope':'Residence is retained station inventory. Chemical work uses actual positive free-energy/process requirements; no atomization or bulk-carbon replacement is assumed.'}
    # Source normalization audit: extreme all-external-faces inventory upper bound.
    v_molar_ml=22414.0
    a_one_cm2=2*(2*2+2*1.5+2*1.5)/100 # mm² -> cm²
    a_ten=a_one_cm2*10
    li=[]
    for dose,tp,vol in [(5e17,680.6,.01684),(3e17,681.7,.01672),(5e16,687.5,.01303)]:
        vmax=dose*a_ten/2/NA*v_molar_ml
        li.append({'dose_N_cm2':dose,'peak_K':tp+273.15,'reported_N2_mL':vol,'ten_cuboid_all_faces_cm2':a_ten,'maximum_implanted_N2_mL_STP_under_that_area':vmax,'reported_to_area_inventory_ratio':vol/vmax})
    out['Li2014_normalization_audit']={'scope':'Uses stated ten2x2x1.5mm specimens, all faces implanted, ideal gas22.414L/mol. Actual exposed area, batch normalization, native-N/background contributions are not fully specified. Ratios diagnose why published integrated signals cannot be used here as implanted-N release fractions.','rows':li,'table3_inferred_area_cm2_from0p00294mL_and0p147L_m2':(.00294/1000)/.147*1e4}
    # Explicit conditional kinetic scale, not a fitted release fraction.
    out['conditional_first_order_barrier_screen']={'scope':'If a specified population follows k=nu exp(-E/kBT), these E values give10% release during a constant-temperature dwell. nu=1e13/s is an explicit attempt-frequency scenario; recombination, migration and trap populations remain to be supplied. No experimental E is inferred from a single TDS peak.','frequency_s':1e13,'rows':[{'T_K':T,'dwell_s':t,'E_eV_for10percent_release':KB_EV*T*math.log(1e13*t/(-math.log(.9)))} for T in [1000,1500] for t in [10,100,3600,86400]]}
    check('constant-temperature first-order kernel',abs(release_kernel([0,10],[1000,1000],1,0)[-1]-(-math.expm1(-10)))<1e-15)
    check('second-order kernel normalization',second_order_release_fraction(2,0.25)==.5)
    check('gas threshold closes',abs(threshold['N_mass_to_cavity_at_threshold_kg_s']/totalmass/fmax-1)<1e-14)
    check('one cuboid external area',abs(a_one_cm2-.2)<1e-15)
    check('annual concentrations nonnegative',all(b['N_per_C_after_year_no_release']>=0 for b in annualBins))
    out['temperature_scope']={'profiles':d['temperature_profiles'],'note':'Use actual chronological histories for future calibrated rate integration. Largest150um-skin sampled material temperature in these histories is1493.56K; separate post resistance raises current layout result to1496.03K. Prompt upper bounds have microsecond durations.1668K is not an hours-long dwell.'}
    (ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n');(ROOT/'validation.json').write_text(json.dumps({'passed':all(x['passed'] for x in checks),'checks':checks},indent=2)+'\n')
    print(json.dumps({'checks':len(checks),'passed':True,'annual':out['annual']['eventual_N_nuclei_cm2'],'oxygen_fraction':annualO/annualN,'max_annual_bin_N_per_C':max(b['N_per_C_after_year_no_release'] for b in annualBins),'process_work_coefficient_W_m2_per_eV':out['process']['electric_work_per_eV_per_N_atom_W_m2_at50pct_efficiency']},indent=2))
if __name__=='__main__':main()
