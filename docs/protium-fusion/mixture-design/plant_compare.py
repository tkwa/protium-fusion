"""Coupled comparison of actual pulse mixtures at declared station constraints.

Uses new actual-composition material and receiver reductions. Rebuilds source,
storage and inventory allocations before calling the frozen thermal mass/power
ledger. Thermal surface proxies are reported separately from qualified limits.
"""
from pathlib import Path
import argparse,importlib.util,itertools,json,math,sys,hashlib
import numpy as np
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path.insert(0,str(ROOT/'component-design'))
import material_cycle as mc
spec=importlib.util.spec_from_file_location('frozen_cold_station',ROOT/'cold-sink-redesign/station.py')
cold=importlib.util.module_from_spec(spec);spec.loader.exec_module(cold)
YEAR=365.25*86400;RREF=1e9
CAT={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
BASE_H={'A':1.,'B1':1.,'B2':2.,'C13':1.}
BASE_DURATION={s:json.loads((ROOT/f'b2-completion/{fn}').read_text())['elapsed_s'] for s,fn in {'A':'regression-A-compression_cap.json','B1':'regression-B1-compression_cap.json','B2':'B2-64-compression-cap.json','C13':'regression-C13-compression_cap.json'}.items()}
BASE=json.loads((ROOT/'open-loop-redesign/operating-cases.json').read_text())['cases']['high_source_rating_surface']
WORK=json.loads((ROOT/'cold-sink-redesign/cases.json').read_text())['cases']['working_200K_60pct_IR11_bypass1']
COMMON=dict(WORK['parameters'])
# These fields are returned diagnostics rather than evaluate keyword arguments.
COMMON.pop('hx_internal_UA_per_hot_W_K',None)
COMMON['hx_internal_UA_per_hot_W_K']=WORK['parameters']['hx_internal_UA_per_hot_W_K']

def load_cases(runs,receivers):
    out=[]
    for path in sorted(Path(runs).rglob('*.json')):
        if any(path.stem.endswith('.'+x) for x in ['initial','inputs','pretrace','receiver']):continue
        try:d=json.loads(path.read_text())
        except (ValueError,UnicodeError):continue
        if not isinstance(d,dict) or 'mixture_inputs' not in d:continue
        rp=Path(receivers)/(path.stem+'.receiver.json')
        if not rp.exists():continue
        r=json.loads(rp.read_text())
        if not math.isclose(r.get('radius_m',0),RREF,rel_tol=1e-12):
            raise ValueError('Receiver cache must be evaluated at radius 1e9 m: '+str(rp))
        if not r.get('usable_ballistic') or not r.get('completed_requested_duration'):continue
        if r['input_sha256'].get(str(path))!=hashlib.sha256(path.read_bytes()).hexdigest():continue
        try:s=mc.load_stage(path)
        except (ValueError,ArithmeticError):continue
        out.append({'name':path.stem,'path':str(path),'record':d,'receiver':r,'stage_data':s})
    return out

def cycle_for(selected,rc=1.,rh=1.):
    return mc.solve_cycle({CAT[s['record']['stage']]:s['stage_data'] for s in selected},carbon_return=rc,hydrogen_return=rh,eta_heat=.6,eta_driver=.6,capture=.99)

def station(selected,cycle,rate,radius=1.1e9,*,cold_parameters=None,retained_CNO_residence_s=YEAR/2,gas_residence_s=3600.,absorber_sigma=5000.):
    pars=dict(COMMON)
    if cold_parameters:pars.update(cold_parameters)
    area=4*math.pi*radius**2
    peak=max(s['receiver']['source_rating']['peak_total_prepared_electric_W']*RREF/radius for s in selected)
    stored=2*max(s['record']['initial_energy_J'] for s in selected)/.6
    m=dict(BASE['mass_breakdown_kg'])
    m['receiver']=absorber_sigma*area;m['insulation']=100*area;m['processing_structure_thermal_service']=1000*area
    m['peak_driver_hardware']=peak/100.;m['two_pulse_storage']=stored/1e5
    m['average_driver_conditioning']=.1*cycle['U_prepared_MeV_per_He']*cold.MEV*rate/.6
    gross=cycle['gross_fuel_kg_per_He']*rate
    hydrogen_output=cycle['gross_H_output_atoms_per_He']*mc.atomic_kg('p')*rate
    cnoflow=sum(cycle['stage_attempts_per_He'][CAT[s['record']['stage']]]*sum(v*mc.atomic_kg(n) for n,v in s['stage_data']['cold_output'].items() if n in mc.DEFAULT_STATES) for s in selected)*rate
    m['working_inventory']=cnoflow*retained_CNO_residence_s+(hydrogen_output+cycle['fresh_hydrogen_kg_per_He']*rate)*gas_residence_s
    m['inflight_fuel_inventory']=sum(cycle['stage_pulse_data'][CAT[s['record']['stage']]]['gross_fuel_kg_per_He']*rate*s['receiver']['timing_at_radius']['mass999']['occupied_s']*radius/RREF for s in selected)
    m['fresh_feed_buffer']=(cycle['fresh_hydrogen_kg_per_He']+cycle['fresh_carbon_kg_per_He'])*rate*30*86400
    b=dict(BASE);b.update({'radius_m':radius,'helium_per_s':rate,'mass_breakdown_kg':m,'tidal_equivalent_mass_kg':BASE['tidal_equivalent_mass_kg']*(radius/BASE['radius_m'])**3})
    ans=cold.evaluate(b,cycle,**pars)
    if 'net_W' not in ans:return ans
    occupation=sum(cycle['stage_pulse_data'][CAT[s['record']['stage']]]['pulses_per_He']*s['receiver']['timing_at_radius']['mass999']['occupied_s'] for s in selected)*radius/RREF*rate
    sourceavg=cycle['U_prepared_MeV_per_He']*cold.MEV*rate/.6
    surface={s['record']['stage']:{'material_K':s['receiver']['thermal']['kinetic_plus_advected_endpoint_internal_peak_rise_K']*(RREF/radius)**2.5,'prompt_surface_K':s['receiver']['prompt']['radial_ray_arrival_surface_screen']['peak_surface_rise_K']*(RREF/radius)**2} for s in selected}
    maxsurface=max(x['material_K']+x['prompt_surface_K'] for x in surface.values())
    ans.update({'selected_cases':{s['record']['stage']:s['name'] for s in selected},'receiver_areal_mass_kg_m2':absorber_sigma,'sequential_occupancy_mass999':occupation,
      'peak_source_electric_W':peak,'source_average_to_peak_ratio':sourceavg/peak,'stored_electric_J':stored,
      'gross_fuel_kg_s':gross,'recovered_H_kg_s':hydrogen_output,'recovered_CNO_kg_s':cnoflow,
      'fresh_H_kg_s':cycle['fresh_hydrogen_kg_per_He']*rate,'fresh_C_kg_s':cycle['fresh_carbon_kg_per_He']*rate,
      'fresh_C_H_mass_ratio':cycle['fresh_C_H_mass_ratio'],'unused_007_C_allowance_kg_s':(.007-cycle['fresh_C_H_mass_ratio'])*cycle['fresh_hydrogen_kg_per_He']*rate,
      'pulse_interval_s':{k:1/(x['pulses_per_He']*rate) for k,x in cycle['stage_pulse_data'].items()},
      'surface_proxies':surface,'max_material_plus_surface_prompt_rise_K':maxsurface,
      'surface_proxy_scope':'Sum of separate maxima in constant-property carbon; prompt mathematical-surface deposition is an intentionally conservative reference, not a qualified physical temperature. Finite spectral penetration and delayed heat require separate treatment.',
      'inventory_residence_assumptions_s':{'retained_CNO':retained_CNO_residence_s,'local_H_loop':gas_residence_s},
      'unburned_H_processing_sensitivities':[{'J_per_kg':cost,'W':cost*hydrogen_output,'fraction_net_export':cost*hydrogen_output/ans['net_W']} for cost in [1e6,1e7,1e8,1e9]],
      'H_processing_J_kg_for_10pct_export':.1*ans['net_W']/hydrogen_output if hydrogen_output else None,
      'processing_cost_scope':'These named H-processing loads are compared with the retained ordinary-auxiliary allowance; they are not additional debits on top of an allowance already intended to cover processing. New requirements exceeding that allowance must be added explicitly.',
      'qualification':'Conditional integrated accounting point. Source formation, receiver recovery and guide/thermal hardware remain engineering candidates.'})
    reasons=list(ans['reasons'])
    if occupation>.8*(1+1e-10):reasons.append('main_window_occupancy')
    if sourceavg>peak*(1+1e-10):reasons.append('source_mean_above_peak')
    if cycle['fresh_C_H_mass_ratio']>.007:reasons.append('fresh_carbon_budget')
    ans['reasons']=reasons;ans['feasible']=not reasons
    return ans

def best_rate(selected,cycle,radius=1.1e9,*,absorber_sigma=5000.):
    occ=sum(cycle['stage_pulse_data'][CAT[s['record']['stage']]]['pulses_per_He']*s['receiver']['timing_at_radius']['mass999']['occupied_s'] for s in selected)*radius/RREF
    cap=.8/occ
    # Fixed installed source peak specification; its mean cannot exceed peak.
    peak=max(s['receiver']['source_rating']['peak_total_prepared_electric_W']*RREF/radius for s in selected)
    cap=min(cap,peak*.6/(cycle['U_prepared_MeV_per_He']*cold.MEV))
    hi=station(selected,cycle,cap,radius,absorber_sigma=absorber_sigma)
    if hi.get('radiator_area_fraction',0)>.9:
        lower=None
        for fraction in [.5,.25,.125,.0625,.03125,.015625]:
            lo=station(selected,cycle,cap*fraction,radius,absorber_sigma=absorber_sigma)
            if 'radiator_area_fraction' in lo and lo['radiator_area_fraction']<.9:
                lower=fraction;break
        if lower is None:return {'feasible':False,'reasons':['no_valid_thermal_rate_bracket']}
        cap=brentq(lambda f:station(selected,cycle,f*cap,radius,absorber_sigma=absorber_sigma)['radiator_area_fraction']-.9,lower,1.,xtol=1e-13)*cap*(1-1e-10)
        hi=station(selected,cycle,cap,radius,absorber_sigma=absorber_sigma)
    return hi

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--receivers',type=Path,default=P/'receiver-results');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--opacity',type=float,default=.001);ap.add_argument('--zones',type=int,default=64);ap.add_argument('--absorber-sigma',type=float,default=5000.);ap.add_argument('--duration-factor',type=float,default=1.);ap.add_argument('--constraint',choices=['fixed_energy','fixed_flow'],default='fixed_energy');a=ap.parse_args()
    allcases=load_cases(a.runs,a.receivers);families={s:[] for s in CAT}
    for c in allcases:
        d=c['record'];i=d['mixture_inputs'];s=d['stage']
        if d['zones']!=a.zones or i.get('temperature_factor',1)!=1 or i.get('prepared_energy_factor',1)!=1 or not math.isclose(i.get('duration_s',0),d['elapsed_s'],rel_tol=1e-10,abs_tol=1e-15):continue
        if not math.isclose(i.get('duration_s',0),BASE_DURATION[s]*a.duration_factor,rel_tol=1e-10,abs_tol=1e-15):continue
        if d['opacity_m2_kg']!=a.opacity:continue
        if i['constraint']!=a.constraint and not(i['hydrogen_per_catalyst']==BASE_H[s] and a.constraint=='fixed_flow'):continue
        if c['stage_data']['input_MeV']<=0:continue
        families[s].append(c)
    if any(not v for v in families.values()):raise ValueError('Incomplete stage family')
    rows=[];rejected=0
    for selected in itertools.product(*families.values()):
        try:cycle=cycle_for(selected)
        except (ValueError,ArithmeticError,np.linalg.LinAlgError):rejected+=1;continue
        if .6*.99*(cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He'])-cycle['U_prepared_MeV_per_He']/.6<=0:rejected+=1;continue
        plant=best_rate(selected,cycle,absorber_sigma=a.absorber_sigma)
        if plant.get('feasible'):rows.append((plant['net_W_kg'],plant,cycle))
        else:rejected+=1
    rows.sort(key=lambda x:x[0],reverse=True)
    baseline=[]
    for s,cc in families.items():
        options=[c for c in cc if c['record']['mixture_inputs']['hydrogen_per_catalyst']==BASE_H[s]]
        if options:baseline.append(options[0])
    controls={}
    if len(baseline)==4:
        cy=cycle_for(baseline);controls['closed_return_same_He_rate']=station(baseline,cy,BASE['helium_per_s'],absorber_sigma=a.absorber_sigma);controls['closed_return_capacity_rate']=best_rate(baseline,cy,absorber_sigma=a.absorber_sigma)
    summary=[]
    for score,p,c in rows:
        summary.append({'mW_kg':score*1000,'selected_cases':p['selected_cases'],'He_s':p['helium_per_s'],'net_W':p['net_W'],'mass_kg':p['mass_kg'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'prepared_MeV_He':c['U_prepared_MeV_per_He'],'fresh_C_H':c['fresh_C_H_mass_ratio'],'H_return_kg_s':p['recovered_H_kg_s'],'surface_proxy_rise_K':p['max_material_plus_surface_prompt_rise_K']})
    output={'scope':__doc__,'comparison':{'radius_m':1.1e9,'receiver_areal_mass_kg_m2':a.absorber_sigma,'radiator_K':200.,'heat_conversion':.6,'driver_efficiency':.6,'maximum_main_window_occupancy':.8,'maximum_radiator_fraction':.9,'eventual_material_return':1.,'opacity_m2_kg':a.opacity,'mixture_constraint':a.constraint,'zones':a.zones,'duration_factor':a.duration_factor},'stage_candidates':{s:[c['name'] for c in cc] for s,cc in families.items()},'feasible_combinations':len(rows),'rejected_combinations':rejected,'controls':controls,'ranked':summary,'top_detailed':[{'plant':p,'cycle':c} for _,p,c in rows[:10]],'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'material_code_sha256':mc.hashfile(mc.__file__)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'feasible':len(rows),'rejected':rejected,'best':summary[0] if summary else None,'out':str(a.out)},indent=2))
if __name__=='__main__':main()
