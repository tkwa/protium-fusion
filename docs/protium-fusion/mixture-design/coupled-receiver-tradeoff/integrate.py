"""Coupled Geant4 response and consistently rebuilt plant accounts.

Nested transport responses describe the all-cold response basis before each
row's explicit plant routing.
Transport scenarios and extra missed-hot sensitivities are not model probability
intervals. Material recovery is a separate, closed-reservoir account.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, sys
HERE=Path(__file__).resolve().parent
MIX=HERE.parent
sys.path.insert(0,str(MIX))
import receiver_tradeoff as rt
import plant_compare as pc
import station_generalized as station
TRANSPORT=MIX.parent/'component-design'/'coupled-photon-transport'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def calculate(inputs_path,plant_path,spectrum_path):
    inputs=read(inputs_path);base=read(plant_path)['top_detailed'][0]
    selected=[inputs['cases'][label] for label in base['summary']['cases'].values()]
    cycle=pc.cycle_for(selected);spectrum=read(spectrum_path)
    checks={
        'all_128_zones':all(s['record']['zones']==128 for s in selected),
        'prepared_energy_same_family':abs(cycle['U_prepared_MeV_per_He']/base['cycle']['U_prepared_MeV_per_He']-1)<1e-12,
        'photon_total_same_family':abs(spectrum['thermal_boundary_MeV_per_He']/(cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He'])-1)<1e-12,
    }
    assert all(checks.values()),checks
    module=rt.load_module('coupled_transport',TRANSPORT/'reweight.py');basis=read(TRANSPORT/'basis.json')
    original=pc.best_rate(selected,cycle,absorber_sigma=5000.)
    pc.cold=station
    default=rt.optimized(selected,cycle,5000.,.99,1.)
    for k in ['net_W','mass_kg','helium_per_s','cold_rejection_W','radiator_area_fraction','driver_W','support_W','replacement_W']:
        checks['default_reproduction_'+k]=bool(abs(default[k]/original[k]-1)<1e-11)
    assert all(checks.values()),checks
    rows=[]
    for sigma in [250.,500.,750.,1000.,2000.]:
        for prompt in ['EM20','all_to_cold']:
            response=module.evaluate(basis,spectrum,sigma,scenario='broad_D1p10',prompt_policy=prompt,cold_intercept_fraction=1.,hot_bypass_ppm=6.4)
            for extra in [0.,.005,.01]:
                loss=response['leakage_fraction']+extra
                assert 0<=loss<1
                for routing in ['cold_equipment','external_escape']:
                    if prompt=='all_to_cold' and routing=='external_escape':continue
                    capture,hot=(1.,1-loss) if routing=='cold_equipment' else (1-loss,1.)
                    p=rt.optimized(selected,cycle,sigma,capture,hot)
                    assert p['feasible'],p.get('reasons')
                    assert max(abs(v) for v in p['closure_relative_residuals'].values())<1e-10
                    assert abs(p['net_W']/p['mass_kg']-p['net_W_kg'])<1e-15
                    rows.append({'receiver_sigma_kg_m2':sigma,'prompt_policy':prompt,'transport_estimated_missed_hot_fraction_HplusU':response['leakage_fraction'],'transport_sampling_80pct_half_width':response['sampling_80pct_half_width'],'additional_missed_hot_fraction_HplusU':extra,'missed_hot_routing':routing,'total_hot_fraction':1-loss,'mW_kg':1000*p['net_W_kg'],'net_W':p['net_W'],'mass_kg':p['mass_kg'],'He_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'heat_in_cold_equipment_W':p['cold_heat_components_W']['captured_source_ineligible_for_hot_engine'],'escape_W':p['escaping_source_W'],'transport_response_before_plant_routing':response,'plant':p})
    return {'scope':__doc__,'family':spectrum['family'],'selected_cases':base['summary']['cases'],'cycle_reference_capture_099':cycle,'default_reproduction_checks':checks,'default_reproduction':default,'rows':rows,'input_sha256':{str(p):sha(p) for p in [inputs_path,plant_path,spectrum_path,TRANSPORT/'basis.json',TRANSPORT/'reweight.py',Path(rt.__file__),Path(pc.__file__),Path(station.__file__)]},'code_sha256':sha(__file__),'interpretation':[
        'The cycle reference contains original 99% heat-capture electrical comparison fields. Actual row exports and heat partitions use each row plant parameters; the nuclear, isotope and prepared-energy quantities are common.',
        'The Geant4 hard-source response is a conditional simulated result; its sampling interval excludes spectrum, geometry, cross-section and apparatus uncertainty.',
        'All-to-cold assigns the entire saved prompt-radiation budget to cold equipment once. Additional 0, 0.5%, and 1% missed-hot fractions are explicit sensitivities. The separate 1 W/m² thermal bypass remains charged.',
        'All missed-hot energy is charged to cold equipment in the primary comparison. External-escape rows are a different geometry scenario.',
        'No ordinary material disappearance is inferred from missed hot heat, a late arrival, implantation, or imperfect first-pass processing.',
        'Each row rebuilds cadence, source rating, radiator, inventories, services and reserve using the same source family. Thin receiver support, pulse surface heating, guide protection and cold-source delivery remain separate qualification conditions.'
    ]}
def calculate_final(family_path,basis_path=None,columns=None):
    family=read(family_path)
    selected=family['plant_selected'];cycle=family['precursor_cycle'];hot_cycle=family['hot_cycle']
    spectrum=family.get('spectrum') or read(Path(family_path).with_name('spectrum.json'))
    hhot=hot_cycle['H_nonnu_MeV_per_He']+hot_cycle['U_prepared_MeV_per_He']
    hpre=cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He']
    prep=family['preparation_escaped_MeV_per_He']
    checks={
        'all_128_zones':all(v['record']['zones']==128 for v in selected),
        'same_nuclear_energy':abs(cycle['H_nonnu_MeV_per_He']/hot_cycle['H_nonnu_MeV_per_He']-1)<1e-12,
        'spectrum_matches_hot_family':abs(spectrum['thermal_boundary_MeV_per_He']/hhot-1)<1e-12,
        'precursor_increment_matches_preparation_escape':abs((hpre-hhot)-prep)<1e-9,
        'source_storage_view_explicit':all(abs(v['record']['initial_energy_J']/v.get('source_prepared_energy_J',v['record']['initial_energy_J'])-1)<1e-12 for v in selected),
    }
    assert all(checks.values()),checks
    pc.cold=station
    module=rt.load_module('coupled_transport_final',TRANSPORT/'reweight.py');basis_path=Path(basis_path) if basis_path else TRANSPORT/'basis.json';basis=read(basis_path)
    rows=[]
    for sigma in (columns or [250.,500.,750.,1000.,2000.]):
        for prompt in ['EM20','all_to_cold']:
            response=module.evaluate(basis,spectrum,sigma,scenario='broad_D1p10',prompt_policy=prompt,cold_intercept_fraction=1.,hot_bypass_ppm=6.4)
            missed=response['leakage_fraction']*hhot+prep
            for extra in [0.,.005,.01]:
                loss=missed/hpre+extra
                assert 0<=loss<1
                for routing in ['cold_equipment','external_escape']:
                    if prompt=='all_to_cold' and routing=='external_escape':continue
                    capture,hot=(1.,1-loss) if routing=='cold_equipment' else (1-loss,1.)
                    p=rt.optimized(selected,cycle,sigma,capture,hot)
                    assert p['feasible'],p.get('reasons')
                    assert max(abs(v) for v in p['closure_relative_residuals'].values())<1e-10
                    known_preparation=prep*station.MEV*p['helium_per_s']
                    rows.append({'receiver_sigma_kg_m2':sigma,'prompt_policy':prompt,'transport_estimated_missed_hot_fraction_HplusUhot':response['leakage_fraction'],'transport_sampling_80pct_half_width_HplusUhot':response['sampling_80pct_half_width'],'additional_missed_hot_fraction_HplusUprecursor':extra,'preparation_radiation_MeV_He':prep,'preparation_radiation_destination':routing,'preparation_radiation_W':known_preparation,'missed_hot_routing':routing,'total_hot_fraction_HplusUprecursor':1-loss,'mW_kg':1000*p['net_W_kg'],'net_W':p['net_W'],'mass_kg':p['mass_kg'],'He_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'heat_in_cold_equipment_W':p['cold_heat_components_W']['captured_source_ineligible_for_hot_engine'],'escape_W':p['escaping_source_W'],'transport_response_before_plant_routing':response,'plant':p})
    return {'scope':__doc__,'family':family.get('family',spectrum['family']),'hot_cycle':hot_cycle,'precursor_cycle':cycle,'preparation_escaped_MeV_per_He':prep,'preparation_energy_increment_MeV_per_He':family['preparation_energy_increment_MeV_per_He'],'family_checks':checks,'rows':rows,'input_sha256':{str(p):sha(p) for p in [family_path,Path(family_path).with_name('spectrum.json'),basis_path,TRANSPORT/'reweight.py',Path(rt.__file__),Path(pc.__file__),Path(station.__file__)]},'code_sha256':sha(__file__),'interpretation':[
        'Pulse transport weights use the actual hot family. The source electrical and storage accounts use the explicit precursor input. Already escaped preparation radiation is allocated exactly once.',
        'The additional 0, 0.5%, and 1% missed-hot sensitivities multiply the complete non-neutrino nuclear plus precursor input energy. They are scenarios, not imposed physical requirements.',
        'Primary rows send all missed hot energy, including preparation radiation, to cold equipment. External-escape rows require a different geometry.',
        'Geant4 sampling intervals do not include uncertain source geometry, spectrum, opacity or component performance.',
        'Whole-plant mass and heat accounts are conditional on source, receiver and support hardware meeting their stated requirements. No material disappearance is inferred from energy leakage.'
    ]}
def main():
    a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path);a.add_argument('--plant',type=Path);a.add_argument('--spectrum',type=Path);a.add_argument('--family',type=Path);a.add_argument('--basis',type=Path);a.add_argument('--columns');a.add_argument('--out',type=Path,required=True);args=a.parse_args()
    if args.family:
        result=calculate_final(args.family,args.basis,[float(x) for x in args.columns.split(',')] if args.columns else None)
    else:
        if not all([args.inputs,args.plant,args.spectrum]):a.error('Supply --family or all of --inputs, --plant, --spectrum')
        result=calculate(args.inputs,args.plant,args.spectrum)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':result.get('family_checks',result.get('default_reproduction_checks')),'rows':[{k:v for k,v in r.items() if k not in ['plant','transport_response_before_plant_routing']} for r in result['rows'] if r['missed_hot_routing']=='cold_equipment' and r.get('additional_missed_hot_fraction_HplusUprecursor',r.get('additional_missed_hot_fraction_HplusU'))==.005]},indent=2))
if __name__=='__main__':main()
