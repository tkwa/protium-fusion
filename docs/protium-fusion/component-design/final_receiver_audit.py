"""Independent exact-family plant balances and schedule-aware thermal carryover.
No hydrodynamics. Component response functions come from preserved audited kernels.
"""
from pathlib import Path
import ast,json,math,csv,hashlib
import numpy as np
P=Path(__file__).parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
source=Path('/home/tkwa/src/protium-fusion/output/final-coupled-tradeoff-extended.json');d=read(source);cycle=d['precursor_cycle'];hotcycle=d['hot_cycle'];MEV=1.602176634e-13;errors={};checks=0
S_hot=hotcycle['H_nonnu_MeV_per_He']+hotcycle['U_prepared_MeV_per_He'];S_pre=cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He'];prep=d['preparation_escaped_MeV_per_He'];assert abs(S_pre-S_hot-prep)<1e-9
for row in d['rows']:
    p=row['plant'];v=p['parameters'];M=p['mass_kg'];A=p['area_m2'];F=p['helium_per_s'];H=cycle['H_nonnu_MeV_per_He']*MEV*F;U=cycle['U_prepared_MeV_per_He']*MEV*F
    loss=(row['transport_estimated_missed_hot_fraction_HplusUhot']*S_hot+prep)/S_pre+row['additional_missed_hot_fraction_HplusUprecursor'];miss=loss*(H+U);hot=(1-loss)*(H+U)-A;gross=.6*hot;driver=U/.6;aux=.1*(gross-driver);net=gross-driver-aux-p['replacement_W']-p['support_W'];escape=miss if row['missed_hot_routing']=='external_escape' else 0
    coldsrc=miss if row['missed_hot_routing']=='cold_equipment' else 0;waste=H-escape-net
    pairs={'mass':(sum(p['mass_breakdown_kg'].values()),M),'hot':(hot,p['engine_hot_W']),'gross':(gross,p['gross_electric_W']),'driver':(driver,p['driver_W']),'aux':(aux,p['ordinary_auxiliary_W']),'net':(net,p['net_W']),'cold_components':(sum(p['cold_heat_components_W'].values()),waste),'radiator':(p['radiator_area_m2']*.9*5.670374419e-8*200**4,waste+11.0254*A),'cold_source':(coldsrc,row['heat_in_cold_equipment_W']),'escape':(escape,p['escaping_source_W']),'preparation_once':(prep*MEV*F,row['preparation_radiation_W']),'receiver':(row['receiver_sigma_kg_m2']*A,p['mass_breakdown_kg']['receiver']),'factory':(M*1e-4/(365.25*86400)*1e8,p['mass_breakdown_kg']['replacement_factory_capacity'])}
    for cat,sd in cycle['stage_pulse_data'].items():pairs['cadence_'+cat]=(p['pulse_interval_s'][cat]*sd['pulses_per_He']*F,1.)
    for k,(x,y) in pairs.items():errors[k]=max(errors.get(k,0),abs(x-y)/max(abs(x),abs(y),1));checks+=1
    assert p['sequential_occupancy_mass999']<=.8 and p['radiator_area_fraction']<=.9000000001
assert max(errors.values())<1e-11
# Reconstruct actual final stage spectra using the existing loader only; stop
# before its original receiver-grid/main/output. No original file is modified.
original=P/'source128_thin_thermal.py';prefix=original.read_text().split('\nrows=[]\nfor k,cp',1)[0];prefix=prefix.replace("P/'source-family-carbon-median128/integration.json'","P/'final-receiver-thermal-input.json'")
ns={'__file__':str(original),'__name__':'final_thermal_inputs_only'};exec(compile(prefix,str(original),'exec'),ns)
stages=ns['stages'];history=ns['history'];tail_bound=ns['tail_bound'];thermal=read(P/'final-receiver-thermal.json');combined=read(P/'final-receiver-combined.json');dispatch=[]
inputs={str(source):sha(source),str(original):sha(original)}
for sigma in [250.,350.]:
    events_path=P/'final-family-dispatch'/('sigma'+str(int(sigma)))/'events.csv';events=list(csv.DictReader(events_path.open()));inputs[str(events_path)]=sha(events_path)
    big=sorted(float(e['focus_s']) for e in events if e['stage']!='B1');G_observed=min(b-a for a,b in zip(big,big[1:]));delays=[float(e['dispatch_delay_s']) for e in events];assert min(delays)>-1e-6 and max(delays)<=86400
    windows={e['stage']:(float(e['incoming_s']),float(e['outgoing_s'])) for e in events};G_serial=min(x[1]+y[0] for i,x in windows.items() if i!='B1' for j,y in windows.items() if j!='B1');assert G_observed>=G_serial
    r=next(r for r in thermal['thermal_rows'] if r['sigma_kg_m2']==sigma);L=r['thickness_m'];k=r['k_W_m_K'];cp=r['cp_J_kg_K'];smallP=stages['B1']['minimum_slot'];smallcut=min(20000,smallP/2)
    bounds=[]
    for label,G in [('observed370day_minimum',G_observed),('serialized_main_window_guarantee',G_serial)]:
        cutoff=min(20000,G/2);oldbig=max(s['fluence'] for n,s in stages.items() if n!='B1')*history(L,k,cp,G,cutoff,'N');oldsmall=stages['B1']['fluence']*history(L,k,cp,smallP,smallcut,'N');tail=(sum(tail_bound(s,G,cutoff) for n,s in stages.items() if n!='B1')+tail_bound(stages['B1'],smallP,smallcut))*L/(3*k)
        delta=oldbig+oldsmall+tail-r['old_big_nonuniform_background_K']-r['old_B1_nonuniform_background_K']-r['late_tails_nonuniform_bound_K']
        temperatures={str(x['skin_k_W_m_K']):x['with_pickup_floor1000_and_bounded1day_dispatch_lateness_K']+delta for x in combined['rows'] if x['sigma_kg_m2']==sigma and x['proxy']=='surface_flux_equivalent_Teff'}
        bounds.append({'basis':label,'large_focus_minimum_s':G,'early_energy_cutoff_s':cutoff,'old_big_background_K':oldbig,'old_B1_background_K':oldsmall,'late_tail_bound_K':tail,'correction_to_saved_daily_guard_K':delta,'Teff_allphase_front_K_by_skin_conductivity':temperatures})
    dispatch.append({'sigma_kg_m2':sigma,'maximum_actual_dispatch_lateness_s':max(delays),'retained_bulk_lateness_bound_s':86400.,'event_count':len(events),'bounds':bounds})
point_path=P.parent/'mixture-design/selected-plant/point.json';point=read(point_path);match=next(r for r in d['rows'] if r['receiver_sigma_kg_m2']==350 and r['prompt_policy']=='EM20' and r['missed_hot_routing']=='cold_equipment' and r['additional_missed_hot_fraction_HplusUprecursor']==0);selected_errors={k:abs(point['plant'][k]/match['plant'][k]-1) for k in ['mass_kg','net_W','helium_per_s','cold_rejection_W','radiator_area_fraction','sequential_occupancy_mass999']};assert max(selected_errors.values())<1e-12;inputs[str(point_path)]=sha(point_path)
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':inputs,'plant_ledger':{'rows':len(d['rows']),'independent_identities':checks,'maximum_relative_error':max(errors.values()),'category_errors':errors,'selected_point_matching_errors':selected_errors},'dispatch_thermal':dispatch,'qualification':['Hot, precursor, and escaped preparation energy boundaries are kept distinct; preparation photons are debited once.','Large-pulse spacing comes from exact main-window serialization, with observed370day spacing separately reported; a one-day large-shot guard is unnecessary for this corrected bound.','The existing<=1day dispatch-lateness storage allowance exceeds actual audited lateness; no new cadence restriction is introduced.','Tail bounds use exact source128 profiles, worst-case full heat retention and the preserved localized carbon thermal properties. Radiative spectrum and material qualification remain separate.']}
(P/'final-receiver-audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
