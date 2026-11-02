#!/usr/bin/env python3
"""Bounded cold-sink heat-engine/HX/insulation requirements; no station optimization.
Run python3 heat_engine.py. Standard library only; reads frozen station inputs.
"""
from pathlib import Path
import hashlib,json,math
P=Path(__file__).resolve().parent
OLD=P.parent/'open-loop-redesign'
SB=5.670374419e-8
SOURCES={
'NASA_Brayton_model_and_mass':'https://ntrs.nasa.gov/api/citations/20050244475/downloads/20050244475.pdf',
'NASA_conversion_review':'https://ntrs.nasa.gov/api/citations/20200000369/downloads/20200000369.pdf',
'NASA_MLI_scope':'https://www.nasa.gov/smallsat-institute/sst-soa/thermal-control/',
'NASA_HX_design_project':'https://techport.nasa.gov/projects/113159'}
def sha(x):return hashlib.sha256(Path(x).read_bytes()).hexdigest()
def envelope(eta,tr,phi=.8,ts=1000.,q=106.61795837196411):
    r=1-eta/phi
    out={'heat_efficiency':eta,'radiator_K':tr,'hot_source_K':ts,'internal_Carnot_fraction_phi':phi,'recovered_heat_flux_W_m2':q,'external_Carnot_limit':1-tr/ts}
    if r*ts-tr<=1e-9:
        return out|{'feasible':False,'reason':'No positive finite exchanger approaches for this phi/eta/temperature combination.'}
    s=math.sqrt((1-eta)*r);th=(tr+s*ts)/(r+s);tc=r*th
    uh=q/(ts-th);uc=(1-eta)*q/(tc-tr)
    residual=eta-phi*(1-tc/th)
    return out|{'feasible':True,'engine_effective_hot_K':th,'engine_effective_cold_K':tc,'hot_approach_K':ts-th,'cold_approach_K':tc-tr,'UA_hot_W_K_m2':uh,'UA_cold_W_K_m2':uc,'UA_external_total_W_K_m2':uh+uc,'UA_per_input_W_per_K':(uh+uc)/q,'dedicated_external_HX_kg_m2_by_conductance_W_K_kg':{str(k):(uh+uc)/k for k in [1,10,100]},'entropy_generation_W_K_m2':q*((1-eta)/tr-1/ts),'efficiency_identity_residual':residual,'scope':'Minimum sum of conductances for isothermal effective engine ports at specified phi. Recuperator, pressure losses, hardware mass, coolant pump and non-isothermal source distribution are additional.'}

def brayton(pr, *, nc=6,nt=6,th=950.,tc=230.,ec=.93,et=.93,er=.95,eg=.98,dp=.03,ts=1000.,tr=200.,q=106.61795837196411,pump_fraction=.01,other_loss_fraction=.01):
    # Constant-cp monatomic gas, gamma5/3; full intercooling and reheat with finite approach.
    # Total flow-path pressure-loss fraction reduces turbine expansion ratio. No free bleed.
    re=pr*(1-dp)
    if re<=1:return None
    top=tc*(1+(pr**(.4/nc)-1)/ec);bot=th*(1-et*(1-re**(-.4/nt)))
    if not tc<top<th or not tc<bot<th:return None
    rec=er*max(0,bot-top);pre=top+rec;coldout=bot-rec
    wc=nc*(top-tc);wt=nt*(th-bot);qh=th-pre+(nt-1)*(th-bot);qc=coldout-tc+(nc-1)*(top-tc)
    shaft=wt-wc
    if qh<=0 or shaft<=0:return None
    eff=eg*shaft/qh
    if eff>=1-tr/ts:return None
    # UA from exact isothermal-reservoir LMTD on each finite-temperature heat exchange.
    cdot=q/qh
    ua_hot=cdot*(math.log((ts-pre)/(ts-th))+(nt-1)*math.log((ts-bot)/(ts-th)))
    ua_cold=cdot*(math.log((coldout-tr)/(tc-tr))+(nc-1)*math.log((top-tr)/(tc-tr)))
    ua_rec=cdot*er/(1-er) if rec>0 else 0
    ua_gen=cdot*(1-eg)*shaft/(tc-tr)
    eta_net=eff-pump_fraction-other_loss_fraction
    uat=ua_hot+ua_cold+ua_rec+ua_gen
    return {'pressure_ratio':pr,'compressor_groups':nc,'turbine_groups':nt,'stage_compressor_isentropic_efficiency':ec,'stage_turbine_isentropic_efficiency':et,'recuperator_effectiveness':er,'generator_efficiency':eg,'total_pressure_loss_fraction':dp,'turbine_inlet_K':th,'compressor_inlet_K':tc,'compressor_stage_outlet_K':top,'turbine_stage_outlet_K':bot,'heater_inlet_after_recuperator_K':pre,'cooler_inlet_after_recuperator_K':coldout,'thermal_to_generator_efficiency':eff,'separate_coolant_pump_fraction_Qhot':pump_fraction,'separate_bearing_windage_misc_fraction_Qhot':other_loss_fraction,'net_converter_efficiency_after_declared_losses':eta_net,'heat_input_per_cp_K':qh,'recuperated_heat_over_external_heat':rec/qh,'external_hot_source_K':ts,'radiator_K':tr,'UA_hot_W_K_m2':ua_hot,'UA_cold_W_K_m2':ua_cold,'UA_recuperator_W_K_m2':ua_rec,'UA_generator_cooling_W_K_m2':ua_gen,'UA_total_W_K_m2':uat,'dedicated_HX_kg_m2_by_conductance_W_K_kg':{str(k):uat/k for k in [1,10,100]},'engine_mass_kg_m2_at10We_kg_before_external_pump':q*eff/10,'thermal_capacity_flow_W_K_m2':cdot,'working_fluid_kg_s_m2_if_cp1000':cdot/1000,'recuperator_NTU_equal_capacity':er/(1-er),'cycle_energy_residual_per_heat_input':(qh-qc-shaft)/qh,'external_entropy_generation_W_K_m2':q*((1-eta_net)/tr-1/ts),'scope':'Conditional constant-cp cycle. Component maps, actual fluid mixture/properties, seals, recuperator axial conduction, bypass/leakage, cooler pressure loss allocation and lifetime need design. Pressure losses are included once via reduced expansion ratio. No whole-station net-power result.'}

def best_cycle(**kwargs):
    good=[brayton(1.05*(80/1.05)**(j/2500),**kwargs) for j in range(2501)]
    return max((x for x in good if x),key=lambda x:x['net_converter_efficiency_after_declared_losses'])

def main():
    files=['station.py','operating-cases.json','station-inputs.json'];base=json.loads((OLD/'operating-cases.json').read_text())['cases']['high_source_rating_surface'];inputs=json.loads((OLD/'station-inputs.json').read_text());q=base['mean_recovered_heat_flux_W_m2'];area=4*math.pi*base['radius_m']**2
    rows=[envelope(e,t,p,q=q) for e in [.4,.5,.6,.65,.7] for t in [197.,200.,225.,250.,300.,450.] for p in [.7,.8,.85,.9]]
    cases={}
    for tr in [197.,200.,250.,300.]:
        for name,n,ec,et,rec,gen,dp in [('three_stage_93pct',3,.93,.93,.95,.98,.03),('six_stage_93pct',6,.93,.93,.95,.98,.03),('six_stage_97pct',6,.97,.97,.98,.99,.02)]:
            cases[f'{name}_radiator{tr:g}K']=best_cycle(nc=n,nt=n,ec=ec,et=et,er=rec,eg=gen,dp=dp,tc=tr+30,tr=tr,q=q)
    leak=[]
    for tr in [197.,200.,250.,300.]:
        black=SB*(1000**4-tr**4)
        for loss in [0.,1.,10.]:
            leak.append({'radiator_K':tr,'bypass_W_m2':loss,'raw_black_leak_W_m2':black,'effective_emissivity_ceiling_if_no_bridges':loss/black,'heat_to_engine_W_m2':q-loss,'converter_gross_at60pct_W_m2':.6*(q-loss),'cold_heat_from_engine_plus_bypass_W_m2':.4*(q-loss)+loss,'gross_penalty_W_m2':.6*loss,'energy_residual_W_m2':.6*(q-loss)+.4*(q-loss)+loss-q})
    shields=[]
    for ep in [.003,.03,.1]:
        for n in [100,900,1000,3000]:
            ee=1/((n+1)*(2/ep-1));black=SB*(1000**4-200**4)
            shields.append({'individual_gray_emissivity_requirement':ep,'floating_shield_count':n,'effective_emissivity_ideal':ee,'radiative_bypass_W_m2':ee*black,'foil_kg_m2_if1micron_iron_density7874':n*1e-6*7874,'foil_kg_m2_if1micron_generic_density3000_to8000':[n*1e-6*3000,n*1e-6*8000],'foil_kg_m2_if1micron_tungsten_density19300_NOT_selected_material':n*1e-6*19300,'scope':'Equal-emissivity diffuse-gray vacuum gaps including matching boundaries; no contact bridges, seams, holes, gas conduction, distortion or lifetime. Refractory foil density is a mass scenario, not validated emissivity or film construction.'})
    bridges=[]
    for k in [1.,10.,100.]:
        for f in [1e-6,1e-4,.001]:
            bridges.append({'bridge_k_W_m_K':k,'solid_area_fraction':f,'length_m':1.,'conductive_bypass_W_m2':k*f*800})
    pulses={k: v['received_energy_J']/area for k,v in inputs['receiver']['stages'].items()};capacity=base['absorber_sigma_kg_m2']*1000
    hx_ratings=[]
    # Simple compact plate core geometry supports interpreting the rating sensitivity.
    for h in [100,1000,10000]:
        for wall in [1e-4,1e-3]:
            rho=8000.;kw=15.;U=1/(2/h+wall/kw);arealm=rho*wall*5
            hx_ratings.append({'film_h_each_side_W_m2_K':h,'wall_thickness_m':wall,'wall_conductivity_W_m_K':kw,'density_kg_m3':rho,'fivefold_mass_multiplier_for_fins_headers_case':5.,'overall_U_W_m2_K':U,'specific_conductance_W_K_kg':U/arealm,'scope':'Geometry requirement with equal fluid films and explicit5xmass. Pressure drop, fouling, joint stresses and distribution remain uncomputed; does not establish an attainable h at chosen pump power.'})
    checks=[abs(x['efficiency_identity_residual'])<1e-12 and x['entropy_generation_W_K_m2']>=0 for x in rows if x['feasible']]
    checks +=[abs(x['cycle_energy_residual_per_heat_input'])<1e-12 and x['external_entropy_generation_W_K_m2']>0 for x in cases.values()]
    checks +=[abs(x['energy_residual_W_m2'])<1e-12 for x in leak]
    checks +=[not envelope(.6,250.,.8)['feasible'],not envelope(.65,197.,.8)['feasible'],envelope(.65,200.,.85)['feasible'],envelope(.7,200.,.9)['feasible']]
    checks +=[x['heat_efficiency']<x['external_Carnot_limit'] for x in rows if x['feasible']]
    out={'scope':'Bounded conditional heat-engine, finite exchanger and hot/cold insulation requirements. Sunlight set absent by user premise; external Dyson thermal background belongs in separate station environment model. No station inputs or article modified.','source_sha256':{f:sha(OLD/f) for f in files},'code_sha256':sha(__file__),'primary_sources':SOURCES,'baseline':{'radius_m':base['radius_m'],'surface_area_m2':area,'captured_heat_W_m2':q,'old_heat_efficiency':.4,'old_radiator_K':450.,'old_engine_equipment_rating_We_kg':10.,'old_insulation_allocation_kg_m2':100.,'old_no_dedicated_HX_mass':True},'finite_approach_envelopes':rows,'conditional_regenerative_Brayton_cycles':cases,'bypass_accounting':leak,'radiation_shield_requirements':shields,'bridge_requirements':bridges,'HX_geometry_rating_sensitivity':hx_ratings,'pulse_buffer':{'stage_received_energy_J_m2':pulses,'absorber_heat_capacity_J_m2_K':capacity,'each_stage_single_pulse_bulk_jump_K':{k:v/capacity for k,v in pulses.items()},'independent_periodic_stage_sawtooth_sum_range_bound_K':sum(pulses.values())/capacity,'scope':'Ideal spatially equilibrated periodic reservoirs with constant mean extraction; sum of each stage energy bounds worst combined sawtooth range. Actual surface conduction, beta-delay heat, phase timing and heat-pipe pickup need transient design. Fixed1000K source is an approximation over this excursion.'},'checks':{'count':len(checks),'passed':sum(checks),'all_passed':all(checks)}}
    (P/'heat-engine.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':out['checks'],'cycles':{k:{f:v[f] for f in ['thermal_to_generator_efficiency','net_converter_efficiency_after_declared_losses','UA_total_W_K_m2','recuperator_NTU_equal_capacity']} for k,v in cases.items()},'buffer':out['pulse_buffer']},indent=2))
if __name__=='__main__':main()
