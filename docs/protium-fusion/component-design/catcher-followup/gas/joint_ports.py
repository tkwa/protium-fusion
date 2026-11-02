"""Joint local gas-port, active sorber, and optical heat interface.

The molecular transmission is a geometric requirement, not a measured baffle.
All admitted 1000 K radiation is charged to the first cold intercept; a separate
specified leakage fraction reaches 5 K. No global chamber capture is credited.
"""
from pathlib import Path
import json, math, hashlib
P=Path(__file__).resolve().parent
d=json.loads((P/'pumping-interfaces.json').read_text())
i=json.loads((P/'inputs.json').read_text())
k=1.380649e-23; sb=5.670374419e-8; T=1000.; A=4*math.pi*i['R_m']**2
m=[2*i['atomic_mass_kg']['p'],i['atomic_mass_kg']['he4']]
v=[math.sqrt(8*k*T/(math.pi*x)) for x in m]
sactive=[10*math.sqrt(T/300),13*math.sqrt(T/273.15)]
rows=[];checks=[]
for Sp in [.0001,.0005,.001]:
 for mouth in [3e-5,1e-4,3e-4]:
  for transmission in [.03,.1,.3]:
   conductance=[mouth*transmission*x/4 for x in v]
   viable=all(c>Sp for c in conductance)
   r={'downstream_speed_target_m3_s_per_wall_m2_each_species':Sp,
      'shared_mouth_area_fraction':mouth,'net_molecular_transmission':transmission,
      'mouth_conductance_H2_He_m3_s_per_wall_m2':conductance,
      'finite_active_area_possible':viable,
      'hot_to_cold_bypass_upper_W_m2':mouth*sb*T**4}
   if viable:
    active_speed=[1/(1/Sp-1/c) for c in conductance]
    area=[s/sa for s,sa in zip(active_speed,sactive)]
    Q200=area[1]*sb*200**4
    r.update({'required_active_panel_area_fraction_Hgetter_Hecryo':area,
      'required_active_panel_area_m2_Hgetter_Hecryo':[x*A for x in area],
      'active_panel_speed_before_port_H2_He_m3_s_per_wall_m2':active_speed,
      'He_200K_to5K_radiation_W_m2_if_full_200K_view':Q200,
      'He_refrigeration_W_m2_etaCarnot0p5_at200K_sink':Q200*78,
      'hot_photon_leakage_into5K_sensitivity':[{
       'fraction_of_admitted1000K_aperture_energy_reaching5K':f,
       'additional5K_heat_W_m2':f*r['hot_to_cold_bypass_upper_W_m2'],
       'additional_refrigeration_W_m2_etaCarnot0p5':78*f*r['hot_to_cold_bypass_upper_W_m2']}
       for f in [0,1e-5,1e-4,1e-3,.01]],
      'active_hardware_mass_sensitivity_kg_per_wall_m2':[
       {'kg_per_active_panel_m2':s,'Hgetter_Hecryo_total_kg_per_wall_m2':sum(area)*s}
       for s in [1,10,100]],
      'scope':'Shared admission mouth, species-specific active capture areas in parallel, molecular entrance and panel resistances in series. Downstream regeneration and compression remain physically isolated from cavity.'})
    for c,s in zip(conductance,active_speed):checks.append(abs((1/(1/c+1/s))/Sp-1)<1e-12)
   rows.append(r)
candidate=next(r for r in rows if r['downstream_speed_target_m3_s_per_wall_m2_each_species']==.001 and r['shared_mouth_area_fraction']==1e-4 and r['net_molecular_transmission']==.1)
pressure=next(r for r in d['continuous_extraction_scenarios'] if r['downstream_pump_speed_m3_s_per_wall_m2_each_species']==.001)
out={'schema':1,'family':i['family'],'species':['H2','He'],'aperture_model':{
 'hot_gas_temperature_K':1000,'first_intercept_design_temperature_K':300,
 'He_final_radiation_shield_K':200,'He_sorber_temperature_K':5,
 'interpretation':'A locally enclosed gas manifold feeds a shared baffled mouth with no open view to the reactor cavity. The upstream manifold remains hot; all 1000 K radiation crossing its mouth is conservatively deposited in the downstream cold system. A 300 K first intercept/getter and 200 K He shield precede isolated 5 K sorbers. Molecular transmission must be demonstrated for the same baffles; .03/.1/.3 are conditional geometry cases. No ray-traced baffle or loaded-sorber performance is claimed.',
 'conductance_formula':'C_i=a_mouth tau_i vbar_i/4; 1/S_i=1/C_i+1/(a_panel_i S_panel_i).',
 'temperature_normalization':'The selected H speed is 10m3/s/m2 at300K. The Day low-load helium graph scale is13m3/s/m2 at273.15K, normalized to cryopanel area; it does not include this new admission path. Warm-equivalent speeds scale sqrt(1000/Tref) only with maintained capture.',
 'optical_energy_rule':'Qhot_to_cold<=a_mouth sigmaSB 1000^4 for a sealed mouth with no additional warm view. All this heat is debited once from engine heat. Radiation from200K to5K is an internal downstream transfer, not an additional source. Refrigeration work is a separate electrical debit. The specified residual hot photon transmission adds5K refrigeration but does not add a second copy of hot bypass heat.'},
 'scenarios':rows,'representative_interface':candidate,
 'representative_pressure':{key:pressure[key] for key in ['mean_total_source_pressure_Pa','peak_pressure_by_stage_Pa','front_equivalent_open_fraction_for_f_0p001_species','channel_plus_pump_speed']},
 'gas_thermodynamics':d['gas_thermodynamic_scales'],
 'requirements':['The retained surface channels and vias must actually supply the series conductance.','The active sorber must sustain its stated speed at the operating mixture, coverage and base pressure.','A real thermal enclosure must intercept1000K view; the5K stage needs a quantified residual hot transmission and200K view.','Local ports do not remove residual gas already in the reactor cavity. Global cleanup is a separate architecture.','Panel-area mass sensitivities exclude manifolds, refrigerators, controls, storage, shields, service and repair.'],
 'checks':{'series_speed_checks':len(checks),'all_pass':all(checks)},
 'source_sha256':{x:hashlib.sha256((P/x).read_bytes()).hexdigest() for x in ['inputs.json','pumping-interfaces.json']},
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert all(checks)
(P/'joint-ports.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(candidate,indent=2))
