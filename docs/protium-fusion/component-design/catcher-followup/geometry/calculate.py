"""Catcher geometry, thermal-resistance and recycled-throughput requirements.
No assumed effective irradiated barrier permeability, repair speed or capture.
"""
from pathlib import Path
import csv,json,math,hashlib
P=Path(__file__).resolve().parent;C=P.parents[1];F=P.parent
read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sp=read(F/'stopping/surface-spectra.json');stop=read(F/'stopping/stopping-results.json');dep=read(F/'stopping/surface-stopping.json');charged=read(F/'stopping/surface-charged-flux.json');hist=read(P/'temperature-histories.json');inp=read(C/'final-receiver-thermal-input.json');point=read(C.parent/'mixture-design/selected-plant/point.json')['plant']
rho=2000.;R=1.1e9;A=4*math.pi*R*R;year=365.25*86400.;layer=150e-6;cap=100e-9;width=100e-9;pitch=20e-6;slotf=width/pitch;gap=50e-6;postf=.25;postk=50.;skin_k=4.;via_r=.5e-3;via_pitch=.01;via_f=math.pi*via_r**2/via_pitch**2;L=.175
# Vias remove backing solid; adding physical thickness retains the full carbon column.
added_thickness=(L+gap*(1-postf)+layer*slotf)/(1-via_f)-L
cat={'A':'n14','B1':'n15','B2':'c12','C13':'c13'};stages=[]
for n,s in sp['stages'].items():
 row=next(x for x in stop['rows'] if x['stage']==n);nom=next(x for x in row['rows'] if x['electronic_stopping_multiplier']==1)
 groups={}
 for name,species in [('H','p'),('He','he4')]:
  d=dep['stages'][n]['species'][species];weights=d['packet_nuclei_m2'];z=d['packet_mean_projected_depth_um'];tot=sum(weights)
  groups[name]={'number_m2':tot,'fraction_beyond_depth_mean_projected_proxy':{str(u):sum(w for w,x in zip(weights,z) if x>u)/tot for u in [50,100,150]},'rear_before_front_absorbing_slab_150um':sum(w*min(x/150,1) for w,x in zip(weights,z))/tot}
 r=inp['families']['source128']['stage_inputs'][cat[n]]['selected']['receiver'];scale=(r['radius_m']/R)**2
 prompt_peak=r['prompt']['radial_ray_arrival_surface_screen']['peak_flux_W_m2']*scale
 peakN=charged['stages'][n]['quadrature64_peak_sampled_number_flux_m2_s']
 stages.append({'stage':n,'period_s':s['period_s'],'pulses_year':year/s['period_s'],'incident_mass_kg_m2':s['total_mass_fluence_kg_m2_baryon'],'groups':groups,'maximum_primary_CSDA_depth_um_by_stopping_multiplier':{str(x['electronic_stopping_multiplier']):x['maximum_endpoint_velocity_CSDA_column_kg_m2']/rho*1e6 for x in row['rows']},'all_charged_sampled_peak_number_m2_s':peakN,'prompt_radial_peak_W_m2':prompt_peak,'cap_lateral_rise_K_if_all_incident_flux_absorbed_k4':prompt_peak*width**2/(8*4*cap),'NRT_proxy_events_m2_Ed25':nom['all_charged_NRT_proxy_events_m2']['25.0'],'NRT_proxy_cap0to100nm_dpa_Ed25':sum(b['NRT_proxy_events_m2']['25.0'] for b in nom['elastic_transfer_depth_bins'] if b['depth_um'][1]<=.1)/(rho*.1e-6/(12.011*1.66053906660e-27))})
profiles=[]
for p in hist['profiles']:
 rows=list(csv.DictReader(Path(p['csv']).open()));q=[float(r['material_W_m2'])+float(r['delayed_decay_W_m2']) for r in rows];base=[float(r['front_bulk_k50_K']) for r in rows]
 # Keep a conservative full150um/k4 resistance above the backing solution; this
 # does not subtract the intact backing layer's resistance already represented.
 resistance=layer/skin_k+gap/(postk*postf)
 thermal=[t+f*resistance for t,f in zip(base,q)]
 profiles.append({'stage':p['stage'],'selection':p['selection'],'maximum_K_backing':max(base),'maximum_K_with150um_k4_and_postgap':max(thermal),'peak_material_plus_decay_W_m2':max(q),'max_added_K':max(q)*resistance,'csv_sha256':sha(p['csv'])})
rate=sum(1/s['period_s'] for s in sp['stages'].values());net_flux=point['net_W']/A if 'net_W' in point else point['net_electric_W']/A
renewal=[]
for depth in [0.1,1,10,50,100,150]:
 mass=rho*depth*1e-6
 for mode,nperyear in [('annual',1.),('monthly',12.),('each_pulse',rate*year)]:
  renewal.append({'depth_um':depth,'schedule':mode,'handled_kg_m2_year':mass*nperyear,'mean_handled_kg_m2_s':mass*nperyear/year,'electric_W_m2_at10kJkg_handling':mass*nperyear/year*1e4,'electric_W_m2_at2MJkg':mass*nperyear/year*2e6,'fraction_selected_net_power_at2MJkg':mass*nperyear/year*2e6/net_flux})
# Approximate hydraulic-diameter long-slit molecular-flow screen, not an exact
# Clausing solution. Gas agent owns exact throughput/pressure/compression.
beta_slot=slotf*8*width/(3*layer);beta_via=via_f*8*via_r/(3*L);beta_series=1/(1/beta_slot+1/beta_via)
gas=[]
for species,mass_u in [('H2',2.),('He',4.)]:
 speed=math.sqrt(8*1.380649e-23*1000/(math.pi*mass_u*1.66053906660e-27));gas.append({'species':species,'mean_speed_1000K_m_s':speed,'slot_C_over_area_m_s':speed/4*beta_slot,'via_C_over_area_m_s':speed/4*beta_via,'series_C_over_area_m_s':speed/4*beta_series})
open_cells=[]
for aspect in [10,100,1000,2000]:
 K=min(1.,4/(3*aspect))
 for s in [1.,.1,.01]:open_cells.append({'length_over_diameter':aspect,'single_return_molecular_transmission':K,'return_capture_probability':s,'eventual_front_escape':K/(K+s*(1-K))})
holes=[]
for stage in stages:
 for effective_radius_nm in [.1,1.,10.]:
  ah=math.pi*(effective_radius_nm*1e-9)**2
  for global_open_fraction in [1e-8,1e-6,1e-4]:holes.append({'stage':stage['stage'],'effective_connected_hole_radius_nm':effective_radius_nm,'global_connected_open_fraction':global_open_fraction,'required_mean_open_lifetime_s_upper':global_open_fraction/(slotf*stage['all_charged_sampled_peak_number_m2_s']*ah),'assumption':'Only hits above slit mouths create through-connected apertures; effective hole area includes probability of such damage.64-point sampled incident peak is a resolved-profile estimate, not an unresolved-structure upper bound. Bulk permeation and seals are separate.'})
checks={'all_nominal_primary_ranges_below150um':all(s['maximum_primary_CSDA_depth_um_by_stopping_multiplier']['1.0']<150 for s in stages),'mass_preserved_by_added_thickness':abs(rho*((L+added_thickness)*(1-via_f)-gap*(1-postf)-layer*slotf)-350)<1e-9,'thermal_resistance_positive':resistance>0,'history_minimum_positive':all(p['maximum_K_backing']>1000 for p in profiles),'flow_series_smaller_than_components':all(g['series_C_over_area_m_s']<min(g['slot_C_over_area_m_s'],g['via_C_over_area_m_s']) for g in gas)}
assert all(checks.values()),checks
out={'scope':__doc__,'inputs_sha256':{str(p):sha(p) for p in [P/'temperature-histories.json',F/'stopping/surface-spectra.json',F/'stopping/surface-stopping.json',F/'stopping/surface-charged-flux.json',F/'stopping/stopping-results.json',C/'final-receiver-thermal-input.json',C.parent/'mixture-design/selected-plant/point.json']},'code_sha256':sha(__file__),'geometry':{'rho_kg_m3':rho,'active_depth_m':layer,'cap_thickness_m':cap,'slot_full_width_m':width,'slot_pitch_m':pitch,'slot_projected_fraction':slotf,'plenum_height_m':gap,'thermal_post_area_fraction':postf,'post_conductivity_W_m_K':postk,'gas_via_radius_m':via_r,'gas_via_pitch_m':via_pitch,'gas_via_fraction':via_f,'retained_graphite_sigma_kg_m2':350.,'additional_physical_depth_m_preserving_mass':added_thickness,'active_solid_kg_m2':rho*layer*(1-slotf),'dogleg_stopping_pad_kg_m2_at_two_via_areas_200um':rho*200e-6*2*via_f,'cap_full_area_kg_m2':rho*cap,'plenum_post_kg_m2':rho*gap*postf,'gas_void_column_kg_m2_equivalent':rho*added_thickness,'new_mass_scope':'Active carbon, posts and vias substitute inside the350kg/m² stock. Added gas seals, compressors, valves and rework tools remain service hardware; no zero hardware-mass credit is claimed.','skin_and_post_extra_resistance_m2_K_W':resistance,'cap_diffusion_time_width2_over_alpha_k4_s':width**2/(4/(rho*1800))},'stages':stages,'thermal_profiles':profiles,'routine_full_layer_rework_sensitivities':renewal,'molecular_conductance_estimate':{'beta_slot':beta_slot,'beta_via':beta_via,'beta_series':beta_series,'rows':gas,'scope':'Diffuse molecular flow; hydraulic-diameter slit approximation plus long circular via. Side plenum, bends, surface accommodation, entrance distribution and pump speed add resistance. No pressure or collection efficiency is assigned.'},'open_cell_sensitivity':open_cells,'connected_hole_acceptance_curve':holes,'net_electric_flux_W_m2':net_flux,'checks':checks,'limits':['NIST projected depths are mean-range proxies, not stochastic retention distributions. Range straggling, defects and channel-aligned transmission remain qualifications.','The0.5% slits require a second stopping floor in the sealed plenum; full backing retains energetic particles, while selective extraction from that floor is separate.','Finite100nm cap/sidechannel first-passage probabilities are supplied by the parent channel-routing packet; geometry has no assumed isotropic retention credit.','The source-energy-conserving temperature histories omit old nonuniform tails bounded separately at about0.1K. Microsecond prompt limits remain separate from material-pulse thermal resistance.','A2MJ/kg dissipation is a declared rework scenario, not a universal thermodynamic cost. Recycled matter remains installed stock; energy and equipment throughput remain chargeable.','NRT counts are total-elastic-energy proxies without recoil partition, recombination, annealing or a mapping to through-connected holes. No probabilistic material-property interval can presently be justified.']}
(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'geometry':out['geometry'],'thermal':profiles,'gas':gas,'checks':checks},indent=2))
