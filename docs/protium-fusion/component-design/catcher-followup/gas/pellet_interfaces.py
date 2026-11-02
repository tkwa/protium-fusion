"""Couple actual gas columns to archived capsule enthalpy and source geometry."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent;R=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');d=json.loads((P/'inputs.json').read_text());ch=json.loads((P/'chamber-dt150.0-q128.json').read_text());cap_path=R/'component-design/capsule-survival/screen.json';cap=json.loads(cap_path.read_text());c=cap['nonmethane']['C13_CH2'];cs=cap['constants'];rho0=c['cold_reference']['rho'];fill0=c['cold_reference']['fill'];u300=c['precooling_20_to300K_enthalpy_bounds_J_kg'];base=ch['cases'][-1];rows=[];eps=8.8541878128e-12
for a in [.01,.02,.05]:
 for fill in [.6,.8,1.]:
  eff=rho0*fill;density_ratio=rho0*fill0/eff
  for q in base['source_profile_pellet_encounters']['c13']:
   v=q['velocity_m_s'];flight=q['flight_s'];IR=c['IR_J_kg_core']*density_ratio*(cs['inner_radius_m']/a)*(flight/cs['flight_time_s']);field=c['electric_J_kg_core']*density_ratio;gas=3*q['maximum_column_kg_m2']*v*v/(8*eff*a);beta=3/(4*eff*a);col=q['maximum_column_kg_m2'];acc={str(E):v*v*eff*a/(6*eps*cs['surface_field_V_m']*E) for E in [1e6,1e7]}
   rows.append({'radius_m':a,'fill_fraction':fill,'core_density_kg_m3':rho0,'core_mass_per_capsule_volume_kg_m3':eff,'mass_velocity_quantile':q['source_mass_quantile'],'v_m_s':v,'flight_s':flight,'full_flight_IR_J_kg_core':IR,'field_energy_J_kg_core':field,'f1g1_prompt_firstflight_gas_J_kg_core':gas,'enthalpy20to300_J_kg_core':u300,'maximum_common_gas_response_multiplier_for_300K':[max(0,(h-IR-field)/gas) for h in u300],'zero_gas_fullflight_reaches300K':[IR+field<=h for h in u300],'gas_scaled_tests':[{'common_gas_response_multiplier':f,'total_specific_energy_J_kg_core':IR+field+f*gas,'maximum_unshielded_path_fraction_for300K': [min(1,max(0,(h-field-f*gas)/IR)) for h in u300]} for f in [0,.01,.1,1]],'approximate20to400K_sensible_allowance_J_kg':[h+1600*(14/15)*100 for h in u300],'additional_full_crystal_melting_latent_J_kg':281000*14/15,'uniform_gas_mean_drag_fraction_Cm1':beta*col,'arrival_delay_uniform_s_Cm1':beta*col*flight/2,'arrival_delay_near_launch_layer_s_Cm1':beta*col*flight,'coherent_transverse_displacement_upper_m_for_3km_s_crossflow':beta*col*3000*flight,'acceleration_length_by_gradient_m':acc,'scope':'Core-specific conservative full incident gas energy; skin/backscatter/re-radiation can alter partition. Multiplier rescales both H/He f1g1 firstflight response, not exactly f/g for other capture times.300K is a selected solid-core comparison; latent/melt values are adjacent polymer screening inputs, not new isotope calorimetry.'})
# Pressure is never the acceptance primitive: link timing and grain tolerances to columns.
limits=[]
for a in [.01,.02,.05]:
 eff=rho0*.8
 for v in [0.31e6,2.4e6,4.15e6]:
  for jitter in [1e-7,1e-6,1e-5]:
   limits.append({'a_m':a,'rho_eff_kg_m3':eff,'v_m_s':v,'differential_arrival_allowance_s':jitter,'column_difference_max_uniform_Cm1_kg_m2':8*eff*a*v*jitter/(3*d['R_m']),'column_difference_max_near_launch_Cm1_kg_m2':4*eff*a*v*jitter/(3*d['R_m']),'scope':'Allowance is a selected comparison, not a measured assembly tolerance. A common modeled drag can be precompensated; spatial/temporal/size fluctuations generate differential error. Sputter recoil multipliesCm separately.'})
# Unresolved directed plume and dust scenarios, all conditional rates.
plumes=[]
A=4*math.pi*d['R_m']**2
for cat,s in d['stages'].items():
 sigma=(s['cold_H_atomic_kg']+s['cold_He_atomic_kg'])/A
 for f in [.001,.01,1]:
  for footprint in [1,1e-3,1e-6]:
   plumes.append({'stage':s['stage'],'cavity_leak_fraction':f,'plume_footprint_fraction':footprint,'mass_column_when_crossed_kg_m2':sigma*f/footprint,'scope':'Uniform column within the specified projected footprint. No bound on unspecified directionality is inferred from chamber-average pressure; timing and path interception must be supplied.'})
dust=[];vp=2.4e6;ap=.01;rp=rho0*.8;mp=4*math.pi/3*rp*ap**3
for rd in [1e-7,1e-6,1e-5,1e-4,1e-3]:
 md=4*math.pi/3*2000*rd**3
 for sig in [1e-12,1e-9,1e-6]:
  N=math.pi*(ap+rd)**2*sig/md;dust.append({'dust_radius_m':rd,'dust_density_kg_m3':2000,'dust_column_kg_m2':sig,'expected_interceptions_per_pellet':N,'probability_at_least_one_Poisson':-math.expm1(-N),'single_impact_energy_per_pellet_mass_J_kg':.5*md*vp*vp/mp,'scope':'Illustrative independent uniform dust grains with small dust velocity; no sputter/dust production rate is assigned. The rare-impact distribution differs from average molecular heating.'})
radial_plumes=[]
Qgas=d['plant']['recovered_H_kg_s']+d['plant']['helium_per_s']*d['atomic_mass_kg']['he4']
for f in [.0001,.001,.01]:
 for theta in [.001,.01,.1]:
  rmin=d['R_m']*theta;u=3000;col=Qgas*f/(4*math.pi*u)*(1/rmin-1/d['R_m'])
  radial_plumes.append({'leaked_mass_fraction':f,'effective_radial_beam_halfangle_rad':theta,'minimum_focusing_radius_m':rmin,'gas_bulk_speed_m_s':u,'radial_column_kg_m2':col,'scope':'Conditional steady inward radial gas flux with finite angular blur, not a measured pore plume. Gas collisions or pressure can broaden or shock the central flow; uniform chamber pressure does not bound this concentration.'})
out={'inward_radial_plume_scenarios':radial_plumes,'source_files_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'inputs.json',P/'chamber-dt150.0-q128.json',cap_path]},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'capsule_reference':c,'capsule_constants':cs,'C13_cases':rows,'differential_timing_column_screens':limits,'directed_plume_scenarios':plumes,'dust_scenarios':dust,'notes':['At fixed surfacefield and acceleratinggradient, increasingradius increases required acceleratorlength linearly; it reducescharge/mass andpacketcount. No sourceperformance credit is taken without updatinghardware.','Ionization/heating does not removeH orHe from the chamber. Promptbackscatteredenergetic ions andsecondaryelectrons are a separate source/charge/energy population, not instantly1000K neutralgas.','Both radiation andgas energy varyacross actualincomingvelocityprofile. Slowgrains canfail300Kfullflight criterion evenwithoutgas; protecteddelivery orcontainedmeltingarephysicaldesignchoices.','Fullgasenergy into coreis anupperbulkenergyallocation. Localstoppingdepth/skinablationmustbechecked separately; a bulkthermalpass doesnotqualify a coating.']};(P/'pellet-interfaces.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([r for r in rows if r['fill_fraction']==.8 and r['mass_velocity_quantile'] in [.1,.5,.999]],indent=2))
