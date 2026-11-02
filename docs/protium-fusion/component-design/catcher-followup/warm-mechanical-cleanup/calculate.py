"""Bounded warm mechanical chamber-cleanup requirements, not a qualified pump.
Warm baffle conductance + molecular pump in series; hot rotor compresses gas
before a small thermally isolated downstream aperture. Separate species cases.
"""
from pathlib import Path
import math,json,hashlib
P=Path(__file__).resolve().parent;C=P.parents[1];kb=1.380649e-23;mu=1.66053906660e-27;sb=5.670374419e-8
read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
b=read(P/'benchmarks/results.json');model=b['benchmarks'][0];T=1000.;Tc=300.;Pin=1e-10;Pf=1e-3;Pc=1e-4;g=.01;K=1e9;u=700.;rho=1800.;radius=1.;f_in=.1;tau=.3;massA=model['mass_kg']/(math.pi*model['nominal_inlet_diameter_m']**2/4)
alpha={s:next(x['effective_aperture_capture_if_same_speed'] for x in b['calculations'] if x.get('model')==model['model'] and x.get('species')==s and x.get('T_K')==300.) for s in ['H2','He']};alpha['H']=alpha['H2']/math.sqrt(2)
mass={'H':1.00784,'H2':2.01588,'He':4.002602};v=lambda s,t:math.sqrt(8*kb*t/(math.pi*mass[s]*mu))
flow_margin=1-Pf/(K*Pin);g0=g/flow_margin;beta_b=f_in*tau;beta_p=1/(1/g0-1/beta_b);fp=beta_p/alpha['H2'];point={'target_g_H2':g,'species_Pin_Pa':Pin,'species_hot_forepressure_Pa':Pf,'chosen_zero_flow_K0_H2':K,'finite_flow_speed_fraction':flow_margin,'hot_entrance_fraction':f_in,'baffle_molecular_transmission':tau,'effective_baffle_beta':beta_b,'required_pump_aperture_fraction':fp,'pump_capture_H2_requirement':alpha['H2'],'pump_capture_He_comparison':alpha['He'],'temperature_K':T,'tip_speed_m_s':u,'rotor_density_kg_m3':rho,'rotor_radius_m':radius,'commercial_whole_pump_mass_per_nominal_aperture_kg_m2':massA,'package_mass_at_benchmark_areal_density_kg_m2_receiver':fp*massA,'extra_baffle_column_if350kg_per_mouth_m2':f_in*350,'scope':'H2 sets this reference geometry. Preserving cold measured light-gas capture at1000K is an engineering requirement. Baffle and pump are connected by a warm molecular-flow volume; conductances combine in series. g is the isotropic-wall-encounter equivalent, not a directional plume guarantee.'}
rows=[]
for s in mass:
 beta=1/(1/(f_in*tau)+1/(fp*alpha[s]));N=g*v(s,T)/4*Pin/(kb*T)
 # Ideal warm outlet -> cooler pump aperture with correct thermal-transpiration
 # reverse flux; pressures on each side are independently specified.
 therm_ratio=Pc/Pf*math.sqrt(T/Tc);aport=N/(tau*v(s,T)/(4*kb*T)*Pf*(1-therm_ratio));S_cold=N*kb*Tc/Pc
 rows.append({'species':s,'assumed_capture_alpha':alpha[s],'mean_speed_1000K_m_s':v(s,T),'g_if_K1e9_for_this_species':beta*flow_margin,'required_gross_speed_per_receiver_area_m_s':g*v(s,T)/4,'number_rate_m2_s_at_target_g':N,'mass_rate_kg_m2_s':N*mass[s]*mu,'ideal_isothermal_compression_W_m2':N*kb*T*math.log(Pf/Pin),'electric_compression_W_m2_at50percent':2*N*kb*T*math.log(Pf/Pin),'cold_pump_speed_m3_s_per_receiver_m2_at300K_P1e4minus':S_cold,'warm_to_cold_aperture_fraction_at_tau0p3':aport,'all_hot_aperture_radiation_to_cold_W_m2_upper':aport*sb*(T**4-Tc**4),'thermal_reverse_flux_fraction':therm_ratio,'scope':'Separate single-species case; Pin andPf are partialpressures. H alpha is an explicit same-S0 light-gas extrapolation, not a measurement. AtomicH also needs its own compression stages or demonstrated recombination.'})
network=[]
for f in [.03,.05,.1,.2,.3]:
 for tr in [.1,.3,.5]:
  C=f*tr
  for af in [.5,1.,2.]:
   possible=C>g0
   ap=1/(1/g0-1/C)/(alpha['H2']*af) if possible else None
   network.append({'entrance_fraction':f,'baffle_transmission':tr,'capture_relative_to_H2_benchmark':af,'finite_pump_possible':possible,'required_pump_aperture_fraction':ap,'whole_pump_mass_at_benchmark_density_kg_m2':None if ap is None else ap*massA})
compression=[]
for s in ['H2','He','H']:
 Kref=model['K0'].get(s,model['K0']['H2']**(1/math.sqrt(2)))
 for speed in [500.,700.,1000.,1500.,2000.]:
  scale=speed/400*math.sqrt(300/T);Ksame=math.exp(math.log(Kref)*scale)
  compression.append({'species':s,'tip_m_s':speed,'reference_tip_assumption_m_s':400.,'same_geometry_allhot_similarity_K0':Ksame,'log_compression_length_factor_for_K1e8':math.log(1e8)/(math.log(Kref)*scale),'log_compression_length_factor_for_K1e9':math.log(1e9)/(math.log(Kref)*scale),'scope':'Similarity of lnK withu/sqrtT. Reference400m/s is a source-scale assumption, not measured tip speed of this exactmodel. Added effective compression length neednotscale wholepackage mass.'})
pressure=[]
for pf in [1e-4,1e-3,.01,1.,4.]:
 for p in [1e-10,1e-11]:
  pressure.append({'Pin_species_Pa':p,'Pfore_species_Pa':pf,'K0_for90percent_speed':pf/(.1*p),'K0_for99percent_speed':pf/(.01*p),'floor_with_K1e9_Pa':pf/1e9})
channel=[]
for s in mass:
 for speed in [500,700,1000,2000]:
  for gg in [1.,2.,3.]:channel.append({'species':s,'tip_m_s':speed,'geometry_factor':gg,'ideal_L_over_h_for_K1e9':math.log(1e9)*v(s,T)*gg/speed,'scope':'Manufacturer idealHolweck exponential channel model. Inter-groove/land backflow, accommodation and tolerances lower achievedK; not a selectedphysical channel.'})
stress=[]
for density in [1500.,1800.,2000.]:
 for speed in [500.,700.,1000.,1500.,2000.]:stress.append({'density_kg_m3':density,'tip_m_s':speed,'thin_ring_hoop_Pa':density*speed**2,'uniform_isotropic_solid_disk_max_Pa_nu0p2':.4*density*speed**2,'scope':'Geometry reference stresses; an anisotropicbladedrotor, hubhole, joins and margins need their ownanalysis.'})
loss=[]
for rad in [.1,1.,10.]:
 for rotor_mA in [10.,100.,300.]:
  E=fp*rotor_mA*.25*u*u;omega=u/rad;nunits=fp/(math.pi*rad*rad)
  for qloss in [.1,1.,5.]:loss.append({'radius_m':rad,'rotor_mass_per_pump_aperture_kg_m2_requirement':rotor_mA,'loss_budget_W_m2_receiver':qloss,'allowed_loss_W_per_pump':qloss/nunits,'allowed_loss_W_per_pump_aperture_m2':qloss/fp,'stored_rotor_J_m2_receiver_disk_model':E,'free_coast_energy_timescale_days':E/qloss/86400,'equivalent_rotational_Q_omegaE_overP':omega*E/qloss,'net_electric_penalty_W_m2_if_fraction0_or1_recovered_at60percent':[qloss,.4*qloss],'scope':'Specified loss and arealmass requirements, not a prediction. Radius/count are varied without asserting same package mass for geometricallyscaled commercial hardware.'})
thermal=[]
for qbudget in [.1,1.,5.]:
 thermal.append({'warm_to_cold_heat_budget_W_m2_receiver':qbudget,'allowed_effective_cold_view_emissivity_times_viewfactor_per_pump_area':qbudget/(fp*sb*(1000**4-300**4)),'allowed_G_W_K_per_pump_aperture_m2':qbudget/(fp*(1000-300)),'support_area_fraction_per_pump_aperture_at_k10_L0p1m':qbudget*.1/(fp*10*(1000-300))})
wind=[]
for s in mass:
 for p in [1e-10,1e-3]:
  shear_power=p/(kb*T)*mass[s]*mu*v(s,T)/4*u*u
  wind.append({'species':s,'local_pressure_Pa':p,'fully_accommodating_shear_power_W_per_moving_surface_m2':shear_power,'scope':'One-sidedincidentMaxwellian at rest, full accommodation, wall speedu; an upper-order local drag estimate for this definedflow, nottotalpumpheat or an intrinsicbearingloss.'})
maintenance=[]
for frac_year in [1e-4,.01,1.]:maintenance.append({'fraction_of_benchmark_package_reworked_per_year':frac_year,'electric_W_m2_at2MJkg':point['package_mass_at_benchmark_areal_density_kg_m2_receiver']*frac_year*2e6/(365.25*86400),'scope':'Conditional same-materialrework dissipation. No fullrotorreplacement rate is inferred; moving or localrepair maycostless.'})
# Concrete downstream gas-only conduit example:larger aperture compensates
# for a long baffled tube; all incoming hotradiation counted once at coldintercept.
unitA=math.pi*radius**2;nunits=fp/unitA;rt=.001;Lt=.1;tt=1e-5;kpipe=10.;bend=.5
fpipe=nunits*math.pi*rt**2;tautube=8*rt/(3*Lt)*bend;required_beta=rows[1]['warm_to_cold_aperture_fraction_at_tau0p3']*.3
conduit={'unit_radius_m':radius,'units_per_receiver_m2':nunits,'pipe_radius_m':rt,'pipe_length_m':Lt,'pipe_wall_m':tt,'pipe_conductivity_W_m_K_scenario':kpipe,'bend_transmission_multiplier_scenario':bend,'diffuse_longtube_transmission_estimate':tautube,'pipe_projected_area_per_receiver_area':fpipe,'net_aperture_beta':fpipe*tautube,'required_beta_at_target_flow':required_beta,'conductance_margin':fpipe*tautube/required_beta,'all_admitted_hot_radiation_to_cold_W_m2_upper':fpipe*sb*(T**4-Tc**4),'pipe_wall_conduction_W_m2':nunits*2*math.pi*rt*tt*kpipe/Lt*(T-Tc),'scope':'Gas-only conduit to300Kintercept; a molecular-flow longtube+bend estimate, not a ray-traced model. Warm support/bearing heatpaths are separate. First300Kintercept blocks broad5Kview.'}
# Independent measured ordinary-temperature pump input, from companion's
# primarypaper review. Size scaling is an acceptance test, not an apparatusclaim.
idle={'primary_measured_module_input_W':48.,'primary_disk_diameter_m':.26,'primary_measured_input_W_per_disk_m2':48/(math.pi*.13**2),'chosen_receiver_idle_budget_W_m2':1.,'required_idle_W_per_pump_aperture_m2':1/fp,'ratio_measured_areal_to_required':48/(math.pi*.13**2)*fp,'module_radius_m_if48W_per_module_held_fixed':math.sqrt(48*fp/math.pi),'scope':'Ordinary-temperature1994 measured completeinput; fixed48W enlargedmodule calculation doesnotpredict itsmass,hotbearings orloss.'}
checks=[]
for s in ['H2','He']:
 S0=fp*alpha[s];C0=f_in*tau;Pmid=(C0*Pin+S0*Pf/K)/(C0+S0);Q1=C0*(Pin-Pmid);Q2=S0*(Pmid-Pf/K);assert math.isclose(Q1,Q2,rel_tol=1e-13);checks.append('network:'+s)
assert math.isclose(rows[1]['g_if_K1e9_for_this_species'],g,rel_tol=1e-13);checks.append('target_g_H2')
for r in rows:assert 0<r['thermal_reverse_flux_fraction']<1 and r['warm_to_cold_aperture_fraction_at_tau0p3']>0;checks.append('cold_flux:'+r['species'])
assert conduit['conductance_margin']>1;checks.append('concrete_conduit_flow_margin')
assert all(x['thin_ring_hoop_Pa']>x['uniform_isotropic_solid_disk_max_Pa_nu0p2'] for x in stress);checks.append('stress_reference')
out={'scope':__doc__,'code_sha256':sha(__file__),'benchmark_sha256':sha(P/'benchmarks/results.json'),'reference':point,'species':rows,'network_requirements':network,'compression_similarity':compression,'pressure_backstream_requirements':pressure,'ideal_drag_channels':channel,'rotor_stress':stress,'bearing_drive_loss_requirements':loss,'thermal_isolation_requirements':thermal,'local_molecular_windage_screens':wind,'rework_energy_sensitivity':maintenance,'concrete_downstream_conduit':conduit,'idle_loss_primary_comparison':idle,'checks_passed':len(checks),'checks':checks,'limits':['No1000K pump, bearing, motor, materialstrength or bladecapture is experimentallyvalidatedbythispacket.','Commercial capture values are a transparent performance target atpreservedmolecularspeedratio; their actualhotperformance requires moving-walltransport.','A useful coldoutlet requires thermaltranspiration/reversefluxaccounting, not merely a volume-flow ratio.','Forelinepartialpressure must be maintained independently of directcatcher pulsepressure.','Motor/bearingidleloss, hotstatorclearance, radiationbaffletransmission and contamination protection are substantive apparatus requirements.','The current aperture/mass calculations do not claim a posteriorprobability distribution; no evidence-backed80% credibleinterval for the completehotpump is available.']}
(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'reference':point,'species':rows,'checks':checks},indent=2))
