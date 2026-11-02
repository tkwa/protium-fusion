"""Read archived station/pulse data; derive bounded sector/support requirements.
No fusion or particle stopping simulation. Requirements and chosen material
properties are explicit, and component masses are subdivisions of old budgets.
"""
import json, math, hashlib
from pathlib import Path
P=Path(__file__).resolve().parents[1]
YEAR=365.25*86400; G=6.67430e-11; SB=5.670374419e-8
files={'thermal':P/'cold-sink-redesign/cases.json', 'pulse':P/'membrane-receiver/pulse-loading.json', 'timing':P/'cold-sink-redesign/direct_recovery.json', 'old':P/'open-loop-redesign/operating-cases.json'}
raw={k:json.loads(v.read_text()) for k,v in files.items()}
c=raw['thermal']['cases']['working_200K_60pct_IR11_bypass1']; old=raw['old']['cases']['high_source_rating_surface']
R=c['radius_m']; A=c['area_m2']; M=c['mass_kg']; sigma=M/A; Mt=old['tidal_equivalent_mass_kg']; v=1e4
selfg=G*M/(2*R*R); tideg=G*Mt/(2*R*R); g=selfg+tideg; pressure=sigma*g
streamM=M*g*R/v**2; streamE=.5*streamM*v*v; period=2*math.pi*R/v
support={'total_areal_mass_kg_m2':sigma,'self_acceleration_m_s2':selfg,'maximum_tidal_acceleration_m_s2':tideg,'conservative_uniform_support_pressure_Pa':pressure,'stream_speed_m_s':v,'stream_mass_kg':streamM,'stream_kinetic_J':streamE,'orbital_period_s':period,'loss_fraction_orbit':.001,'support_W':.001*streamE/period,'tangential_drag_over_normal_force':.001/(4*math.pi),'stream_areal_mass_kg_m2':streamM/A}
# Isotropic homogeneous-shell reference, using receiver material alone as an
# unrealistically perfect continuous load-bearing skin. It diagnoses the need
# for bracing/support; it does not test a designed sandwich shell.
E=1e10; nu=.2; rhoC=2000.; t=5000/rhoC
support['passive_shell_screen']={'graphite_modulus_Pa_assumed':E,'poisson_ratio_assumed':nu,'thickness_m':t,'classical_buckling_pressure_Pa':2*E/(math.sqrt(3*(1-nu*nu)))*(t/R)**2,'compressive_membrane_stress_Pa':pressure*R/(2*t),'ideal_equal_face_sandwich_separation_m':pressure*R*R/(2*E*t),'scope':'Perfect homogeneous shell and ideal sandwich scaling; real imperfections/core shear reduce buckling strength.'}
# Local rail accounting assumes an isotropically distributed collection of
# closed near-great-circle guides. Total length density 2/pitch is a chosen
# locally two-directional reference, not a proved globe tiling.
rails=[]
for pitch in [100.,1000.,10000.]:
 length_density=2/pitch; lam=streamM/A/length_density
 rails.append({'local_rail_pitch_m':pitch,'length_per_area_m_m2':length_density,'combined_counterstream_linear_kg_m':lam,'normal_force_per_rail_length_N_m':pressure/length_density,'counterstream_mass_flow_each_kg_s':lam*v/2,'kinetic_J_per_rail_m':.5*lam*v*v,'linear_drag_N_m':pressure/length_density*support['tangential_drag_over_normal_force']})
support['local_rails']=rails
# Tangent departure after complete loss of the guide's inward deflection.
# The actual gravity on a stream is at most GM/R²+tide; include it here.
a_depart=v*v/R-G*M/R**2-tideg
support['guide_fault']=[{'radial_clearance_m':gap,'time_to_boundary_s':math.sqrt(2*gap/a_depart),'stream_path_before_boundary_m':v*math.sqrt(2*gap/a_depart)} for gap in [.1,1.,10.,100.]]
pulses={}
for n,s in raw['pulse']['stages'].items():
 impulse=s['totals']['momentum_fluence_exact_nodal_Pa_s']; J=s['exact_nodal']['momentum_kg_m_s']/A
 assert abs(impulse-J)<1e-10
 dv=impulse/sigma
 pulses[n]={'momentum_fluence_Pa_s':impulse,'peak_pressure_Pa':s['smooth']['peak_momentum_flux_Pa'],'interval_s':raw['timing']['stages'][n]['pulse_interval_s_at_reference_rate'],'bulk_radial_delta_v_m_s':dv,'ideal_stop_time_if_all_outward_support_removed_s':dv/g,'impulsive_displacement_to_stop_m':dv*dv/(2*g),'scope':'Outgoing capture only; whole installed areal mass responds, impulse applied instantly, outward support removed. Finite arrival generally reduces peak displacement; localized footprints/launch recoil need their own loads.'}
support['outgoing_capture']=pulses
support['mean_capture_pressure_Pa']=sum(x['momentum_fluence_Pa_s']/x['interval_s'] for x in pulses.values())
support['unsupported_sector_fall']=[{'allowed_differential_displacement_m':d,'time_s':math.sqrt(2*d/g)} for d in [.01,.1,1.]]
# Local two-pipe helium loop. A full-load supply plus full-load return for 2L
# deliberately bounds a distributed header; terminal channels/HX are excluded.
# Darcy f is a declared requirement, pressure losses include 6 velocity heads.
def loop(L,q,T,deltaT,u,p=2e6,stress=1e8,f=.015,K=6.,eta=.8):
 cp=2.5*8.314462618/.004002602; Rs=8.314462618/.004002602; rho=p/(Rs*T)
 Q=q*L*L; mdot=Q/(cp*deltaT); flow_area=mdot/(rho*u); D=math.sqrt(4*flow_area/math.pi); ell=2*L
 dp=(f*ell/D+K)*rho*u*u/2
 wall_t=p*D/(2*stress); wall_volume=math.pi*((D/2+wall_t)**2-(D/2)**2)*ell
 return {'sector_side_m':L,'thermal_W':Q,'mean_T_K':T,'gas_deltaT_K':deltaT,'pressure_Pa':p,'gas_speed_m_s':u,'gas_density_kg_m3':rho,'gas_flow_kg_s':mdot,'main_pipe_bore_m':D,'wall_thickness_m':wall_t,'wall_allowable_Pa':stress,'pipe_wall_density_kg_m3':8000.,'header_wall_kg_m2':wall_volume*8000/(L*L),'gas_inventory_kg_m2':flow_area*ell*rho/(L*L),'pressure_loss_Pa':dp,'pump_W':mdot/rho*dp/eta,'pump_fraction_of_thermal':dp/(rho*eta*cp*deltaT),'scope':'Two full-flow 2L headers; pressure wall only, no channels/HX/supports. 100 MPa allowable at hot temperature is a selected material requirement.'}
qhot=c['hot_engine_input_W_m2']; qcold=c['cold_rejection_W_m2']+c['environment_absorbed_W']/A
loops=[]
for L in [100.,1000.,10000.]:
 hot=loop(L,qhot,980.,20.,10.)
 cold=loop(L,qcold,215.,20.,5.)
 loops.append({'sector_side_m':L,'hot':hot,'cold':cold,'combined_header_wall_and_gas_kg_m2':sum(x['header_wall_kg_m2']+x['gas_inventory_kg_m2'] for x in [hot,cold]),'additional_external_pump_W_m2':(hot['pump_W']+cold['pump_W'])/(L*L),'service_allowance_kg_m2':1000.})
# Thermal bridges carry worst support + largest outgoing pressure (conservative
# static peak sign envelope), with a structural factor 4. End temperatures and
# length explicit; these are axial strength floors, not a complete tie design.
pload=pressure+max(x['peak_pressure_Pa'] for x in pulses.values()); ties=[]
for allowable,k,length in [(1e7,1.,1.),(1e7,10.,1.),(1e8,20.,1.),(1e8,100.,1.)]:
 farea=4*pload/allowable
 ties.append({'allowable_Pa':allowable,'conductivity_W_m_K':k,'length_m':length,'solid_cross_section_fraction':farea,'bridge_heat_W_m2':k*farea*800/length,'tie_mass_kg_m2_at_density3000':3000*farea*length,'scope':'Axial strength/thermal lower bound. Buckling, lateral stability, joints, seismic/stream impact loads require extra structure.'})
black=SB*(1000**4-200**4)
def lmtd(a,b):
    assert a>0 and b>0
    return (a-b)/math.log(a/b) if abs(a-b)>1e-10 else a
loop_UA={'receiver_to_hot_loop':qhot/lmtd(30.,10.),
         'hot_loop_to_950K_isothermal_cycle_bound':qhot/lmtd(40.,20.),
         '230K_isothermal_cycle_to_cold_loop_bound':(.4*qhot)/lmtd(25.,5.),
         'cold_loop_to_200K_radiator':qcold/lmtd(25.,5.)}
thermal={'hot_input_W_m2':qhot,'cold_plus_environment_W_m2':qcold,'loops':loops,'ties':ties,'black_hot_cold_flux_W_m2':black,'max_unbaffled_hot_opening_fraction_for_1W_m2':1/black,'mean_through_receiver_gradient_K_at_k100':qhot*t/100,'receiver_diffusion_time_s_at_alpha5e_5':t*t/5e-5,'unbaffled_circular_ports':[{'count':1000,'radius_m':a,'opening_fraction':1000*math.pi*a*a/A,'bypass_W_m2':black*1000*math.pi*a*a/A} for a in [1e4,1e5,1e6]]}
thermal['added_loop_UA_W_K_m2']=loop_UA
thermal['sum_loop_UA_W_K_m2']=sum(loop_UA.values())
thermal['additional_UA_over_old_external_minimum_W_K_m2']=sum(loop_UA.values())-c['external_exchanger_minimum_UA_W_K']/A
thermal['actual_candidate_source_ports']=[{'count':N,'radius_m':100.,'opening_fraction':N*math.pi*100**2/A,'unbaffled_bypass_W_m2':black*N*math.pi*100**2/A} for N in [100000,1000000]]
# Processing is gross handled material, not irreversible loss. Existing factory
# mass/flow coefficient retained, no lower electric cost -> free smaller factory.
shots={k:YEAR/vv['interval_s'] for k,vv in pulses.items()}; large=sum(vv for k,vv in shots.items() if k!='B1'); allshots=sum(shots.values()); factorycoeff=c['factory_mass_per_repaired_kg_s']; renew=[]
for event,count in [('large',large),('all',allshots)]:
 for depth in [1e-7,1e-6,1e-5,1e-4,1e-3]:
  flow=rhoC*depth*count/YEAR
  renew.append({'events':event,'events_per_year':count,'depth_processed_m':depth,'gross_flow_kg_s_m2':flow,'gross_kg_year_m2':flow*YEAR,'factory_mass_kg_m2_at_retained_coefficient':flow*factorycoeff,'electric_W_m2_at_2MJkg':flow*2e6,'fraction_of_reference_net_power':flow*2e6/c['net_export_W_m2'],'scope':'Hypothetical processed depth, not predicted erosion or imported makeup. Actual process energy and factory throughput remain recipe-specific.'})
maint={'shots_per_year':shots,'baseline_gross_handled_kg_year_m2':c['routine_repaired_kg_s']*YEAR/A,'retained_factory_mass_per_kg_s':factorycoeff,'baseline_factory_kg_m2':c['mass_breakdown_kg']['replacement_factory_capacity']/A,'processed_depth_scenarios':renew,'depth_for_10pct_net_electric_at_2MJkg':{k:.1*c['net_export_W_m2']/(2e6*rhoC*n/YEAR) for k,n in [('large',large),('all',allshots)]},'scope':'2 MJ/kg is effective total cost sensitivity, consistent with accepted high manufacturing efficiency. No atomization charge; feed/product free energy must be declared separately.'}
checks={'support_mass_relative':streamM/c['mass_breakdown_kg']['support_streams']-1,'support_power_relative':support['support_W']/c['support_W']-1,'pressure_force_relative':streamM*v*v/R/(M*g)-1,'mean_capture_pressure_fraction_of_static_support':support['mean_capture_pressure_Pa']/pressure}
assert max(abs(checks[k]) for k in ['support_mass_relative','support_power_relative','pressure_force_relative'])<1e-12
out={'scope':'Conditional repeatable sector first checkpoint; no completed global guide topology or recovery design. Scalar choices do not become measured hardware.', 'input_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in files.items()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reference':{'radius_m':R,'mass_kg':M,'area_m2':A,'net_W_m2':c['net_export_W_m2'],'receiver_kg_m2':5000.,'service_kg_m2':1000.,'insulation_kg_m2':100.,'radiator_kg_m2_emitter':20.},'support':support,'thermal':thermal,'manufacturing':maint,'checks':checks}
Path(__file__).with_name('sector-account.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checks':checks,'support':{k:support[k] for k in ['conservative_uniform_support_pressure_Pa','stream_mass_kg','stream_kinetic_J','support_W','passive_shell_screen','guide_fault']},'outgoing':pulses,'thermal_loops':[{k:x[k] for k in ['sector_side_m','combined_header_wall_and_gas_kg_m2','additional_external_pump_W_m2']} for x in loops],'maintenance_10pct_depth':maint['depth_for_10pct_net_electric_at_2MJkg']},indent=2))
