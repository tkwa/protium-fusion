"""Conditional channel-pitch/diffusivity trade; preserves all original PDE records."""
from pathlib import Path
import sys,json,math,hashlib
import numpy as np
P=Path(__file__).resolve().parent;C=P.parent;F=C.parent
sys.path.insert(0,str(C));from calculate import solve
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
src=read(F/'stopping/surface-stopping.json');fine=read(C/'fine-candidates.json');geo=read(F/'geometry/results.json');gas=read(F/'gas/pumping-interfaces.json');snap=read(F/'barrier-review/input-snapshot.json');ports=read(F/'gas/joint-ports.json');inp=read(F/'gas/inputs.json')
checks=[]
def ck(name,value):checks.append({'name':name,'pass':bool(value)})
def reduce(z,u,t,eta):
 stages={}
 for name,x in src['stages'].items():
  stages[name]={}
  for species,y in x['species'].items():
   dep=np.array(y['packet_mean_projected_depth_um']);n=np.array(y['packet_nuclei_m2']);f=1-np.interp(dep,np.r_[0,z,150.],np.r_[0,u,u[-1]]);times=np.interp(dep,np.r_[0,z,150.],np.r_[0,t,t[-1]])
   stages[name][species]={'front_fraction':float(n@f/n.sum()),'expected_exit_s_if_Dnormal_1um2_s':float(n@times/n.sum()),'minimum_Dparallel_m2_s_for_mean_exit_below1s':float(n@times/n.sum())*eta*1e-12}
 out={}
 for species in ['p','he4']:
  flow={name:sum(x['species'][species]['packet_nuclei_m2'])/x['period_s']for name,x in src['stages'].items()};out[species]=sum(flow[name]*stages[name][species]['front_fraction']for name in flow)/sum(flow.values())
 return {'stages':stages,'flow_weighted_front_fraction':out}
ref=next(x for x in fine['rows']if x['cap_um']==.1 and x['anisotropy']==1e4);strong=next(x for x in fine['rows']if x['cap_um']==.1 and x['anisotropy']==1e5)
# A scaled cell multiplies every matrix and RHS entry by the same factor.
# This is an exact finite-cap similarity, independent of source depth spectrum.
similar=[]
for pitch in [5.,10.,20.]:
 similar.append({'pitch_um':pitch,'eta_to_equal_20um_1e4':1e4*(pitch/20)**2,'eta_to_equal_20um_1e5':1e5*(pitch/20)**2,'reference_front_fraction':ref['resolutions'][-1]['flow_weighted_front_exit_fraction'],'stronger_reference_front_fraction':strong['resolutions'][-1]['flow_weighted_front_exit_fraction']})
z0,u0,t0=solve(20.,.1,150.,1e4,160,64)
for pitch in [5.,10.]:
 z,u,t=solve(pitch,.1,150.,1e4*(pitch/20)**2,160,64);ck(f'exact_similarity_pitch{pitch}',np.max(abs(u-u0))<1e-9 and np.max(abs(t-t0))<1e-8)
rows=[]
for pitch,eta in [(5.,1000.),(10.,1000.)]:
 versions=[]
 for nz,ny in [(320,128),(640,256)]:
  z,u,t=solve(pitch,.1,150.,eta,nz,ny);r=reduce(z,u,t,eta);r.update({'depth_resolution_parameter':nz,'half_pitch_cells':ny});versions.append(r);print(pitch,eta,nz,r['flow_weighted_front_fraction'],flush=True)
 diff=max(abs(versions[0]['stages'][s][k]['front_fraction']-versions[1]['stages'][s][k]['front_fraction'])for s in src['stages']for k in ['p','he4']);ck(f'grid_refinement_pitch{pitch}',diff<1e-5)
 rows.append({'pitch_um':pitch,'eta':eta,'cap_um':.1,'depth_um':150.,'versions':versions,'maximum_stage_species_absolute_refinement_difference':diff})
# Flow requirements: same 100nm slit width,150um depth,50um post-filled plenum,
# .5mm-radius vias at1cm pitch through.175m backing. Slit hydraulic approximation.
kb=1.380649e-23;sb=5.670374419e-8;T=1000.;geometry=geo['geometry'];q=gas['mean_throughput_Pa_m3_s_per_wall_m2'];m=[2*inp['atomic_mass_kg']['p'],inp['atomic_mass_kg']['he4']];v=np.sqrt(8*kb*T/(math.pi*np.array(m)));names=['H2','He'];width=1e-7;depth=150e-6;via_r=.5e-3;via_pitch=.01;via_f=math.pi*via_r**2/via_pitch**2;beta_via=via_f*8*via_r/(3*.175);port=ports['representative_interface'];cold=[];gr=[]
cadence=sum(1/x['period_s']for x in src['stages'].values());rho=2000.;f_leak=.001
for pitch in [20.,10.,5.]:
 mouth=width/(pitch*1e-6);beta_slot=mouth*8*width/(3*depth);Cslot=v/4*beta_slot;Cvia=v/4*beta_via;Cchannel=1/(1/Cslot+1/Cvia)
 for Sp in [.001,.004]:
  S=1/(1/Cchannel+1/Sp);beta=4*S/v;phi=f_leak/(1-f_leak)*beta;pressure=np.array(q)/S
  peak={stage:kb*T*(x['H_peak_flux_m2_s']/2/S[0]+x['He_peak_flux_m2_s']/S[1])for stage,x in snap['stages'].items()}
  void=50e-6*.75+depth*mouth;added=(.175+50e-6*.75+depth*mouth)/(1-via_f)-.175
  r={'pitch_um':pitch,'slot_fraction':mouth,'downstream_speed_m_s':Sp,'slit_C_m_s_H2_He':Cslot.tolist(),'via_C_m_s_H2_He':Cvia.tolist(),'channel_series_C_m_s_H2_He':Cchannel.tolist(),'full_series_speed_m_s_H2_He':S.tolist(),'mean_source_partial_pressure_Pa_H2_He':pressure.tolist(),'mean_source_total_pressure_Pa':float(sum(pressure)),'immediate_release_peak_total_pressure_Pa':peak,'equivalent_front_open_fraction_for_added_leak0p001_H2_He':phi.tolist(),'per_slit_cap_area_equivalent_transmission_for_added_leak0p001_H2_He':(phi/mouth).tolist(),'local_plenum_plus_slit_void_m3_m2':void,'local_void_relaxation_s_H2_He':(void/S).tolist(),'added_physical_thickness_m_preserving350kg_m2':added,'unsupported_caps_mass_kg_m2':rho*1e-7*mouth,'full_face_cap_mass_kg_m2':rho*1e-7,'unsupported_caps_replacement_each_event_W_m2_at2MJkg':rho*1e-7*mouth*cadence*2e6,'connected_hole_production_times_lifetime_relative_to20um_Sp0p001':None}
  gr.append(r);ck(f'series_pitch{pitch}_Sp{Sp}',all(S<Cchannel)&all(S<Sp));ck(f'leak_partition_pitch{pitch}_Sp{Sp}',np.max(abs(phi/(phi+beta)-f_leak))<1e-14)
base=gr[0]
for r in gr:r['connected_hole_production_times_lifetime_relative_to20um_Sp0p001']=(np.array(r['per_slit_cap_area_equivalent_transmission_for_added_leak0p001_H2_He'])/np.array(base['per_slit_cap_area_equivalent_transmission_for_added_leak0p001_H2_He'])).tolist()
# One possible direct cold-panel implementation. Warm compression is a separate
# apparatus option and must not inherit these panel loads or free pumping speed.
for Sp in [.001,.004]:
 cp=v/4*port['shared_mouth_area_fraction']*port['net_molecular_transmission'];sactive=np.array([10*math.sqrt(T/300),13*math.sqrt(T/273.15)]);required=1/(1/Sp-1/cp);area=required/sactive;Q5=area[1]*sb*200**4
 cold.append({'Sp_m_s':Sp,'same_downstream_shared_mouth_fraction':port['shared_mouth_area_fraction'],'same_molecular_transmission':port['net_molecular_transmission'],'port_C_H2_He_m_s':cp.tolist(),'active_area_fraction_Hgetter_Hecryo':area.tolist(),'panel_only_mass_kg_m2_receiver_at10kg_m2_active':float(sum(area)*10),'He_200K_to5K_heat_W_m2':Q5,'He_refrigeration_W_m2_at0p5Carnot_200K_sink':Q5*78,'hot_to_cold_bypass_upper_W_m2':port['shared_mouth_area_fraction']*sb*T**4});ck(f'cold_port_viable{Sp}',all(cp>Sp))
warm={'inlet_speed_m_s':.004,'existing_downstream_line_speed_m_s':.001,'ideal_pressure_ratio_at_equal_throughput':4.,'extra_isothermal_work_W_m2_at1000K_reversible':float(sum(q)*math.log(4)),'extra_electric_W_m2_at50percent_free_energy_efficiency':float(sum(q)*math.log(4)*2),'minimum_unit_admission_area_fraction_H2_He':(4*.004/v).tolist(),'immediate_release_peak_gas_work_W_m2_at50percent':{name:2*math.log(4)*kb*T*(x['H_peak_flux_m2_s']/2+x['He_peak_flux_m2_s'])for name,x in snap['stages'].items()},'scope':'Conditional warm booster between reduced-pressure source line and original downstream line. Same throughput; temperature unchanged. Thermodynamic gas-work floor only: moving surfaces, gas admission, seals, drive losses, hot compatibility, pump mass and peak throughput require apparatus. Existing downstream sorber/optical area can remain if its actual inlet state and load remain unchanged. No measured warm pumping speed or efficiency inferred.'}
capstopped={}
for species in ['p','he4']:
 total=sum(sum(x['species'][species]['packet_nuclei_m2'])/x['period_s']for x in src['stages'].values());num=sum(float(np.array(x['species'][species]['packet_nuclei_m2'])[np.array(x['species'][species]['packet_CSDA_depth_um'])<=.1].sum())/x['period_s']for x in src['stages'].values());capstopped[species]=num/total
finite={'incident_fraction_stopping_in100nm_cap_CSDA_H_He':capstopped,'scope':'Normal straight-CSDA component only. Over slit projections the stopped cap component can return directly, while transmitted projectiles need the previously proposed stopping floor and sealed extraction. This screen does not assign reflection, cap trapping or a floor capture efficiency. Whole-surface mobile front fraction is (1-fslot)*Fsolid + fslot*Fslot, before independent reflection; Fslot includes roof stops and eventual floor return.','floor_or_roof_front_return_sensitivities':[{'pitch_um':p,'slot_fraction':width/(p*1e-6),'extra_front_fraction_if_all_cap_stops_return':{k:width/(p*1e-6)*v for k,v in capstopped.items()},'additional_fraction_if_transmitted_floor_stream_returns_fraction':{str(f):width/(p*1e-6)*f for f in [.001,.01,.1]}}for p in [20.,10.,5.]]}
out={'scope':__doc__,'inputs_sha256':{str(p):sha(p)for p in[Path(__file__),C/'calculate.py',C/'fine-candidates.json',F/'stopping/surface-stopping.json',F/'geometry/results.json',F/'gas/pumping-interfaces.json',F/'gas/joint-ports.json',F/'gas/inputs.json',F/'barrier-review/input-snapshot.json']},'exact_similarity':similar,'new_refinements':rows,'geometry_and_pumping':gr,'finite_slit_component':finite,'direct_cold_panel_comparison':cold,'warm_booster_requirement':warm,'limits':['Already-mobile particles, constant diffusion tensor, mean projected depth with no straggling; zero-width mathematical channels and reflecting rear.','Physical slit width remains100nm. Finite cap over slits, direct floor impacts, surface damage, seals and finite sink backpressure are separate. Relative cap-leak and repair conditions refer only to the slit-supported area; the full face cap remains0.0002kg/m².','Gas conductance uses inherited hydraulic-diameter/slender-channel approximation. Side-plenum distribution and doglegs add resistance. No lower chamber pressure is credited until effective apparatus speed is supplied.','The.001 added plenum-leak fraction is a comparison budget, additional to PDE front loss and primary reflection.','First-exit times scale inversely with actual Dnormal. Required anisotropy alone does not establish liberation, rate or hot irradiated material properties.'],'checks':checks}
(P/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'geometry':gr,'cold':cold,'warm':warm,'checks':checks},indent=2));assert all(c['pass']for c in checks)
