"""Support consequences of receiver mass changes at the same archived pulse rate.
This is a support parameter study, not a new selected plant performance claim.
"""
from pathlib import Path
import ast,json,math,hashlib,sys
import numpy as np
P=Path(__file__).parent;ROOT=P.parent
sys.path.insert(0,str(ROOT/'mixture-design'))
import plant_compare as pc
src=Path('/home/tkwa/src/protium-fusion/output/mixture-runs/plant-scan-003.json');pk=json.loads(src.read_text())['top_detailed'][0]
selected=[];hashes={str(src):pc.mc.hashfile(src)}
for name,case in pk['plant']['selected_cases'].items():
 rp=ROOT/'mixture-design/receiver-results'/f'{case}.receiver.json';r=json.loads(rp.read_text());path=Path(r['source_paths']['record']);d=json.loads(path.read_text())
 assert r['radius_m']==1e9 and r['completed_requested_duration']
 selected.append({'name':case,'record':d,'receiver':r,'stage_data':pc.mc.load_stage(path)})
 for f in [rp,path]:hashes[str(f)]=pc.mc.hashfile(f)
modepath=P/'support_modes.py';tree=ast.parse(modepath.read_text())
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['energy2','hessian','case']]
for f in [modepath,Path(pc.__file__),Path(pc.cold.__file__),Path(pc.mc.__file__),P/'sector-account.json',P/'support-topology.json']:hashes[str(f)]=pc.mc.hashfile(f)
base=json.loads((P/'sector-account.json').read_text());R=base['reference']['radius_m'];n2=base['support']['maximum_tidal_acceleration_m_s2']/(3*R);G=6.67430e-11
v=1e4;omega=v/R;threshold=8*n2*R**3/G
old_eps=max(json.loads((P/'support-topology.json').read_text())['geometry']['axis_surface_stream_kg_m2'].values())/base['support']['total_areal_mass_kg_m2']
rows=[]
for sigreceiver in [5000.,2000.,1000.]:
 plant=pc.station(selected,pk['cycle'],pk['plant']['helium_per_s'],absorber_sigma=sigreceiver)
 mass=plant['mass_kg'];sigma=mass/plant['area_m2'];a=G*mass/(2*R**3)
 stream=sigma*np.array([(a+4*n2)/2,(a-2*n2)/2,(a-4*n2)/2])/omega**2
 valid=bool(min(stream)>0);modes=[]
 if valid:
  ns={'np':np,'math':math,'R':R,'a':a,'n2':n2,'n':math.sqrt(n2),'sigma':sigma}
  exec(compile(ast.Module(body=nodes,type_ignores=[]),str(modepath)+'::definitions','exec'),ns)
  for days in [30.,100.,300.]:
   r=ns['case'](1e8,True,days)
   modes.append({k:r[k] for k in ['active_attitude_response_days','positive_growth_detected','largest_growth_per_s','H_eigenvalues_in_units_a']})
  assert not any(x['positive_growth_detected'] for x in modes)
 # For1microradian/100day illustrative correction, use differential
 # counterstream angular momentum as a reaction-wheel reservoir.
 tau=100*86400;angle=1e-6;Lbody=(2/3)*mass*R*R*angle/tau
 Lcap=(2/3)*stream*plant['area_m2']*omega*R*R
 delta=np.divide(Lbody,Lcap,out=np.full(3,np.nan),where=Lcap>0)
 row={'receiver_sigma_kg_m2':sigreceiver,'mass_kg_at_same_rate99pct_capture_reference':mass,'total_installed_kg_m2':sigma,'self_coefficient_a_s2':a,'a_over_solar_n2':a/n2,
      'positive_three_family_stream_solution':valid,'axis_stream_kg_m2_diagnostic':stream.tolist(),'heaviest_local_stream_to_body_ratio':max(stream)/sigma,'ratio_to_old_local_control_epsilon':max(stream)/sigma/old_eps,
      'static_support_pressure_upper_Pa':sigma*(a+3*n2)*R,'global_affine_control_cases':modes,
      'illustrative_attitude_correction':{'angle_rad':angle,'response_s':tau,'body_angular_momentum_change_kg_m2_s':Lbody,'per_axis_counterstream_maximum_net_L':Lcap.tolist(),'required_differential_speed_fraction':[(float(x) if math.isfinite(x) else None) for x in delta]},
      'minimum_uniform_prestress_N_m_to_restore_nonnegative_family_densities':max(0.,(4*n2-a)*sigma*R*R/2),
      'scope':'Same archived source/storage/service budgets and rate.99% capture is retained only to obtain the comparison mass; photon qualification is separate. Negative diagnostic stream density invalidates this topology rather than representing negative material.'}
 rows.append(row)
out={'scope':__doc__,'input_sha256':hashes,'code_sha256':pc.mc.hashfile(__file__),'reference_cases':pk['plant']['selected_cases'],'reference_He_s':pk['plant']['helium_per_s'],'radius_m':R,'solar_n2_s2':n2,'mass_threshold_for_positive_three_family_streams_kg':threshold,'extensional_modulus_N_m':1e8,'rows':rows,
 'limitations':['The global9DOF model and local Fourier feedback screen are separate reductions; a full discrete network remains unclosed.',
 'The attitude correction is an amplitude normalization, not a claimed required pointing tolerance. Startup angular momentum must be supplied during construction; internal wheels only redistribute it.',
 'The prestress value is the algebraic minimum for the reference force construction; a positive margin and its membrane/guide rebalancing are required before using a below-threshold plant.',
 'Fixed guide package mass is retained. Smaller support-stream mass is not credited as a free reduction of cryostats, conductor, controls, source or factory hardware.']}
Path(__file__).with_name('receiver-support-family.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'threshold':threshold,'rows':[{k:r[k] for k in ['receiver_sigma_kg_m2','mass_kg_at_same_rate99pct_capture_reference','a_over_solar_n2','positive_three_family_stream_solution','ratio_to_old_local_control_epsilon','minimum_uniform_prestress_N_m_to_restore_nonnegative_family_densities']} for r in rows]},indent=2))
