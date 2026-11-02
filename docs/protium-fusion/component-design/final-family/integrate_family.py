"""Compose actual delivered burns, source preparation, material and receiver ledgers.
Hot records remain physical records. A separate plant allocation view exposes the
precursor energy and profile without changing any burn endpoint or isotope yield.
"""
from pathlib import Path
import argparse,copy,json,hashlib,importlib.util,sys,math,types
import numpy as np
CAT={'n14':'A','n15':'B1','c12':'B2','c13':'C13'}
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def mod(name,p):s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def preparation_view(spec,entry,research,receiver):
 from fast_eos import FastEOS
 import pulse_model as pm
 resolve=lambda p:Path(p) if Path(p).is_absolute() else research/p
 record=resolve(spec['record']);initial=resolve(spec['initial']);trace=resolve(spec['trace']);d=read(record);q=read(initial);hot=entry['record'];u=hot['assembly_energy_ledger'];tr=np.load(trace);mass=np.asarray(tr['cell_mass_kg']);m=hot['mass_kg']*np.asarray(q['mass_fraction_of_total']);assert np.max(abs(m/mass-1))<1e-10
 ri=np.r_[0.,q['outer_radius_m']];vi=np.r_[0.,np.abs(q['velocity_m_s'])];assert np.all(np.diff(ri)>0) and np.all(np.diff(vi)>0)
 Y=np.array([entry['stage_data']['initial_Y'].get(k,0.) for k in pm.SPECIES]);assert abs(Y@pm.AA-1)<1e-12
 rho=m/(4*np.pi/3*np.diff(ri**3));Ui=m*FastEOS().state(rho/1000,np.asarray(q['T_K']),np.full(len(m),Y@pm.ZZ),np.full(len(m),sum(Y))).e*1e-4
 K=float(.25*np.sum(m*(vi[:-1]**2+vi[1:]**2)));I=float(sum(Ui));U=K+I;assert math.isclose(U,u['preparation_initial_FastEOS_J'],rel_tol=1e-10);assert math.isclose(U,d['initial_energy_J'],rel_tol=1e-10)
 rc,vc=receiver.reps(ri,vi);R=1e9;launch=-(R-rc)/vc;edges=np.r_[launch[0]-.5*(launch[1]-launch[0]),.5*(launch[1:]+launch[:-1]),launch[-1]+.5*(launch[-1]-launch[-2])];assert np.all(np.diff(edges)>0)
 centroid=.5*m*vc**2;centroid*=K/sum(centroid);local=.25*m*(vi[:-1]**2+vi[1:]**2)
 rating={'peak_kinetic_electric_W':float(max(centroid/np.diff(edges))/.6),'peak_total_prepared_electric_W':float(max((centroid+Ui)/np.diff(edges))/.6),'peak_total_exact_local_halfcell_allocation_W':float(max((local+Ui)/np.diff(edges))/.6),'electric_storage_J':U/.6,'initial_electrons_per_baryon_min_max':[float(Y@pm.ZZ)]*2,'initial_ions_per_baryon_min_max':[float(sum(Y))]*2,'scope':'Actual common FastEOS precursor profile; same60% source-rating convention and centroid/local comparison as existing receiver reducer. This is an auxiliary spherical source boundary, not qualified cold-grain assembly.'}
 ro=np.r_[0.,hot['profile']['outer_radius_m']];vo=np.r_[0.,hot['profile']['velocity_m_s']];timing=receiver.timing(ri,vi,ro,vo,m,R)
 names=d['state_variables'];state=tr['state'].reshape(int(tr['zones']),int(tr['variables_per_zone']),-1);t=tr['time_s'];cum=m@state[:,names.index('radiated_energy/uunit'),:]*float(tr['energy_unit_erg_g'])*1e-4;assert math.isclose(float(cum[-1]),u['preparation_escaped_radiation_J'],rel_tol=1e-10);assert min(np.diff(cum))>-1e-7*max(cum[-1],1.)
 rad=float(cum[-1]);burn=hot['initial_energy_J'];delta=U-burn;closure=U-burn-rad
 return {'source_prepared_energy_J':U,'source_prepared_internal_J':I,'source_prepared_kinetic_J':K,'burn_start_energy_J':burn,'preparation_escaped_radiation_J':rad,'preparation_input_increment_J':delta,'preparation_handoff_residual_J':closure,'relative_handoff_residual':closure/U,'source_rating':rating,'timing_at_radius':timing,'reference_radius_m':R,'preparation_duration_s':float(t[-1]),'radiation_waveform':{'time_relative_burn_start_s':(t-t[-1]).tolist(),'cumulative_escaped_J':cum.tolist(),'outer_radius_m':(state[-1,names.index('r/R0'),:]*d['R0_m']).tolist(),'scope':'Preparation-only saved grey luminosity accumulator; separate from burn prompt. No spectrum or receiver recovery assigned.'},'profile':q,'initial_Y':dict(entry['stage_data']['initial_Y']),'source_paths':{'record':str(record),'initial':str(initial),'trace':str(trace)},'source_sha256':{str(p):sha(p) for p in [record,initial,trace]},'scope':'Physical hot record unchanged; precursor source allocation/timing is an explicit alternative boundary.'}

def integrate(stage_map,research,work,out,family,thermal_points=65536):
 sys.path.insert(0,str(research));sys.path.insert(0,str(research/'component-design'));import material_cycle as mc
 vi=mod('family_variant',research/'mixture-design/variant_integrate.py');rv=mod('family_receiver',research/'mixture-design/receiver.py');pc=mod('family_plant',research/'mixture-design/plant_compare.py');general=mod('family_general_station',research/'mixture-design/station_generalized.py')
 definition=read(stage_map);rv.self_check();selected,archive,checks,sources=vi.extract_explicit(stage_map,research,work,rv,mc,thermal_points=thermal_points)
 hotcycle=pc.cycle_for(selected);precursor=copy.deepcopy(hotcycle);plantselected=copy.deepcopy(selected);preparations={};delta=escaped=0.
 for e in plantselected:
  cat=pc.CAT[e['record']['stage']];spec=definition['stage_map'][cat];e['source_prepared_energy_J']=e['record']['initial_energy_J']
  if 'preparation' not in spec:continue
  p=preparation_view(spec['preparation'],e,research,rv);preparations[cat]=p;factor=hotcycle['stage_pulse_data'][cat]['pulses_per_He']/mc.MEV;delta+=factor*p['preparation_input_increment_J'];escaped+=factor*p['preparation_escaped_radiation_J'];sources.update(p['source_sha256'])
  e['source_prepared_energy_J']=p['source_prepared_energy_J'];e['storage_allocation_override']={'original_burn_initial_energy_J':e['record']['initial_energy_J'],'source_prepared_energy_J':p['source_prepared_energy_J'],'scope':'Only the legacy plant calculator storage read uses this derived view. The unmodified physical burn record remains in hot_selected/stage_inputs.'};e['record']['initial_energy_J']=p['source_prepared_energy_J'];e['receiver']['hot_start_source_rating']=e['receiver']['source_rating'];e['receiver']['hot_start_timing_at_radius']=e['receiver']['timing_at_radius'];e['receiver']['source_rating']=p['source_rating'];e['receiver']['timing_at_radius']=p['timing_at_radius'];e['receiver']['preparation_source_override']=p['source_paths'];checks.append({'name':cat+' precursor source energy/profile/trace reconstruction','passed':abs(p['relative_handoff_residual'])<1e-10})
 precursor['U_hot_burn_start_MeV_per_He']=hotcycle['U_prepared_MeV_per_He'];precursor['U_prepared_MeV_per_He']+=delta;precursor['preparation_escaped_MeV_per_He']=escaped;precursor['preparation_energy_increment_MeV_per_He']=delta;precursor['energy_partition_residual_MeV_per_He']+=escaped-delta;precursor['preaux_net_MeV_per_He']=hotcycle['preaux_net_MeV_per_He']-delta/.6;precursor['preaux_electric_J_per_He']=precursor['preaux_net_MeV_per_He']*mc.MEV;precursor['preparation_policy']='Driver uses precursor input; preparation radiation gets no useful electric-recovery credit. Its cold/external destination is selected separately. Existing energy_capture_fraction applies only to hot-burn heat boundary.'
 H=precursor['H_nonnu_MeV_per_He'];Uhot=hotcycle['U_prepared_MeV_per_He'];Upre=precursor['U_prepared_MeV_per_He'];effective_capture=.99*(H+Uhot)/(H+Upre)
 precursor['hot_stage_energy_capture_fraction']=.99;precursor['energy_capture_fraction']=effective_capture;precursor['energy_capture_fraction_scope']='Comparison policy:99% of the hot-burn boundary, preparation input increment escaping externally. Other preparation destinations need explicit cold/source partition.'
 frozen=pc.cold
 def adapter(b,c,**pars):return general.evaluate(b,c,energy_capture_fraction=effective_capture,driver_efficiency=.6,**pars)
 pc.cold=types.SimpleNamespace(**{k:getattr(frozen,k) for k in dir(frozen) if not k.startswith('__')});pc.cold.evaluate=adapter
 plants={str(S):pc.best_rate(plantselected,precursor,absorber_sigma=S) for S in [5000.,2000.]}
 for S,p in plants.items():checks+=vi.result_checks(plantselected,precursor,p)
 for k in ['nuclear_Q_MeV_per_He','neutrino_MeV_per_He','H_nonnu_MeV_per_He','fresh_C_H_mass_ratio','gross_H_output_atoms_per_He','bulk_kinetic_MeV_per_He','hot_internal_MeV_per_He','escaped_prompt_MeV_per_He','delayed_decay_heat_MeV_per_He']:checks.append({'name':'unchanged hot nuclear/material/output '+k,'passed':precursor[k]==hotcycle[k]})
 checks.append({'name':'precursor total energy partition with separate preparation escape','value':precursor['energy_partition_residual_MeV_per_He'],'passed':abs(precursor['energy_partition_residual_MeV_per_He'])<1e-8})
 for c in checks:c['passed']=bool(c['passed'])
 assert all(x['passed'] for x in checks),[x for x in checks if not x['passed']]
 result={'schema':1,'family':family,'hot_selected':selected,'plant_selected':plantselected,'stage_inputs':archive,'hot_cycle':hotcycle,'precursor_cycle':precursor,'preparation_escaped_MeV_per_He':escaped,'preparation_energy_increment_MeV_per_He':delta,'preparations':preparations,'capacity_plants_preparation_external_no_recovery':plants,'baseline_plant_policy':{'hot_capture':.99,'effective_total_capture':effective_capture,'driver_efficiency':.6,'heat_efficiency':.6,'preparation_radiation_destination':'external escape','scope':'Comparison point only; parent coupled photon ledger selects its own hot capture and cold/external partition. No preparation-energy credit is silently added.'},'checks':checks,'source_sha256':sources,'code_sha256':{str(p):sha(p) for p in [Path(__file__),research/'mixture-design/variant_integrate.py',research/'mixture-design/receiver.py',research/'mixture-design/plant_compare.py',research/'mixture-design/station_generalized.py',research/'component-design/material_cycle.py']},'scope':__doc__}
 # Hot-account interchange for source spectrum and downstream physical receiver tools.
 cases={};labels={}
 for e in selected:
  label=e['name'];e=copy.deepcopy(e);e['label']=label;cases[label]=e;labels[e['record']['stage']]=label
 write(out/'hot-reduced-inputs.json',{'schema':1,'family':family,'cases':cases,'scope':'Unmodified hot-burn physical records; precursor allocation is separate in integration.json.'})
 firstplant=plants['2000.0'];view={'top_detailed':[{'summary':{'cases':labels},'cycle':hotcycle,'plant':firstplant}]};write(out/'hot-spectrum-interface.json',view)
 photon=mod('family_photon',research/'component-design/current-photon128/current_photon128.py');sp=photon.extract(research,{'cases':cases},view,read(research/'component-design/current-photon128/sources.json'));sp['family']=family;sp['preparation_metadata']={'precursor_prepared_MeV_per_He':Upre,'hot_prepared_MeV_per_He':Uhot,'preparation_escaped_MeV_per_He':escaped,'preparation_energy_increment_MeV_per_He':delta,'scope':'This spectrum remains normalized to H+U_hot. Preparation radiation is separate; no spectrum or hot recovery credited.'};write(out/'spectrum.json',sp);result['spectrum']=sp;result['spectrum_path']=str(out/'spectrum.json');write(out/'integration.json',result)
 print(json.dumps({'family':family,'U_hot':Uhot,'U_precursor':Upre,'preparation_escape_MeV_He':escaped,'fresh_C_H':hotcycle['fresh_C_H_mass_ratio'],'checks':len(checks),'plants':{s:{k:p.get(k) for k in ['feasible','net_W_kg','helium_per_s','sequential_occupancy_mass999']} for s,p in plants.items()}},indent=2));return result
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--stage-map',type=Path,required=True);ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--family',required=True);a=ap.parse_args();integrate(a.stage_map.resolve(),a.research_root.resolve(),a.work_dir.resolve(),a.out_dir.resolve(),a.family)
