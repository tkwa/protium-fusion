"""Extract exact final350 schedule and per-cell endpoint H/He inventories."""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
R=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');P=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'component-design'));import material_cycle as mc
F=R/'component-design/final-family/delivered-B1-C13H2/integration.json';D=R/'component-design/final-family-dispatch/sigma350';f=json.loads(F.read_text());d=json.loads((D/'input.json').read_text());ev=list(csv.DictReader((D/'events.csv').open()));names={'A':'n14','B1':'n15','B2':'c12','C13':'c13'};stages={};hashes={str(p):mc.hashfile(p) for p in [F,D/'input.json',D/'events.csv',D/'schedule.json']}
for s in f['plant_selected']:
 cat=names[s['stage_data']['label']];rec=s['record'];trpath=Path(s['receiver']['source_paths']['trace']);tr=np.load(trpath);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);z=tr['state'].reshape(nz,nv,-1)[:,:,-1];nn=len(tr['species']);sp=tr['species'].tolist();m=tr['cell_mass_kg'];Y=np.maximum(z[:,:nn],0.)
 # Normalize each auxiliary composition just as the saved material model does.
 Y/= (Y@np.array([mc.nuc(x).A for x in sp]))[:,None]
 r=np.r_[0,z[:,rec['state_variables'].index('r/R0')]*rec['R0_m']];v=np.r_[0,z[:,rec['state_variables'].index('v/cs0')]*float(tr['velocity_unit_cm_s'])/100]
 assert np.all(np.diff(r)>0) and np.all(np.diff(v)>0)
 masses={q:(m*Y[:,sp.index(q)]/mc.MU*mc.atomic_kg(q)).tolist() for q in ['p','he4']}
 # Neutron beta decay contributes its actual tiny proton number.
 if 'n' in sp:masses['p']=(np.array(masses['p'])+m*Y[:,sp.index('n')]/mc.MU*mc.atomic_kg('p')).tolist()
 sdata=s['stage_data'];expected={q:sdata['cold_output'].get(q,0)*sdata['attempts_per_pulse']*mc.atomic_kg(q) for q in ['p','he4']};errors={q:sum(masses[q])/expected[q]-1 for q in masses};print(cat,errors);assert abs(errors['p'])<1e-8
 for q in masses:
  correction=abs(sum(masses[q])-expected[q]);assert correction<1e-8*sum(expected.values());masses[q]=(np.array(masses[q])*expected[q]/sum(masses[q])).tolist()
 # Source profile actually used by source-rating adapter, including B1 precursor.
 if cat=='n15':
  prep=f['preparations'][cat];print('preparation keys',list(prep))
  ip=R/'component-design/assembly-audit/records/B1-forward-k0p001.initial.json'
 else:ip=Path(s['receiver']['source_paths']['initial'])
 initial=json.loads(ip.read_text());vi=np.r_[0,np.abs(initial['velocity_m_s'])];mi=rec['mass_kg']*np.array(initial['mass_fraction_of_total']);speedrms=float(np.sqrt(.5*np.sum(mi*(vi[:-1]**2+vi[1:]**2))/sum(mi)))
 stages[cat]={'stage':sdata['label'],'pulse_mass_kg':rec['mass_kg'],'cold_H_atomic_kg':expected['p'],'cold_He_atomic_kg':expected['he4'],'HHe_cell_mass_kg':masses,'cell_mass_kg':m.tolist(),'radius_bounds_m':r.tolist(),'velocity_bounds_m_s':v.tolist(),'endpoint_time_s':float(tr['time_s'][-1]),'gas_cell_total_relative_error':errors,'source_rms_speed_m_s':speedrms,'source_max_speed_m_s':float(max(vi)),'source_initial_profile':str(ip),'source_mass_weights':(mi/sum(mi)).tolist(),'source_velocity_bounds_m_s':vi.tolist(),'receiver_main_window':d['timing'][cat],'cold_output_counts':sdata['cold_output'],'attempts_per_pulse':sdata['attempts_per_pulse']}
 hashes[str(trpath)]=mc.hashfile(trpath);hashes[str(ip)]=mc.hashfile(ip)
for e in ev:
 for k in e:
  if k not in ['cat','stage']:e[k]=float(e[k])
out={'family':f['family'],'plant':d['plant'],'cycle':f['precursor_cycle'],'events':ev,'stages':stages,'atomic_mass_kg':{q:mc.atomic_kg(q) for q in ['p','he4']},'R_m':d['plant']['radius_m'],'source_sha256':hashes,'code_sha256':mc.hashfile(__file__),'scope':'Current final350 actual composition and chronological focus times; no old gas source reused. Per-cell final H/He mass and monotone velocities support mass-conserving arrival projections. Local release, trapping, chemical state and pump capture remain explicit parameters.'};(P/'inputs.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:{x:v[x] for x in ['cold_H_atomic_kg','cold_He_atomic_kg','source_rms_speed_m_s','source_max_speed_m_s','gas_cell_total_relative_error']} for k,v in stages.items()},indent=2))
