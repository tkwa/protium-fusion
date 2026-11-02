"""Bounded whole-dispatch audit of the frozen representative-k comparison.
All75largepulse windows are resolved at0.25s. Remainingtimes have a positive
bulk-plus-source-tail bound using exacteventspacing and20sbackgroundsampling.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
P=Path(__file__).resolve().parent;main=P/'calculate.py';src=main.read_text().split('rows=[evaluate_k',1)[0];ss={'__file__':str(main),'__name__':'all_event_import'};exec(compile(src,str(main),'exec'),ss);hist=ss['ns'];d=json.loads((P/'results.json').read_text());ks=[x['required_effective_k_W_m_K'] for x in d['chosen_screens']]+[1000.];events=hist['ev'];et=hist['et'];names=hist['enames'];big=[i for i,n in enumerate(names) if n!='B1'];dt=.25;step=20.;rel=np.arange(.25,4000.+dt,dt);tf=np.concatenate([et[i]+rel for i in big]);owner=np.repeat(big,len(rel));assert np.all(np.diff(tf)>0);tend=float(hist['tend']);tc=np.arange(0,tend+step,step);j=np.searchsorted(et[big],tc,side='right')-1;outside=(j<0)|(tc-et[np.array(big)[np.maximum(j,0)]]>4000);rows=[];checks=[]
for k in ks:
 rep=ss['evaluate_k'](k);floor=rep['new_initial_bulk_offset_K'];off=ss['qmean']*ss['L']/(6*k);R=ss['Rtot'];b,n,m,de=hist['evaluate'](tf);T=b+floor+n+off+(m+de)*R;im=int(np.argmax(T));ei=int(owner[im]);maxf=float(T[im]);bc,nc,mc,dc=hist['evaluate'](tc);Tc=bc+floor+nc+off+(mc+dc)*R
 # Positive source-supremum bound outsideevery large pulse's first4000s.
 times=hist['curves']['A']['t'];B1times=np.array([et[i] for i,nm in enumerate(names) if nm=='B1']);gap=float(min(np.diff(B1times)));bound_nonuniform=0.;bound_heat=0.;components={}
 for name,c in hist['curves'].items():
  S=c['nonuniform']+R*(c['q_material']+c['q_decay']);Q=c['q_material']+c['q_decay'];revS=np.maximum.accumulate(S[::-1])[::-1];revQ=np.maximum.accumulate(Q[::-1])[::-1]
  if name=='B1':
   ids=[int(np.searchsorted(times,z)) for z in np.arange(0,20000,gap)];v=float(sum(revS[ii] for ii in ids));q=float(sum(revQ[ii] for ii in ids))
  else:
   tt=np.array([et[i] for i,nm in enumerate(names) if nm==name]);assert min(np.diff(tt))>20000;ix=int(np.searchsorted(times,4000.));v=float(revS[ix]);q=float(revQ[ix])
  bound_nonuniform+=v;bound_heat+=q;components[name]={'front_nonuniform_plus_surface_bound_K':v,'heat_flux_bound_W_m2':q}
 bulk_sample=float(max((bc+floor)[outside]));bulk_sampling_guard=step*bound_heat/ss['CA'];tailguard=.1;outside_bound=bulk_sample+bulk_sampling_guard+off+bound_nonuniform+tailguard
 checks.append({'name':f'outside_sample_under_bound_k{k}','passed':bool(max(Tc[outside])<=outside_bound)});assert max(Tc[outside])<=outside_bound
 chosen_ids={x[0]['event_index'] for x in ss['profiles']};detail=[]
 for pr,tt,ab,bb,mm,dd,old in ss['profiles']:
  b0,n0,m0,d0=hist['evaluate'](ab);V=b0+floor+n0+off+(m0+d0)*R;mask=tt>=.25;iall=int(np.argmax(V));imask=int(np.flatnonzero(mask)[np.argmax(V[mask])]);detail.append({'stage':pr['stage'],'selection':pr['selection'],'unmasked_argmax_index':iall,'masked_argmax_index':imask,'detail_index_matches':iall==imask,'unmasked_time_s':float(tt[iall]),'masked_time_s':float(tt[imask]),'unmasked_max_K':float(V[iall]),'masked_max_K':float(V[imask])})
 rows.append({'k_W_m_K':k,'representative_max_K':rep['maximum_material_front_over_profiles_K'],'all_large_windows_max_K':maxf,'difference_from_representative_K':maxf-rep['maximum_material_front_over_profiles_K'],'controlling_event_index':ei,'controlling_stage':names[ei],'controlling_focus_s':float(et[ei]-hist['R']/299792458.),'controlling_relative_prompt_time_s':float(tf[im]-et[ei]),'controller_was_in_old7':ei in chosen_ids,'coarse_global_max_K':float(max(Tc)),'maximum_sampled_outside_large_windows_K':float(max(Tc[outside])),'outside_bound_K':outside_bound,'outside_bound_components':components,'outside_bulk_sampling_guard_K':bulk_sampling_guard,'retained_old_tail_guard_K':tailguard,'outside_sampled_bulk_max_K':bulk_sample,'minimum_B1_focus_gap_s':gap,'profile_detail_index_check':detail,'scope':'Finequartersecondalllargeevents pluspositiveoutsidebound. The0.1Khistoricaltailguard appliesoutside; largewindow valuesretain theiroriginal separatelyqualifiedtailtruncation.'})
table_detail=[]
for row in d['conductivity_cases']:
 for pp in row['profiles']:
  ar=np.genfromtxt(pp['history_file'],delimiter=',',names=True);temp=ar['skin_plus_posts_front_K'];tt=ar['relative_s'];ia=int(np.argmax(temp));ii=int(np.flatnonzero(tt>=.25)[np.argmax(temp[tt>=.25])]);table_detail.append({'k_W_m_K':row['bulk_effective_k_W_m_K'],'isothermal':row['isothermal_limit'],'stage':pp['stage'],'selection':pp['selection'],'unmasked_argmax_time_s':float(tt[ia]),'masked_argmax_time_s':float(tt[ii]),'detail_matches':ia==ii,'reported_phase_max_K':pp['max_front_material_phase_K'],'independent_masked_max_K':float(temp[ii])})
out={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'main_code_sha256':hashlib.sha256(main.read_bytes()).hexdigest(),'main_results_sha256':hashlib.sha256((P/'results.json').read_bytes()).hexdigest(),'event_count':len(events),'large_events_finely_resolved':len(big),'fine_samples':len(tf),'fine_dt_s':dt,'coarse_dt_s':step,'cases':rows,'frozen_table_detail_index_all_k':table_detail,'checks':checks,'limits':['Thisauditsresolvedmaterial/decaytemperature,notmicrosecondphotondeposition.','Source-tailboundusesminimumactualB1spacingandpositivekernel; it deliberatelyoverallocatesnonuniformsource contributions.','New k values do not replace the frozenrepresentativepacket; anymissedcontrollerisreportedexplicitly.']};(P/'all-event-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([{k:r[k]for k in ['k_W_m_K','representative_max_K','all_large_windows_max_K','difference_from_representative_K','controlling_event_index','controlling_stage','controller_was_in_old7','outside_bound_K']}for r in rows],indent=2))
