import sys,json,pathlib,math
base=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(base));from simulate import trace,case
rows=[]
for theta in [20,45,60]:
 for u in [0.,700.]:
  for direction in [-1,1]:
   a=trace(n=3000,gas='He',theta=theta,u=u,L_pitch=3.,direction=direction,seed=3456,accelerate=True)
   b=trace(n=3000,gas='He',theta=theta,u=u,L_pitch=3.,direction=direction,seed=3456,accelerate=False,max_events=300000)
   row={'theta':theta,'u':u,'direction':direction,'n':3000,'outcome_match':a['outcome_sha256']==b['outcome_sha256'],'reflections_match':a['reflection_count_sha256']==b['reflection_count_sha256'],'unresolved':a['unresolved']+b['unresolved'],'work_mean_difference_J':a['mean_wall_mechanical_work_per_incident_J']-b['mean_wall_mechanical_work_per_incident_J'],'maximum_reflections':max(a['maximum_reflections'],b['maximum_reflections'])}
   assert row['outcome_match'] and row['reflections_match'] and row['unresolved']==0
   rows.append(row)
p=base/'records';(p/'path-validation.json').write_text(json.dumps(rows,indent=2)+'\n');print('path comparisons passed',len(rows),flush=True)
selected=[]
for gas in ['He','H2']:
 for alpha in [0.,.5,1.]:
  selected.append(case(dict(gas=gas,u=700.,theta=45.,L_pitch=1.,open_fraction=.65,alpha=alpha,T=1000.),1000000,1300000+len(selected)*2))
# High-stat independent rerun of the most deviant initial reciprocity control.
selected.append(case(dict(gas='H2',u=0.,theta=45.,L_pitch=2.,open_fraction=.65,alpha=1.,T=1000.),1000000,1400000))
(p/'selected-million.json').write_text(json.dumps({'n_per_direction':1000000,'cases':selected},indent=2)+'\n');print('million-particle comparisons',len(selected))
for x in selected: print(x['parameters'],x['K0'],x['net_long_cascade_capture'])
