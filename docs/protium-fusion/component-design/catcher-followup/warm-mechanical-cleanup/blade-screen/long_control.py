import sys,json,pathlib
base=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(base));from simulate import case
rows=[]
for gas in ['He','H2']:
 rows.append(case(dict(gas=gas,u=700.,theta=45.,L_pitch=60.,open_fraction=.65,alpha=0.,T=1000.),1000000,1500000+2*len(rows)))
(base/'records'/'continuous-long-control.json').write_text(json.dumps({'scope':'A single long specular corridor; no interstage thermalization. This tests why multiplying short-channel K factors needs additional physical reservoirs/stators.','cases':rows},indent=2)+'\n')
for r in rows:print(r['parameters']['gas'],r['K0'],r['net_long_cascade_capture'])
