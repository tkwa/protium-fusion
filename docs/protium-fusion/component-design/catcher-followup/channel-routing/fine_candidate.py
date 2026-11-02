"""Refine selected ideal channel geometries and export packet-wise mobile routing."""
from pathlib import Path
import json,hashlib
import numpy as np
from calculate import solve,sha
P=Path(__file__).resolve().parent
src=P.parent/'stopping/surface-stopping.json';d=json.loads(src.read_text())
rows=[];packets={}
for cap,eta in [(.1,1e4),(.1,1e5),(.01,1e4),(.01,1e5)]:
 versions=[]
 for nz,ny in [(320,128),(640,256)]:
  z,u,t=solve(20.,cap,150.,eta,nz,ny); stages={};global_f={}
  for name,x in d['stages'].items():
   st={}
   for sp,y in x['species'].items():
    dep=np.asarray(y['packet_mean_projected_depth_um']);n=np.asarray(y['packet_nuclei_m2']);prob=1-np.interp(dep,np.r_[0,z,150.],np.r_[0,u,u[-1]])
    tt=np.interp(dep,np.r_[0,z,150.],np.r_[0,t,t[-1]])
    st[sp]={'mobile_front_exit_fraction':float(n@prob/n.sum()),'expected_first_exit_s_if_Dnormal_1um2_s':float(n@tt/n.sum()),'minimum_Dparallel_m2_s_for_mean_first_exit_below1s':float(n@tt/n.sum())*eta*1e-12}
    if nz==640:packets[f'cap{cap}_eta{eta}_{name}_{sp}']={'packet_mobile_front_exit_fraction':prob.tolist(),'source_packet_nuclei_m2':n.tolist(),'period_s':x['period_s']}
   stages[name]=st
  for sp in ['p','he4']:
   flow={s:sum(x['species'][sp]['packet_nuclei_m2'])/x['period_s'] for s,x in d['stages'].items()}
   global_f[sp]=sum(stages[s][sp]['mobile_front_exit_fraction']*f for s,f in flow.items())/sum(flow.values())
  versions.append({'depth_resolution_parameter':nz,'half_pitch_cells':ny,'stages':stages,'flow_weighted_front_exit_fraction':global_f})
 diff=max(abs(versions[0]['stages'][s][sp]['mobile_front_exit_fraction']-versions[1]['stages'][s][sp]['mobile_front_exit_fraction']) for s in d['stages'] for sp in ['p','he4'])
 rows.append({'cap_um':cap,'pitch_um':20.,'active_depth_um':150.,'anisotropy':eta,'resolutions':versions,'maximum_absolute_refinement_difference':diff})
 print(cap,eta,versions[-1]['flow_weighted_front_exit_fraction'],diff,flush=True)
result={'scope':'Ideal already-mobile routing only: finite cap + perfect zero-width side collection, uniform lateral source, reflecting rear, per-energy mean-depth approximation. Total physical capture includes independent prompt reflection, traps, straggling, pores, backpressure and collector return.',
 'source_sha256':sha(src),'solver_sha256':sha(P/'calculate.py'),'code_sha256':sha(__file__),'rows':rows,'packets':packets,
 'checks':{'grid_refinement_absolute_below1e-5':max(r['maximum_absolute_refinement_difference'] for r in rows)<1e-5}}
(P/'fine-candidates.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(result['checks']);assert all(result['checks'].values())
