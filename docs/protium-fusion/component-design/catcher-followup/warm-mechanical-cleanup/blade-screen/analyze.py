"""Checks and finite-cascade accounting for the archived moving-blade runs."""
from pathlib import Path
import json,math,hashlib
import numpy as np
p=Path(__file__).resolve().parent
checks=[]
def check(name,c):
 checks.append({'name':name,'passed':bool(c)})
 if not c:raise AssertionError(name)
def load(name):return json.loads((p/'records'/name).read_text())
files=['grid.json','controls.json','refine.json','selected-million.json','continuous-long-control.json']
allrows=[x for f in files for x in load(f)['cases']]
for i,r in enumerate(allrows):
 for key in ['forward','reverse']:
  v=r[key];check(f'{i}:{key}:resolved mass',v['transmitted']+v['returned']==v['n'] and v['unresolved']==0)
  check(f'{i}:{key}:specular energy',v['max_specular_relative_energy_error']<1e-12)
  check(f'{i}:{key}:work and heat close',v['max_work_heat_energy_closure_relative_error']<1e-12)
controls=[]
for r in load('controls.json')['cases']:
 f=r['forward'];b=r['reverse'];tf=f['T'];tr=b['T'];se=math.sqrt(tf*(1-tf)/f['n']+tr*(1-tr)/b['n'])
 z=(tf-tr)/se if se else 0
 check('reciprocity:'+str(r['parameters']),abs(z)<5)
 controls.append({'parameters':r['parameters'],'reciprocity_z':z,'K0':r['K0']})
for i,x in enumerate(load('path-validation.json')):check('exact skip versus explicit path:'+str(i),x['outcome_match'] and x['reflections_match'] and x['unresolved']==0)
selected=[]
for r in load('selected-million.json')['cases']:
 par=r['parameters'];f=r['forward'];b=r['reverse'];tf=f['T'];tr=b['T'];n=f['n'];K=tf/tr;o=par['open_fraction'];D=o*(tf-tr)
 se=math.sqrt(tf*(1-tf)/n+tr*(1-tr)/n);se_log=math.sqrt((1-tf)/(n*tf)+(1-tr)/(n*tr))
 item={'parameters':par,'Tf':tf,'Tr':tr,'K0':K,'K0_log_normal95':[K*math.exp(-1.95996398454*se_log),K*math.exp(1.95996398454*se_log)],'K0_log_normal80':[K*math.exp(-1.28155156554*se_log),K*math.exp(1.28155156554*se_log)],'net_capture':D,'net_capture_normal95':[D-o*1.95996398454*se,D+o*1.95996398454*se],'net_capture_normal80':[D-o*1.28155156554*se,D+o*1.28155156554*se]}
 if par['u']>0:
  for goal in [1e8,1e9]:
   N=math.ceil(math.log(goal)/math.log(K));c=o*(tf-tr)/(1-K**(-N));item['stages_for'+str(goal)]={'N':N,'K_total':K**N,'zero_backpressure_capture':c}
   # Independent dense linear network: q=a*p_i-b*p_(i+1).
   # p0=1,pN=0; solve N-1 node flux-continuity equations.
   A=np.zeros((N-1,N-1));rhs=np.zeros(N-1);a=o*tf;bb=o*tr
   for j in range(N-1):
    A[j,j]=a+bb
    if j>0:A[j,j-1]=-a
    else:rhs[j]=a
    if j<N-2:A[j,j+1]=-bb
   pi=np.concatenate(([1.],np.linalg.solve(A,rhs),[0.]));q=a*pi[:-1]-bb*pi[1:]
   check('network formula:'+str(par)+str(goal),np.max(np.abs(q-c))<1e-10)
   check('network continuity:'+str(par)+str(goal),np.max(q)-np.min(q)<1e-10)
   # Actual backpressure ratio1e7; K>=1e8 retains>=90% speed.
   item['stages_for'+str(goal)]['capture_at_pressure_ratio1e7']=c*(1-1e7/K**N)
  # Countercheck simply adding series stationary conductances isn't the pump.
  item['same_baffle_required_aperture_fraction_for_g0p01_at_Ktotal1e9']=1/(D*(.99/.01-1/.03))
  item['same_aperture0p244834_effective_g_at_Ktotal1e9']=.99/(1/.03+1/(.244834*D))
 selected.append(item)
summary={'scope':'Statistical intervals describe independent Monte Carlo sampling only. Wall-scattering probabilities, thermalized interstage reservoirs, tip seals, stator throughput and fabrication are unqualified physical assumptions.','cases':len(allrows),'particles_total':sum(r['forward']['n']+r['reverse']['n'] for r in allrows),'selected_million_per_direction':selected,'zero_speed_and_zero_tilt_controls':controls,'maximum_reciprocity_abs_z':max(abs(x['reciprocity_z']) for x in controls),'max_specular_relative_energy_error':max(r[k]['max_specular_relative_energy_error'] for r in allrows for k in ['forward','reverse']),'max_work_heat_closure_relative_error':max(r[k]['max_work_heat_energy_closure_relative_error'] for r in allrows for k in ['forward','reverse']),'single_continuous_long_corridor':load('continuous-long-control.json'),'hashes':{f:hashlib.sha256((p/'records'/f).read_bytes()).hexdigest() for f in files},'network_equations':{'row':'q=f*(Tf*p_i-Tr*p_(i+1)) times vbar/4, for identical1000K stationary Maxwellian reservoirs','cascade':'K=Tf/Tr; Ktot=K**N; S_eff/(A*vbar/4)=f*(Tf-Tr)/(1-K**(-N))','finite_backpressure':'q/A=(vbar/4)*capture*(p_in-p_out/Ktot)','scope':'Reservoir densities/pressures are wellmixed stationary Maxwellians at sameT; physical stators/thermalizing surfaces and any extra resistance must be supplied.'}}
(p/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(p/'validation.json').write_text(json.dumps({'passed':True,'checks':checks},indent=2)+'\n')
print('checks',len(checks),'cases',len(allrows),'particles',summary['particles_total']);print('max reciprocity z',summary['maximum_reciprocity_abs_z'],'max energy',summary['max_specular_relative_energy_error'])
for r in selected:
 print(r['parameters']['gas'],r['parameters']['alpha'],r['parameters']['u'],'Tf Tr',r['Tf'],r['Tr'],'K',r['K0'],'capture95',r['net_capture_normal95'],'N1e8',r.get('stages_for100000000.0',{}).get('N'),'N1e9',r.get('stages_for1000000000.0',{}).get('N'),'fp',r.get('same_baffle_required_aperture_fraction_for_g0p01_at_Ktotal1e9'))
