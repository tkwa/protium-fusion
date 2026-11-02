"""Same-cohort finite-dispatch atom closure using frozen routed arrival arrays.
No hot-buffer replay. The exact exponential cavity return is integrated between
saved impulse times; every emitted atom is either held in G or recaptured C.
"""
from pathlib import Path
import argparse,json,math,hashlib
import numpy as np
from scipy.signal import lfilter
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--gas-dir',type=Path,default=P.parent);a=ap.parse_args();gas=a.gas_dir
files=['inputs.json','route-inputs.json','chamber-dt150.0-q128-cap0.1_eta100000.0.npz','chamber-dt150.0-q128-cap0.1_eta100000.0.json']
d=json.loads((gas/files[0]).read_text());route=json.loads((gas/files[1]).read_text())['cases']['cap0.1_eta100000.0'];z=np.load(gas/files[2]);meta=json.loads((gas/files[3]).read_text());t=z['time_s'];J=z['source_H_He_kg'];dt=float(t[1]-t[0]);A=4*math.pi*d['R_m']**2
I=np.zeros_like(J);checks=[]
def ck(name,ok,value=None):
 checks.append({'name':name,'passed':bool(ok),'value':value});assert ok,(name,value)
for j,sp in enumerate(['p','he4']):
 for cat,s in route.items():
  amount=np.array(s['species'][sp]['packet_nuclei_m2'])*A*d['atomic_mass_kg'][sp];arrival=np.array(s['arrival_after_ignition_s'])
  for e in d['events']:
   if e['cat']!=cat:continue
   ix=np.ceil((e['focus_s']+arrival)/dt).astype(int);valid=(ix>=0)&(ix<len(t));np.add.at(I[j],ix[valid],amount[valid])
 ck(f'local branch nonnegative {sp}',np.min(I[j]-J[j])>=-1e-3,float(np.min(I[j]-J[j])))
L=I-J;Q=np.array(meta['cavity_source_mass_rate_atomic_kg_s']);m=np.array([2*d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']]);v=np.sqrt(8*1.380649e-23*1000/(math.pi*m));tau=4*d['R_m']/(3*v)
rows=[]
for g in [[.01,.01],[.01,.01359546360457702],[.1,.1]]:
 stages=[]
 for j,sp in enumerate(['H','He']):
  rate=g[j]/tau[j];e=math.exp(-rate*dt);G0=Q[j]/rate
  G=G0*np.exp(-rate*t)+lfilter([1.],[1.,-e],J[j])
  C=float(np.sum(G[:-1])*(-math.expm1(-rate*dt)))
  si=float(sum(I[j]));sj=float(sum(J[j]));sl=float(sum(L[j]));delta=float(G[-1]-G0);err=(sl+C+delta-si)/si
  ck(f'cohort closure {sp} g{g[j]}',abs(err)<2e-12,err)
  ck(f'local+cavity partition {sp} g{g[j]}',abs((sl+sj)/si-1)<2e-12)
  atoms=d['atomic_mass_kg']['p' if j==0 else 'he4']
  stages.append({'species':sp,'incident_I_kg':si,'direct_local_L_kg':sl,'cavity_injection_J_kg':sj,'initial_cavity_G_kg':G0,'final_cavity_G_kg':float(G[-1]),'delta_cavity_G_kg':delta,'captured_from_cavity_integral_C_kg':C,'physical_total_local_input_integral_L_plus_C_kg':sl+C,'incident_nuclei':si/atoms,'incident_mean_cavity_fraction':sj/si,'combined_closure_relative_to_incident':err,'false_extra_mass_if_I_plus_C_is_called_same_source_recovery_kg':C,'mean_local_input_minus_incident_kg_s':-delta/(t[-1]-t[0])})
 rows.append({'effective_capture_H_He':g,'species':stages})
out={'schema':1,'scope':'Numerical same-cohort check of the frozen ideal mobile routing only. This finite-record cavity replay uses all4096 packet nodes, exact saved chronological source bins and existing steady initial cavity stock; no new pulse or hot-buffer evolution. No total physical retention probability is assigned.','time_start_s':float(t[0]),'time_end_s':float(t[-1]),'duration_s':float(t[-1]-t[0]),'dt_s':dt,'source_rule':'I is total arrival on the exact routed quadrature. J is the archived routed cavity impulse. L=I-J. Capture is the integral of lambda*G between each impulse, evaluated analytically.','equations':['J=f_packet I; L=(1-f_packet)I','dG/dt=J-C; dB/dt=L+C-O','integral(L+C)+deltaG=integral(I)','d(G+B)/dt=I-O; steady O=I','With a held solid cohort S: dS/dt=L-R; dB/dt=R+C-O, so d(G+B+S)/dt=I-O. Regeneration and delayed release replace immediate release of those atoms.'],'cases':rows,'buffer_capacity_interface':'The frozen100%-incident hot-buffer response is a maximum-throughput capacity case. A separately added global-return envelope can be used for extra capacity margin only. It is not a same-cohort physical input history. Captured earlier steady cavity stock can make finite-interval returns differ from current arrivals; deltaG records it.','routing_roundoff_rule':meta['arrival_quadrature_scope'],'checks':checks,'source_sha256':{f:hashlib.sha256((gas/f).read_bytes()).hexdigest() for f in files},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'lineage.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'time_s':float(t[-1]),'reference':rows[0]},indent=2))
