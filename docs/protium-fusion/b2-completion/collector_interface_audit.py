"""Independent reduced-input/source/timing interface audit (no thermal rerun)."""
import json,sys,hashlib,math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
BASE=Path(__file__).resolve().parent.parent;sys.path.insert(0,str(BASE))
from fast_eos import FastEOS
from independent_cycle_audit import one_stage,STAGES,MU,MEV
R=1e9

def representative(profile):
 r=np.r_[0.,np.asarray(profile['outer_radius_m'])];v=np.r_[0.,np.abs(profile['velocity_m_s'])]
 center=[];speed=[]
 for j in range(len(r)-1):
  rc=.75*(r[j+1]**4-r[j]**4)/(r[j+1]**3-r[j]**3)
  center.append(rc);speed.append(v[j]+(v[j+1]-v[j])*(rc-r[j])/(r[j+1]-r[j]))
 return np.array(center),np.array(speed)

def v_quantile(profile,f,energy=False):
 v=np.r_[0.,np.abs(profile['velocity_m_s'])];dm=np.asarray(profile['mass_fraction_of_total']);sl=dm/np.diff(v)
 denom=sum(sl*(v[1:]**3-v[:-1]**3)/6) if energy else sum(dm)
 def fraction(x):
  hi=np.minimum(v[1:],np.maximum(v[:-1],x))
  val=sum(sl*(hi**3-v[:-1]**3)/6) if energy else sum(sl*(hi-v[:-1]))
  return val/denom
 return brentq(lambda x:fraction(x)-f,0,float(v[-1]),xtol=1e-9)

def run():
 home=BASE/'open-loop-redesign';inp=json.loads((home/'collector-inputs.json').read_text());report=json.loads((home/'collector.json').read_text());process=json.loads((home/'process.json').read_text())['selected'];rows={};occ={m:0. for m in ['mass999','energy99','full_resolved']}
 for cat,item in inp['stages'].items():
  d=item['record'];p=item['initial'];published=report['stages'][cat];stage,protons=STAGES[cat]
  mass=d['mass_kg']*np.asarray(p['mass_fraction_of_total']);rc,v=representative(p);edges_r=np.r_[0.,p['outer_radius_m']];vb=np.abs(p['velocity_m_s']);node=.5*np.r_[mass[:-1]+mass[1:],mass[-1]]
  ke=sum(.5*node*vb**2);ecell=.5*mass*v**2;ecell*=ke/sum(ecell)
  baryons=sum(({'c12':12,'c13':13,'n14':14,'n15':15}[cat],protons));charge=({'c12':6,'c13':6,'n14':7,'n15':7}[cat]+protons)/baryons;yi=(1+protons)/baryons
  rho=mass/(4*np.pi/3*np.diff(edges_r**3))/1000
  Ucell=mass*FastEOS().state(rho,np.asarray(p['T_K']),charge,yi).e*1e-4
  times=-(R-rc)/v;edges=np.r_[times[0]-(times[1]-times[0])/2,(times[:-1]+times[1:])/2,times[-1]+(times[-1]-times[-2])/2];width=np.diff(edges)
  peak=max((ecell+Ucell)/width)/.6
  timing={}
  for name,f,en in [('mass999',.001,False),('energy99',.01,True)]:
   incoming=R/v_quantile(p,f,en);outgoing=R/v_quantile(d['profile'],f,en);timing[name]=incoming+outgoing
  timing['full_resolved']=R/min(representative(p)[1])+R/min(representative(d['profile'])[1])
  for name,t in timing.items():occ[name]+=t*process['stage_pulse_data'][cat]['pulses_per_He']
  # The trace records energy/uunit per gram; convert directly, without fitting
  # its final point to reported radiation energy.
  trace=Path(item['source_paths']['raw_trace']);tr=np.load(trace);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);nn=len(tr['species']);z=tr['state'].reshape(nz,nv,-1)
  rawrad=np.asarray(tr['cell_mass_kg'])@z[:,nn+4,:]*float(tr['energy_unit_erg_g'])*1e-4
  cached=np.asarray(item['prompt_cumulative_J'])
  cold=one_stage(d,{cat:1,'p':protons});coldheat=(cold['qcold']-cold['coldnu'])*MEV*d['mass_kg']/(baryons*MU)
  received=d['final_kinetic_J']+d['final_internal_J']+d['radiation_escaped_J']+coldheat
  differences={'initial_K':ke/d['initial_kinetic_energy_J']-1,'initial_U':sum(Ucell)/d['initial_internal_energy_J']-1,
   'source_peak':peak/published['source_rating']['peak_total_prepared_electric_W']-1,
   'source_integral':sum(((ecell+Ucell)/width)*width)/d['initial_energy_J']-1,
   'prompt_raw_integral':float(np.max(abs(cached-rawrad))/max(abs(rawrad[-1]),1)),
   'received_energy':received/published['received_energy_J']-1}
  differences.update({'timing_'+k:t/published['timing_at_R1e9'][k]['occupied_s']-1 for k,t in timing.items()})
  checks={k:bool(abs(v)<1e-8) for k,v in differences.items()}
  rows[cat]={'relative_differences':differences,'checks':checks,'source_peak_W':float(peak),'timing_s':timing,'trace_sha256_matches':hashlib.sha256(trace.read_bytes()).hexdigest()==item['source_sha256']['raw_trace']}
 totals={k:occ[k]/report['occupied_s_per_He_at_R1e9'][k]-1 for k in occ};totals['peak_source']=max(x['source_peak_W'] for x in rows.values())/report['peak_source_electric_W_at_R1e9']-1
 return {'scope':'Independent source-energy binning, raw prompt integral, timing quantiles, aggregate occupancy, and cold-heat interface. No independent thermal convolution or physical driver/collector validation.',
 'input_sha256':{n:hashlib.sha256((home/n).read_bytes()).hexdigest() for n in ['process.json','collector.json','collector-inputs.json']},
 'stages':rows,'aggregate_relative_differences':totals,
 'all_checks_pass':all(all(x['checks'].values()) and x['trace_sha256_matches'] for x in rows.values()) and all(abs(v)<1e-8 for v in totals.values()),
 'checks_passed':int(sum(sum(x['checks'].values())+int(x['trace_sha256_matches']) for x in rows.values())+sum(abs(v)<1e-8 for v in totals.values())),
 'limits':['Velocity-bin quantiles inherit a chosen continuous interpolation and cannot establish clearing of the central zero-speed tail.',
 'Source total prepared energy includes initial EOS heat and kinetic energy, but does not model their physical delivery, switching, focusing or protection.',
 'Cold heat uses the same conservative neutrino convention as the cycle. Heat deposition and stopping-induced transmutation remain separate.']}
if __name__=='__main__':
 out=run();Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['all_checks_pass','checks_passed','aggregate_relative_differences']},indent=2))
