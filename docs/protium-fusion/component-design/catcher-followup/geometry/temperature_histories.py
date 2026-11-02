"""Exact-family representative350 receiver T(t), with conserved source heat.
Linear constant-property backing, chronological main-window schedule and finite
bulk storage. Prompt microsecond peaks are separate interval bounds, not dwell.
"""
from pathlib import Path
import json,math,hashlib,csv,sys
import numpy as np
from scipy.signal import fftconvolve
P=Path(__file__).parent;COMP=P.parents[1];OUT=Path('/home/tkwa/src/protium-fusion/output/catcher-followup/geometry');OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(COMP))
read=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
loader=COMP/'source128_thin_thermal.py';prefix=loader.read_text().split('\nrows=[]\nfor k,cp',1)[0].replace("P/'source-family-carbon-median128/integration.json'","P/'final-receiver-thermal-input.json'");ns={'__file__':str(loader),'__name__':'catcher_history_inputs'};exec(compile(prefix,str(loader),'exec'),ns)
stages=ns['stages'];R=ns['R'];AREA=ns['AREA'];dt=.25;CUT=20000.;rho=2000.;cp=1800.;k=50.;sigma=350.;L=sigma/rho;CA=sigma*cp;alpha=k/(rho*cp)
plantpath=COMP.parent/'mixture-design/selected-plant/point.json';plant=read(plantpath)['plant'];evpath=COMP/'final-family-dispatch/sigma350/events.csv';ev=list(csv.DictReader(evpath.open()));ev.sort(key=lambda x:float(x['focus_s']));tend=read(COMP/'final-family-dispatch/sigma350/schedule.json')['dispatch_end_s'];edges=np.arange(round(CUT/dt)+1)*dt;time=edges[1:];kernel=np.diff(ns['front_step'](edges,L,k,cp,'N'))
# Exact impulse response of nonuniform Neumann modes at t>=dt.
n=np.arange(1,600);lam=alpha*(n*np.pi/L)**2
hp=np.zeros_like(time)
for i in range(0,len(time),4096):hp[i:i+4096]=2/CA*np.sum(np.exp(-lam[:,None]*time[i:i+4096]),axis=0)
curves={};energy_errors={}
for name,s in stages.items():
 qK,qU,qD,_=ns['pulse_bins'](s,edges,False);q=qK+qU+qD
 cumulative=s['prompt_fluence']+np.cumsum(q)*dt
 nonuniform=fftconvolve(q,kernel)[:len(q)]+s['prompt_fluence']*hp
 curves[name]={'t':np.r_[0.,time],'delivered':np.r_[s['prompt_fluence'],cumulative],'nonuniform':np.r_[0.,nonuniform],'q_material':np.r_[0.,qK+qU],'q_decay':np.r_[0.,qD]}
 energy_errors[name]=(s['fluence']-cumulative[-1])/s['fluence'];assert -1e-10<=energy_errors[name]<1e-6
cat={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
for name,s in stages.items():s['period']=plant['pulse_interval_s'][cat[name]]
qmean=sum(s['fluence']/s['period'] for s in stages.values())
et=np.array([float(e['focus_s'])+R/299792458. for e in ev]);enames=[e['stage'] for e in ev];E=np.array([stages[n]['fluence'] for n in enames]);prefixE=np.r_[0.,np.cumsum(E)]
def evaluate(ts):
 ts=np.asarray(ts);j=np.searchsorted(et,ts,side='right');bulk=(prefixE[j]-qmean*ts)/CA;front=np.zeros_like(ts);material=np.zeros_like(ts);decay=np.zeros_like(ts)
 # Truncated histories are replaced by their complete energy in uniform storage.
 # Lost nonuniform tails remain bounded separately by final-receiver-audit.
 for t0,name in zip(et,enames):
  lo=np.searchsorted(ts,t0,side='left');hi=np.searchsorted(ts,t0+CUT,side='right')
  if hi<=lo:continue
  ages=ts[lo:hi]-t0;c=curves[name];bulk[lo:hi]+=(np.interp(ages,c['t'],c['delivered'])-stages[name]['fluence'])/CA
  front[lo:hi]+=np.interp(ages,c['t'],c['nonuniform']);material[lo:hi]+=np.interp(ages,c['t'],c['q_material']);decay[lo:hi]+=np.interp(ages,c['t'],c['q_decay'])
 return bulk,front,material,decay
coarse_dt=60.;tc=np.arange(0,tend+coarse_dt,coarse_dt);bulk,nonuniform,qmat,qdec=evaluate(tc)
# Since heating is nonnegative, bulk cannot fall faster than qmean/CA between
# samples. This positive correction makes the1000K minimum conservative.
floor=1000.-float(min(bulk))+qmean*coarse_dt/CA+qmean*L/(2*k)
steady_back_offset=qmean*L/(6*k);T=bulk+floor+nonuniform+steady_back_offset
np.savez_compressed(OUT/'scheduled-coarse.npz',time_s=tc,bulk_K=bulk+floor,front_K=T,front_skin4_K=T+(qmat+qdec)*1e-4/4,q_material_W_m2=qmat,q_decay_W_m2=qdec)
proxy=read(COMP/'final-receiver-prompt-proxy.json');profiles=[]
for name in ['A','B1','B2','C13']:
 indices=[i for i,n in enumerate(enames) if n==name]
 # Highest warm-background event of each stage; representative of its chronology.
 chosen=max(indices,key=lambda i:float(np.interp(et[i]-.25,tc,T)))
 selections=[('warmest_background',chosen)]+([('ordinary_midyear',min(indices,key=lambda i:abs(et[i]-tend/2)))] if name=='B1' else[])
 if name=='C13':selections.extend((f'event_{j+1}',i) for j,i in enumerate(indices) if i!=chosen)
 for label,i in selections:
  relative=np.arange(-600.,4000.+dt,dt);absolute=et[i]+relative;b,f,m,de=evaluate(absolute);front=b+floor+f+steady_back_offset;skin=front+(m+de)*1e-4/4;skin150=front+(m+de)*1.5e-4/4
  filename=f'{name}-{label}.csv';fp=OUT/filename
  with fp.open('w',newline='') as stream:
   writer=csv.writer(stream);writer.writerow(['relative_to_prompt_arrival_s','absolute_s','bulk_K','front_bulk_k50_K','front_with100um_skin_k4_K','material_W_m2','delayed_decay_W_m2','front_with150um_skin_k4_K']);writer.writerows(zip(relative,absolute,b+floor,front,skin,m,de,skin150))
  pre=float(np.interp(-.25,relative,front));case=proxy['stages'][name]['cases']['surface_flux_equivalent_Teff'];duration=stages[name]['prompt_duration_s'];flash={'duration_s':duration,'start_relative_s':0.,'upper_increment_bulk_K':case['bulk_carbon_k50_cp1800_reference_rise_K'],'upper_increment_skin4_K':case['shallow_lowk4_cp1800_reference_rise_K'],'preflash_background_K':pre,'interpretation':'Conservative rectangular interval at the finite-depth Teff proxy peak. This interval is microseconds long; no peak is held through the material pulse. Radiation/thermal transfer remains a proxy.'}
  short=np.r_[np.linspace(0,duration,40),np.geomspace(duration*(1+1e-9),.25,240)]
  qtail=stages[name]['prompt_fluence']/np.sqrt(math.pi*k*rho*cp*np.maximum(short-duration,1e-30));bulkflash=np.minimum(case['bulk_carbon_k50_cp1800_reference_rise_K'],qtail)
  skintail=stages[name]['prompt_fluence']/np.sqrt(math.pi*4*rho*cp*np.maximum(short-duration,1e-30));skinflash=np.minimum(case['shallow_lowk4_cp1800_reference_rise_K'],skintail)
  shortpath=OUT/f'{name}-{label}-prompt-upper.csv'
  with shortpath.open('w',newline='') as stream:
   writer=csv.writer(stream);writer.writerow(['relative_to_prompt_arrival_s','bulk_front_upper_K','skin4_front_upper_K']);writer.writerows(zip(short,pre+bulkflash,pre+skinflash))
  flash['short_time_upper_csv']=str(shortpath);flash['short_time_upper_sha256']=sha(shortpath);flash['short_time_scope']='Until0.25s, upper history is min(finite-depth peak, wholepromptenergy surfaceimpulse at its finalemissiontime). Use this rather than linearly interpolating the coarse prompttail; the low-k surface bound is conservative for the finite100um skin.'
  profiles.append({'stage':name,'selection':label,'event_index':i,'focus_s':float(ev[i]['focus_s']),'prompt_arrival_s':float(et[i]),'csv':str(fp),'sha256':sha(fp),'maximum_bulk_front_K':float(max(front)),'maximum_skin_front_K':float(max(skin)),'maximum_150um_skin_front_K':float(max(skin150)),'duration_above_K_s':{str(th):float(sum(front>th)*dt) for th in [1100,1200,1300,1400,1500,1600]},'prompt_interval_upper':flash})
out={'scope':__doc__,'code_sha256':sha(__file__),'inputs_sha256':{str(p):sha(p) for p in [loader,COMP/'final-receiver-thermal-input.json',plantpath,evpath,COMP/'final-receiver-prompt-proxy.json',COMP/'final-receiver-audit.json']},'model':{'sigma_kg_m2':sigma,'rho_kg_m3':rho,'cp_J_kg_K':cp,'k_W_m_K':k,'mean_full_source_pickup_W_m2':qmean,'selected_rate_He_s':plant['helium_per_s'],'coarse_background_spacing_s':coarse_dt,'local_profile_spacing_s':dt,'initial_bulk_offset_K':floor,'minimum_sampled_bulk_K':float(min(bulk+floor)),'bulk_missed_minimum_bound_K':qmean*coarse_dt/CA,'steady_back_extraction_front_offset_K':steady_back_offset,'history_cutoff_s':CUT,'pulse_energy_not_yet_arrived_at_cutoff_fraction':energy_errors,'maximum_coarse_front_K':float(max(T)),'time_reference':'Absolute dispatch focus times plus common R/c for prompt arrival; material arrival already references that clock.'},'profiles':profiles,'coarse_history':{'npz':str(OUT/'scheduled-coarse.npz'),'sha256':sha(OUT/'scheduled-coarse.npz')},'limits':['Heat-deposition curves conserve stage energy. Delayed heat is assigned its exponential radioactive time law without the larger envelope used for worst-case thermal bounds; full angular/relativistic arrival smearing is not solved.','All source heat is retained for this hot receiver history, a conservative reference rather than the coupled hot/cold capture partition.','The initial bulk offset is chosen to maintain minimum1000K during the finite archived schedule; startup heat and initial nonuniform state are explicit operating assumptions, not generated energy.','Nonuniform contributions older than20000s are omitted from the representative curves. The separate serialized-window thermal audit bounds total old tails/background corrections at about0.1K for350; add that bound when a strict upper curve is needed.','Microsecond prompt intervals are exported separately. The0.25s curves must not be used to infer an unresolved prompt temperature history.','The k4 skin is a quasi-steady resistance on the material/decay pulse; its millisecond transient, contact resistance and damaged-material kinetics are unresolved.','These are physical-model temperature histories for rate integration, not measured temperatures or a credible interval.']}
(P/'temperature-histories.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'model':out['model'],'profiles':profiles},indent=2))
