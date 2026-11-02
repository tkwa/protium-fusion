"""Passive tagged heavy nuclei on archived spherical T/rho/light-species histories.

The tagged abundance does not change the saved temperature, density, proton or
alpha background. A tagged heavy nucleus follows one-heavy-reactant channels;
two-heavy-reactant and heavy fission channels are separately exposure-bounded.
The native-isotope reconstruction tests trace interpolation, not hydro feedback.
"""
import argparse,json,hashlib,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from numba import njit
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--substeps',type=int,default=2);ap.add_argument('--stage',default='all');ap.add_argument('--cache',type=Path);ap.add_argument('--history-root',type=Path);ap.add_argument('--zero-neutron',action='store_true')
a=ap.parse_args();ROOT=a.root.resolve();sys.path.insert(0,str(ROOT));import pulse_model as pm
HEAVY=np.where(pm.AA>4)[0];NH=len(HEAVY);HM={int(v):i for i,v in enumerate(HEAVY)};TAGS=['c12','c13','n14','n15']
COLS=[HM[pm.SPECIES.index(x)] for x in TAGS];CO=pm.COEFF;RI=pm.RATE_INDEX;PREF=pm.PREF;DEXP=pm.DEXP;REACT=pm.REACT;NR=pm.NR
ONE=[];MULTI=[];PRODUCT=[]
for j,r in enumerate(pm.RATES):
 hs=[i for i in REACT[j] if i>=0 and pm.AA[i]>4]
 ps=[pm.NUCS.index(n) for n in r.products if n.A>4]
 if len(hs)==1 and len(ps)<=1:ONE.append((j,HM[hs[0]],HM[ps[0]] if ps else NH))
 elif hs:MULTI.append(j)
ONE=np.asarray(ONE,dtype=np.int64);MULTI=np.asarray(MULTI,dtype=np.int64)
@njit(cache=True)
def generator(T,rho,y):
 t=T/1e9;tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)]);rate=np.zeros(NR)
 for s in range(len(CO)):
  expo=0.
  for k in range(7):expo+=CO[s,k]*tf[k]
  rate[RI[s]]+=np.exp(expo)
 G=np.zeros((NH+4,NH+4));omitted=np.zeros(NH)
 for jj in range(len(ONE)):
  j,src,dst=ONE[jj];coef=rate[j]*PREF[j]*rho**DEXP[j]
  for rr in REACT[j]:
   if rr>=0 and pm_AA[rr]<=4:coef*=max(y[rr],0.)
  G[src,src]-=coef;G[dst,src]+=coef
  for l in range(3):G[NH+1+l,src]+=STOICH[l,j]*coef
 for j in MULTI:
  for tag in REACT[j]:
   if tag<0 or pm_AA[tag]<=4:continue
   coef=rate[j]*PREF[j]*rho**DEXP[j];skipped=False
   for rr in REACT[j]:
    if rr<0:continue
    if rr==tag and not skipped:skipped=True
    else:coef*=max(y[rr],0.)
   omitted[heavy_lookup[tag]]+=coef
 return G,omitted
pm_AA=pm.AA;STOICH=pm.STOICH;heavy_lookup=np.full(pm.NN,-1,dtype=np.int64)
for i,k in enumerate(HEAVY):heavy_lookup[k]=i
mapping=json.loads((ROOT/'open-loop-redesign/corrected-stage-map.json').read_text());loading=json.loads((ROOT/'membrane-receiver/pulse-loading.json').read_text())
DECAY={'n13':'c13','o14':'n14','o15':'n15','f17':'o17','f18':'o18','ne18':'o18','ne19':'f19','na21':'ne21','na22':'ne22','mg22':'ne22','mg23':'na23'}
answer={};checks=[];sourcehash={};start=time.monotonic()
for native,path in mapping.items():
 meta=json.loads((ROOT/path).read_text());tag=meta['stage']
 if a.stage!='all' and tag!=a.stage:continue
 raw=Path(loading['stages'][tag]['sources']['raw_trace']['path'])
 if a.history_root:
  hp=a.history_root/(tag+'-tracer-history.npz');d=np.load(hp);times=d['time_s'];temp=d['T_K'];rho=d['rho_g_cm3'];yy=d['Y'];weights=d['weights'];spec=d['species'].tolist();nz=len(weights)
  sourcehash[str(hp)]=hashlib.sha256(hp.read_bytes()).hexdigest()
 else:
  d=np.load(raw);spec=d['species'].tolist();nz=int(d['zones']);nv=int(d['variables_per_zone']);times=d['time_s'];z=d['state'].T.reshape((-1,nz,nv));weights=d['cell_mass_kg']/np.sum(d['cell_mass_kg'])
  temp=z[:,:,pm.NN]*meta['T0_K'];rr=np.concatenate([np.zeros((len(times),1)),z[:,:,pm.NN+1]*meta['R0_m']],axis=1)
  rho=d['cell_mass_kg'][None,:]/(4*np.pi/3*np.diff(rr**3,axis=1))/1000;yy=z[:,:,:pm.NN]
 assert spec==pm.SPECIES
 if a.cache:
  a.cache.mkdir(parents=True,exist_ok=True);np.savez_compressed(a.cache/(tag+'-tracer-history.npz'),time_s=times,T_K=temp,rho_g_cm3=rho,Y=yy,weights=weights,species=spec)
 total=np.zeros((NH+4,len(TAGS)));worstprob=0.;omitted_int=np.zeros((nz,NH));zoneout=[]
 for cell in range(nz):
  X=np.zeros((NH+4,len(TAGS)));X[COLS,np.arange(len(TAGS))]=1.
  for it in range(len(times)-1):
   dt=(times[it+1]-times[it])/a.substeps
   for ss in range(a.substeps):
    f=(ss+.5)/a.substeps;T=np.exp((1-f)*np.log(temp[it,cell])+f*np.log(temp[it+1,cell]));rh=np.exp((1-f)*np.log(rho[it,cell])+f*np.log(rho[it+1,cell]));y=(1-f)*yy[it,cell]+f*yy[it+1,cell]
    if a.zero_neutron:y[0]=0.
    G,omit=generator(T,rh,y)
    if not np.all(np.isfinite(G)):raise ValueError('nonfinite tagged rates')
    X=expm(G*dt)@X;omitted_int[cell]+=omit*dt
  total+=weights[cell]*X;zoneout.append(X)
  if cell%16==0:print(json.dumps({'stage':tag,'cell':cell,'network_source_sha256':hashlib.sha256(Path(pm.__file__).read_bytes()).hexdigest(),'rate_coefficients_sha256':hashlib.sha256(CO.tobytes()).hexdigest(),'pynucastro_version':pm.p.__version__,'elapsed_s':time.monotonic()-start}),flush=True)
 nativecol=TAGS.index(native);cat0=yy[0,0,pm.SPECIES.index(native)]
 ref=(weights@yy[-1][:,HEAVY])/cat0
 refdict={pm.SPECIES[i]:float(v) for i,v in zip(HEAVY,ref)}
 cold={}
 for k in range(len(TAGS)):
  dst={}
  for i,v in zip(HEAVY,total[:NH,k]):
   name=pm.SPECIES[i];out=DECAY.get(name,name);dst[out]=dst.get(out,0)+float(v)
  reusable=sum(dst.get(x,0.) for x in TAGS)
  light={pm.SPECIES[l]:float(total[NH+1+l,k]) for l in range(3)}
  allout=dst.copy();allout['p']=light['p']+light['n'];allout['he4']=light['he4']
  baryon=float(pm.AA[pm.SPECIES.index(TAGS[k])]-sum(pm.AA[pm.SPECIES.index(n)]*v for n,v in allout.items()))
  cold[TAGS[k]]={'cold_output_nuclei_per_tag':dst,'light_delta_per_tag':light,'cold_all_species_net_output_per_tag':allout,'baryon_residual_per_tag':baryon,'reusable_CNO_probability':reusable,'noncycle_probability':float(1-reusable),'sink_probability':float(total[NH,k])}
 normerr=float(np.max(np.abs(total[:NH+1].sum(axis=0)-1)));nativeerr=float(np.sum(np.abs(total[:NH,nativecol]-ref)))
 checks.append({'name':tag+' tagged number conservation','passed':normerr<1e-8,'error':normerr})
 checks.append({'name':tag+' positivity','passed':float(np.min(total[:NH+1]))>-1e-12,'minimum':float(np.min(total[:NH+1]))})
 answer[tag]={'native':native,'saved_samples':len(times),'substeps':a.substeps,'temperature_range_K':[float(temp.min()),float(temp.max())],
 'native_hot_reconstruction_L1_per_initial_catalyst':nativeerr,'native_original_hot_output':refdict,'native_reconstructed_hot_output':{pm.SPECIES[i]:float(v) for i,v in zip(HEAVY,total[:NH,nativecol])},
 'cold_tracer':cold,'omitted_multiheavy_max_integrated_hazard_by_tag':{name:float(omitted_int[:,COLS[k]].max()) for k,name in enumerate(TAGS)},
 'omitted_multiheavy_massweighted_hazard_by_tag':{name:float(weights@omitted_int[:,COLS[k]]) for k,name in enumerate(TAGS)},
 'zone_mass_fraction':weights.tolist(),'zone_cold_noncycle_probability':{name:[float(1-sum(X[i,k] for i,n in enumerate([pm.SPECIES[x] for x in HEAVY]) if DECAY.get(n,n) in TAGS)) for X in zoneout] for k,name in enumerate(TAGS)}}
 sourcehash[str(raw)]=loading['stages'][tag]['sources']['raw_trace']['sha256'];sourcehash[path]=hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
result={'scope':__doc__,'tagged_species':TAGS,'heavy_species':[pm.SPECIES[i] for i in HEAVY],'zero_neutron_background_control':a.zero_neutron,'one_heavy_reactions':len(ONE),'omitted_multiheavy_reaction_names':[str(pm.RATES[j]) for j in MULTI],
 'stages':answer,'checks':checks,'all_checks_pass':all(x['passed'] for x in checks),'source_sha256':sourcehash,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'network_source_sha256':hashlib.sha256(Path(pm.__file__).read_bytes()).hexdigest(),'rate_coefficients_sha256':hashlib.sha256(CO.tobytes()).hexdigest(),'pynucastro_version':pm.p.__version__,'elapsed_s':time.monotonic()-start}
a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'completed':list(answer),'elapsed_s':result['elapsed_s'],'checks':result['all_checks_pass']},indent=2))
