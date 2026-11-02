"""Steady isothermal carbon flights between fully hood-covered spherical walls.
Local specular/diffuse billiards + irreversible sticking; sphericalflight angle
is preserved on return, with randomcellphase/azimuth. No isotropicreturnshortcut.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
P=Path(__file__).resolve().parent;d=json.loads((P/'results.json').read_text());f=1-(.01/.01001)**2;direct=d['direct_flight_radial_column_bounds'][0]['direct_radial_column_weight_per_floor_emission'];seed=994402

def cosdir(rng,n,axis,sgn):
 a=np.sqrt(rng.random(n));ph=rng.uniform(0,2*np.pi,n);b=np.sqrt(1-a*a);v=np.empty((n,3));v[:,axis]=a*sgn;v[:,(axis+1)%3]=b*np.cos(ph);v[:,(axis+2)%3]=b*np.sin(ph);return v

def run(N,s,alpha,source,seed):
 rng=np.random.default_rng(seed);idx=np.arange(N);pos=np.column_stack([rng.random(N),rng.random(N),np.full(N,0. if source=='floor' else 10.-1e-10)]);v=cosdir(rng,N,2,1.);w=np.zeros(N);firstw=np.zeros(N);dw=np.zeros(N);flights=np.zeros(N,dtype=int);coll=np.zeros(N,dtype=int);it=0;captured=0;backcaptures=0
 while len(idx):
  bound=np.where(v>0,np.array([1.,1.,10.]),np.zeros(3));ds=np.divide(bound-pos,v,out=np.full_like(v,np.inf),where=np.abs(v)>1e-15);ds=np.where(ds>1e-12,ds,np.inf);ax=np.argmin(ds,axis=1);dl=ds[np.arange(len(idx)),ax];pos+=v*dl[:,None];out=(ax==2)&(v[:,2]>0)
  oo=np.flatnonzero(out)
  if len(oo):
   ii=idx[oo];mu=v[oo,2];aw=np.arcsin(mu)/(2*np.sqrt(np.maximum(1-mu*mu,1e-30)));w[ii]+=aw;first=flights[ii]==0;firstw[ii[first]]=aw[first];directmask=first&(coll[ii]==0);dw[ii[directmask]]=aw[directmask];flights[ii]+=1
   # Exactspheregeometry returnsatthesamecosine; cellphase andazimuth mix.
   ph=rng.uniform(0,2*np.pi,len(oo));tr=np.sqrt(1-mu*mu);v[oo]=np.column_stack([tr*np.cos(ph),tr*np.sin(ph),-mu]);pos[oo,:2]=rng.random((len(oo),2));pos[oo,2]=10.-1e-10
   hit=rng.random(len(oo))<f;hitidx=oo[hit];capt=rng.random(len(hitidx))<s;dead=hitidx[capt];captured+=len(dead);backcaptures+=len(dead)
   # Unstucktip events reemit into cavity andcreate anotherrealflight.
   rest=hitidx[~capt];rr=rng.random(len(rest))<alpha;v[rest,2]*=-1
   if sum(rr):v[rest[rr]]=cosdir(rng,sum(rr),2,1.)
  else:dead=np.array([],dtype=int)
  # Local solidwall/floor encounters.
  wall=np.flatnonzero(~out);coll[idx[wall]]+=1;ab=rng.random(len(wall))<s;kill=wall[ab];captured+=len(kill);live=wall[~ab];diff=rng.random(len(live))<alpha
  for axis in range(3):
   aa=live[(ax[live]==axis)&diff]
   if len(aa):v[aa]=cosdir(rng,len(aa),axis,-np.sign(v[aa,axis]))
  aa=live[~diff];v[aa,ax[aa]]*=-1;pos[live]+=1e-10*v[live]
  mask=np.ones(len(idx),bool);mask[np.r_[dead,kill]]=False;idx=idx[mask];pos=pos[mask];v=v[mask];it+=1
  if it>10000:raise RuntimeError('residual')
 assert captured==N
 correction=direct if source=='floor' else 1.;control=w-dw+correction if source=='floor' else w-firstw+1
 return {'source':source,'rays':N,'sticking_probability':s,'diffuse_fraction':alpha,'raw_total_radial_column_factor':float(np.mean(w)),'raw_firstflight_radial_column_factor':float(np.mean(firstw)),'direct_firstflight_MC_factor':float(np.mean(dw)),'analytic_direct_firstflight_factor':correction,'control_variate_total_radial_column_factor':float(np.mean(control)),'control_variate_sampling_SE':float(np.std(control,ddof=1)/math.sqrt(N)),'mean_global_flights':float(np.mean(flights)),'probability_at_least_one_global_flight':float(np.mean(flights>0)),'mean_extra_global_flights':float(np.mean(np.maximum(0,flights-1))),'maximum_global_flights':int(max(flights)),'all_atoms_captured':captured,'incoming_tip_captures':backcaptures,'iterations':it}
rows=[]
for s,a in [(.1,1.),(.25,1.),(.5,1.),(1.,1.),(.25,.5),(.25,0.)]:
 for source in ['floor','tip']:rows.append(run(1_000_000,s,a,source,seed+int(s*1000)+int(a*10)+(9999 if source=='tip' else 0)))
combined=[]
for st,al in [(.1,1.),(.25,1.),(.5,1.),(1.,1.),(.25,.5),(.25,0.)]:
 floor=next(r for r in rows if r['sticking_probability']==st and r['diffuse_fraction']==al and r['source']=='floor');tip=next(r for r in rows if r['sticking_probability']==st and r['diffuse_fraction']==al and r['source']=='tip')
 for ratio in [1.,3.,10.]:
  c=(1-f)*floor['control_variate_total_radial_column_factor']+f*ratio*tip['control_variate_total_radial_column_factor'];se=math.hypot((1-f)*floor['control_variate_sampling_SE'],f*ratio*tip['control_variate_sampling_SE']);combined.append({'s':st,'alpha_diffuse':al,'tip_to_floor_yield_ratio':ratio,'wholeface_steady_radial_column_relative_to_bare_unitcapture_cosine_case':c,'sampling_standard_error':se,'radial_column_suppression':1/c})
out={'scope':__doc__,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((P/'results.json').read_bytes()).hexdigest(),'wall_tip_fraction':f,'rows':rows,'combined':combined,'scope_limits':['Thecarbonwalltemperatureandemissionspeedlawareidenticalthroughout; changedtemperature/energyrequiresE1/vweighting.','Thisisasteadycolumnfactor,not theactualpulsedcentralcollisionhistory.','Allreceiverareaishoodcoveredwithoneorientation-independentstickinglaw; exposedpumpports,bafflesandmaterialgaspathsmustbeadded.','Analyticdirectfirstflightcontrolvariatereducesrareaxialrayvariance; samplingSEdoesnotrepresentphysicaluncertainty.','Noinstantglobalisotropizationismade: successive sphericalflightanglespersistuntilalocaldiffusereflection.','Everycarbonatomiscapturedinternallyeventually; thisisretainedinventory,notfreshcarbonloss.']};(P/'cavity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(combined,indent=2))
