"""Convolve exact current-family He packets with conditional diffusion/geometry.
Requires stopping/surface-stopping.json. Mean projected ranges are proxies,
not sampled straggling depths. No retained fraction or anisotropy is selected.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import brentq
import helium_release as h
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 sf=P.parent/'stopping/surface-spectra.json';df=P.parent/'stopping/surface-stopping.json'
 s=json.loads(sf.read_text());d=json.loads(df.read_text());assert d['source_surface_spectra_sha256']==sha(sf)
 stages={};checks=[]
 for name,row in s['stages'].items():
  w=np.array(row['species']['he4']['packet_nuclei_m2']);w=w/w.sum()
  dr=d['stages'][name]['species']['he4'];x=np.array(dr['packet_mean_projected_depth_um'])*1e-6
  assert len(w)==len(x);assert np.allclose(w,np.array(dr['packet_nuclei_m2'])/sum(dr['packet_nuclei_m2']))
  arrival=np.array(row['packet_arrival_after_ignition_s']);period=row['period_s'];D=float(h.diffusivity(1000))
  out={'projected_depth_um_quantiles':{str(q):h.weighted_quantile(x*1e6,w,q) for q in [.001,.01,.1,.5,.9,.99,.999]},'constant1000K_mobile_cases':[],'optimistic_strip_cases':[]}
  for mult in [1e-6,1e-3,1.]:
   def F(t):return float(w@h.released(x,D*mult*np.maximum(t-arrival,0)))
   def ageF(t):return float(w@h.released(x,D*mult*t))
   times={str(q):float(brentq(lambda t:F(t)-q,0,1e13)) for q in [.5,.9,.99]}
   age_times={str(q):float(brentq(lambda t:ageF(t)-q,0,1e13)) for q in [.5,.9,.99]}
   out['constant1000K_mobile_cases'].append(dict(diffusivity_multiplier=mult,
     mobile_fraction_released_at_next_same_stage_ignition=F(period),
     incident_fraction_not_yet_arrived_at_next_same_stage_ignition=float(w[arrival>period].sum()),
     release_time_after_ignition_s=times,
     single_age_ensemble_post_implantation_release_time_s=age_times,
     release_curve=[{'time_after_ignition_s':float(t),'fraction_of_incident_mobile_population_released':F(t)} for t in [100,200,400,1000,period,86400]]))
  for pitch in [20e-6,100e-6]:
   def front(eta):return float(w@h.strip_front_probability(x,pitch,eta))
   req={str(q):float(10**brentq(lambda lg:front(10**lg)-q,-6,20)) for q in [.05,.01,.001]}
   # All-front cap region has no side sinks; gambler's-ruin lower bound on loss.
   caps={str(cap*1e6):float(w@np.maximum(1-x/cap,0)) for cap in [10e-9,100e-9,1e-6]}
   out['optimistic_strip_cases'].append(dict(pitch_m=pitch,
     front_fraction_by_anisotropy={str(e):front(e) for e in [1,10,100,1000,1e4,1e6]},
     required_anisotropy_for_front_fraction=req,
     isotropic_normal_diffusion_cap_front_loss_lower_bound_by_cap_nm={str(float(k)*1000):v for k,v in caps.items()}))
  out['eventual_rear_fraction_isotropic_175mm_slab']=float(w@np.minimum(x/.175,1))
  out['period_s']=period;stages[name]=out
  checks.append({'name':name+' packet counts','passed':bool(abs(w.sum()-1)<1e-13)})
  for z in out['constant1000K_mobile_cases']:
   checks.append({'name':name+' fraction '+str(z['diffusivity_multiplier']),'passed':bool(0<=z['mobile_fraction_released_at_next_same_stage_ignition']<=1)})
  checks.append({'name':name+' monotone D','passed':all(a['mobile_fraction_released_at_next_same_stage_ignition']<=b['mobile_fraction_released_at_next_same_stage_ignition'] for a,b in zip(out['constant1000K_mobile_cases'],out['constant1000K_mobile_cases'][1:]))})
 assert all(c['passed'] for c in checks)
 result=dict(scope='Exact source spectra + NIST mean projected-depth proxies. Conditional mobile graphite reference at1000K; D multipliers are scenarios. Strip sides perfectly absorbing up tofront; buried cap adds losses and invalidates direct transfer of strip efficiency.',
   surface_spectra_sha256=sha(sf),surface_stopping_sha256=sha(df),code_sha256=sha(__file__),kernel_sha256=sha(P/'helium_release.py'),stages=stages,checks=checks)
 (P/'helium-spectrum.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'checks':len(checks),'B1':stages['B1']},indent=2))
if __name__=='__main__':main()
