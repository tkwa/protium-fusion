"""Conditional geometry scan; every case is constrained by the same buckling screen."""
from bundle_solenoid import *
from time import monotonic
import gzip
start=monotonic();rows=[]
# Geometry is scaled to reduce finite-current-density and shield overhead.
R=1000
for sh,sc in [(4e9,2e9),(8e9,4e9),(8e9,8e9),(16e9,8e9),(16e9,16e9)]:
 for B in [40,60,80,100,120,140,160,180,200]:
  w=winding(B,R+5,sigma_hoop=sh)
  ac=w['field_energy_J_per_m']/sc
  rc=sqrt(w['outer_radius_m']**2+ac/pi)
  ioa=(rc*rc+w['outer_radius_m']**2)/4
  # 20% margin above the factor-of-two conventional Euler criterion.
  aspect=pi*sqrt(1e12*ioa/(2*sc*1.2))/R
  for heat in [.5,.6,.7]:
   for T in [1.,1.3,1.6,1.9,2.2]:
    r=plant(core_radius=R,aspect=aspect,B=B,T=T,heat=heat,sigma_hoop=sh,sigma_comp=sc)
    rows.append(r)
  print('done',sh/1e9,sc/1e9,B,'sec',round(monotonic()-start,1),flush=True)
rows.sort(key=lambda r:r['net_W_kg'],reverse=True)
(ROOT/'bundle-scan.json.gz').write_bytes(gzip.compress(json.dumps({'scope':'Fixed radial-profile family; all sources remain conditional, 1000m radius, geometry at conventional Euler-screen margin1.2.','cases':rows},separators=(',',':'),allow_nan=False).encode(),mtime=0))
for sh,sc in [(4e9,2e9),(8e9,4e9),(8e9,8e9),(16e9,8e9),(16e9,16e9)]:
 for heat in [.5,.6,.7]:
  good=[r for r in rows if r['assumptions']['sigma_hoop_Pa']==sh and r['assumptions']['sigma_comp_Pa']==sc and r['assumptions']['heat_efficiency']==heat and r['buckling_checks_pass']]
  if not good:print('NONE',sh,sc,heat);continue
  r=good[0];print('BEST',sh/1e9,sc/1e9,heat,'B',r['assumptions']['nominal_B_T'],'T',r['assumptions']['T_MeV'],'mW',r['net_W_kg']*1e3,'busmW',r['plant_bus_net_W_kg']*1e3,'length',r['length_m'],'mass',r['total_mass_kg'],'P',r['net_W'],'local',r['local_buckling_margin'],flush=True)
