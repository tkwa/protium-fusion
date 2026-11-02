"""Explicit high-voltage end-transfer branch; neither extraction path is engineered."""
from bundle_solenoid import *
rows=[]
for B in [100,120,140,160,180]:
 w=winding(B,1005,sigma_hoop=16e9);ac=w['field_energy_J_per_m']/8e9
 rc=sqrt(w['outer_radius_m']**2+ac/pi)
 aspect=pi*sqrt(1e12*(rc*rc+w['outer_radius_m']**2)/4/(2*8e9*1.2))/1000
 for heat in [.6,.7]:
  for eta in [.99,.999,.9999,.99999]:
   r=plant(core_radius=1000,aspect=aspect,B=B,T=1.6,heat=heat,sigma_hoop=16e9,sigma_comp=8e9,end_barrier_kT=25,end_transfer_in=eta,end_transfer_out=eta)
   rows.append(r)
   print(B,heat,eta,'endMeV',r['additional_end_transfer_loss_MeV_per_pn'],'mWkg',r['net_W_kg']*1e3,'bus',r['plant_bus_net_W_kg']*1e3,flush=True)
(ROOT/'end-transfer-sensitivity.json').write_text(json.dumps({'scope':'Conditional full-mass branch adding920MeV of coherent end-barrier circulation per pn. Specified recovery, no RF/electrode/selection design.','cases':rows},indent=2,allow_nan=False)+'\n')
