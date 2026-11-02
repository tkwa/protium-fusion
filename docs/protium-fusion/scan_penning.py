from penning_equilibrium import *
from scipy.optimize import minimize_scalar
import time
start=time.monotonic();out=[];failed=0
for xi in [.1,.2,.3,.4,.5,.6,.7,.8,.9]:
 for frac in [.25,.5,.75,1,1.1,1.2,1.35,1.5]:
  S=2*KAPPA*xi*(1-xi)*frac
  for ph in [.1,.25,.5,.75,.9]:
   eq=equilibrium(S*ph,S*(1-ph)/6,xi,60,rtol=1e-7,nout=2001)
   if not eq['success']:
    failed+=1;continue
   def loss(x):return -at_wall(eq,x)['field_mass_only_net_W_kg']
   opt=minimize_scalar(loss,bounds=[12,60],method='bounded',options={'xatol':.02})
   row=at_wall(eq,opt.x);row['charge_density_fraction']=frac;row['proton_charge_fraction']=ph
   out.append(row)
 print('xi',xi,'sec',round(time.monotonic()-start,1),'success',len(out),'failed',failed,flush=True)
out.sort(key=lambda a:a['field_mass_only_net_W_kg'],reverse=True)
(ROOT/'penning-scan.json').write_text(json.dumps({'scope':'Coarse grid plus scalar wall-radius optimization, specified .6 MeV,100 T,95/95/40 efficiency,50 MeV side cost; field only.','successful':len(out),'failed_profile':failed,'elapsed_s':time.monotonic()-start,'best':out[:30]},indent=2,allow_nan=False)+'\n')
for r in out[:10]:print({k:round(r[k],6) for k in ['xi','charge_density_fraction','proton_charge_fraction','xwall','net_MeV_per_pn','field_mass_only_net_W_kg','mean_pair_rotation_MeV','axis_to_wall_potential_V','rotation_speed_wall_over_c']})
