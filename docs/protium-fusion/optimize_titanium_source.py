"""Refine Ti49 radial equilibria at the higher heat-conversion operating point."""
from catalyst_penning import *
from time import monotonic
start=monotonic();z,a=22,49;M=nucmass(z,a)/nucmass(1,1);S=z/(2*M)
rows=[];bad=0
for xi in [.35,.4,.45,.5,.55,.6,.65,.7]:
 for fraction in [.5,.7,.85,.95,1.,1.05,1.1,1.15,1.2]:
  for share in [.4,.5,.6,.7,.8,.9,.95]:
   eq=equilibrium(z,a,S*fraction*share,S*fraction*(1-share)/z,xi)
   if not eq['success']:bad+=1;continue
   fn=lambda x:-at_wall(eq,x,'49Ti',T=1.7,heat=.6)['field_mass_only_net_W_kg']
   opt=minimize_scalar(fn,bounds=(10,60),method='bounded',options={'xatol':.005})
   r=at_wall(eq,opt.x,'49Ti',T=1.7,heat=.6);r['central_charge_fraction']=fraction;r['proton_charge_share']=share;rows.append(r)
rows.sort(key=lambda r:r['field_mass_only_net_W_kg'],reverse=True)
out={'scope':'Refined conditional Ti49 source equilibrium at1.7MeV,60% heat; all original source assumptions remain.','attempted':len(rows)+bad,'density_blowups':bad,'best':rows[:40]}
(ROOT/'titanium-source-refined.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print('time',monotonic()-start,'attempted',len(rows)+bad);print(json.dumps(rows[0],indent=2))
