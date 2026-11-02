"""Independent geometry, absorption and complete-ledger checks for the audit."""
from pathlib import Path
import json,gzip,math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from capture_competition import read_reference,competition,geometry,allocated_plant,mixture,ELEMENTS,COATINGS,ROOT,HERE,MEV,U,YEAR,SB
checks=[]
def check(name,ok,detail=None):
 checks.append({'name':name,'pass':bool(ok),'detail':detail})
 if not ok:print('FAIL',name,detail)
def close(name,a,b,rtol=1e-10,atol=1e-10):check(name,abs(a-b)<=atol+rtol*abs(b),{'actual':a,'expected':b,'relative_error':a/b-1 if b else None})
r=read_reference();base=competition(r);eta=base['hydrogen_capture_fraction_after_other_losses']
# Force the reference's total useful capture while preserving its original material.
c=competition(r,other_neutron_loss=1-.995/eta);p=allocated_plant(r,c)
for key in ['net_W','total_mass_kg','net_W_kg','plant_bus_net_W','capture_penalty_W','nuclear_after_capture_W','waste_heat_W','cold_radiation_W','refrigeration_W']:
 close('saved reference recovery '+key,p[key],r[key],rtol=3e-11)
# Cylindrical volume integrations do not use the implementation's polynomial.
for t in [1e-7,1e-5,1e-4]:
 g=geometry(r,t,'Ti30Zr30V40');coat=0.;carbon=0.;water_displaced=0.
 for row in r['channel_allocations']:
  rr,tw,N=row['radius_m'],row['wall_m'],row['channels'];ll=r['active_length_m'];newtw=tw*(rr+t)/rr
  coat+=quad(lambda x:2*math.pi*x,rr,rr+t,epsabs=1e-18)[0]*N*ll
  carbon+=quad(lambda x:2*math.pi*x,rr+t,rr+t+newtw,epsabs=1e-18)[0]*N*ll
  water_displaced+=quad(lambda x:2*math.pi*x,rr+tw,rr+t+newtw,epsabs=1e-18)[0]*N*ll
 close('coating integral '+str(t),g['coating_volume_m3'],coat,rtol=1e-9)
 close('backing integral '+str(t),g['carbon_wall_mass_kg'],carbon*r['assumptions']['wall_density'],rtol=1e-10)
 close('water displacement '+str(t),g['displaced_water_volume_m3'],water_displaced,rtol=1e-9)
# Direct mass-fraction calculation of each film's N sigma.
for name in COATINGS:
 c=competition(r,1e-6,name);m=c['geometry']['coating_mass_kg'];mix=mixture(name)
 direct=sum(m*f/(ELEMENTS[e]['A']*U)*ELEMENTS[e]['abs_barn']*1e-28 for e,f in mix['mass_fractions'].items())
 close('alloy independent atom count '+name,c['absorption_inventory_m2']['coating'],direct)
 close('absorption probability sum '+name,sum(c['absorption_fractions'].values()),1.)
# Recover the H number inventory from known water atomic mass independently.
from catalyst_catalog import atomicmass
from nonneutral_screen import C,Q
h_atoms=2*r['water_volume_m3']*r['assumptions']['warm_water_density']/((atomicmass(8,16)+2*atomicmass(1,1))*MEV/C**2)
close('H atom number',base['absorption_inventory_m2']['hydrogen'],h_atoms*.3326e-28)
# Full reference implementation with different useful capture fraction.
from network_plant import plant
energy=json.loads((ROOT/'network-ion-budget.json').read_text());family=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
for capture in [.997,.95,.9]:
 c=competition(r);c['hydrogen_capture_fraction_after_other_losses']=capture
 ours=allocated_plant(r,c);direct=plant(energy,family,neutron_capture=capture)
 for key in ['net_W','total_mass_kg','nuclear_after_capture_W','cold_radiation_W','refrigeration_W','plant_bus_net_W']:
  close('original plant independent capture '+str(capture)+' '+key,ours[key],direct[key],rtol=3e-10)
# Numerical mass feedback root, independent of its algebraic solution.
result=json.loads((HERE/'capture-competition.json').read_text())
for c in result['cases']+result['impurity_cases']:
 p=c['plant'];a=r['assumptions'];fixed=sum(v for k,v in p['mass_breakdown_kg'].items() if k not in ['radiators','replacement_manufacturing_equipment','replacement_inventory'])
 radkg=a['radiator_areal_mass']/(a['radiator_emissivity']*SB*(a['radiator_K']**4-a['ambient_K']**4));k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR
 def f(m):
  net=p['plant_bus_net_W']-k*m
  return fixed+(p['nuclear_after_capture_W']-net)*radkg+a['manufacturing_kg_W']*k*m+a['annual_replacement_fraction']*a['replacement_inventory_years']*m-m
 sol=brentq(f,.1*r['total_mass_kg'],3*r['total_mass_kg'])
 label=c['geometry']['coating']+'/'+str(c['geometry']['thickness_m'])+'/'+str(c['water_boron_ppm'])+'/'+str(c['wall_boron_ppm'])+'/'+str(c['wall_nitrogen_ppm'])
 close('numerical maintenance root '+label,p['total_mass_kg'],sol)
 close('mass ledger '+label,p['mass_ledger_relative_error'],0.,atol=2e-15)
 close('energy ledger '+label,p['energy_ledger_relative_error'],0.,atol=2e-15)
# Limits independently evaluated after reading the output.
for d in result['coating_tolerances']:
 for label in ['capture_99_5','specific_power_1mW']:
  t=d[label+'_thickness_m']
  if t is None:continue
  c=competition(r,t,d['coating'])
  got=c['hydrogen_capture_fraction_after_other_losses'] if label=='capture_99_5' else allocated_plant(r,c)['net_W_kg']
  close('coating limit '+d['coating']+'/'+label,got,.995 if label=='capture_99_5' else .001,rtol=1e-8)
for d in result['impurity_limits']:
 for label in ['capture_99_5','specific_power_1mW']:
  c=competition(r,**{d['impurity']:d[label+'_ppm']})
  got=c['hydrogen_capture_fraction_after_other_losses'] if label=='capture_99_5' else allocated_plant(r,c)['net_W_kg']
  close('impurity limit '+d['impurity']+'/'+label,got,.995 if label=='capture_99_5' else .001,rtol=1e-8)
report={'checks':checks,'count':len(checks),'passed':sum(x['pass'] for x in checks),'all_pass':all(x['pass'] for x in checks)}
(HERE/'capture-competition-validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print({k:v for k,v in report.items() if k!='checks'})
if not report['all_pass']:raise SystemExit(1)
