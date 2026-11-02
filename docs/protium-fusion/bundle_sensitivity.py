"""One-at-a-time assumptions around a conditional 120T shared-field allocation."""
from bundle_solenoid import *
base=dict(core_radius=1000,B=120,T=1.6,heat=.6,sigma_hoop=16e9,sigma_comp=8e9)
w=winding(base['B'],1005,sigma_hoop=16e9);ac=w['field_energy_J_per_m']/8e9
rc=sqrt(w['outer_radius_m']**2+ac/pi)
base['aspect']=pi*sqrt(1e12*(rc*rc+w['outer_radius_m']**2)/4/(2*8e9*1.2))/1000
sweeps={'side_cost':[0,10,50,100,230], 'heater':[.90,.93,.95,.97,.99],
 'direct':[.85,.90,.95,.98,.99], 'potential_in':[.99,.995,.999,.9999,1.],
 'potential_out':[.99,.995,.999,.9999,1.], 'heat':[.4,.5,.6,.7,.8],
 'extraction_MeV_per_pn':[0,.2,.5,1.,2.], 'conversion_kg_W':[.1,1.,10.],
 'gamma_leak':[1e-6,1e-5,1e-4,1e-3], 'neutron_leak':[1e-8,1e-6,1e-4,1e-3],
 'annual_replacement_fraction':[0,1e-5,1e-4,1e-3,.01], 'embodied_J_kg':[1e7,1e8,1e9],
 'J':[1e7,1e8,1e9], 'end_system_fraction':[.15,.3,.5,1.],
 'radiator_areal_mass':[.1,1.,10.,100.], 'neutron_capture':[.95,.98,.995,.999]}
rows=[]
for key,values in sweeps.items():
 for value in values:
  args=dict(base);args[key]=value
  try:r=plant(**args)
  except ValueError as e:r={'error':str(e)}
  rows.append({'parameter':key,'value':value,'result':r})
 out=[(x['value'],round(x['result'].get('net_W_kg',float('nan'))*1e3,4)) for x in rows if x['parameter']==key]
 print(key,out,flush=True)
(ROOT/'bundle-sensitivity.json').write_text(json.dumps({'scope':'One-at-a-time sensitivity, retaining the specified radial profile; no feasibility established.','base_arguments':base,'reference':plant(**base),'cases':rows},indent=2,allow_nan=False)+'\n')
