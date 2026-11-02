from penning_equilibrium import *
from scipy.optimize import minimize_scalar
best=json.loads((ROOT/'penning-scan.json').read_text())['best'][0]
eq=equilibrium(best['h0'],best['c0'],best['xi'],60,nout=4001)
rows=[]
for th in [.4,.6,.7,.8]:
 for side in [0.,10.,50.,230.]:
  cases=[]
  for T in np.linspace(.4,2.,33):
   fn=lambda x:-at_wall(eq,x,T=float(T),thermal_eff=th,side_cost=side,potential_in_eff=.999,potential_out_eff=.999,radiation_escape_fraction=1.)['field_mass_only_net_W_kg']
   sol=minimize_scalar(fn,bounds=(12,60),method='bounded',options={'xatol':.005})
   row=at_wall(eq,sol.x,T=float(T),thermal_eff=th,side_cost=side,potential_in_eff=.999,potential_out_eff=.999,radiation_escape_fraction=1.)
   row['heat_conversion_efficiency']=th;row['net_side_restoration_cost_MeV']=side;cases.append(row)
  high=max(cases,key=lambda a:a['field_mass_only_net_W_kg'])
  rows.append({'heat_efficiency':th,'side_cost_MeV':side,'best':high,'cases':cases})
  print('heat',th,'side',side,'T',round(high['T_MeV'],2),'net',round(high['net_MeV_per_pn'],3),'mW/kg',round(high['field_mass_only_net_W_kg']*1000,3),flush=True)
(ROOT/'thermal-recovery-sensitivity.json').write_text(json.dumps({'scope':'Specified 100 T equilibrium,95% heater/ion recovery,99.9% potential transfer, full vacuum ion radiation. Heat conversion and all-inclusive side-restoration cost are assumptions. No complete plant or converter.','scenarios':rows},indent=2,allow_nan=False)+'\n')
