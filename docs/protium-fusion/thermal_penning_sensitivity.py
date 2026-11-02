from penning_equilibrium import *
from scipy.optimize import minimize_scalar
best=json.loads((ROOT/'penning-scan.json').read_text())['best'][0]
eq=equilibrium(best['h0'],best['c0'],best['xi'],60,nout=4001)
rows=[]
for T in [.4,.5,.6,.7,.8,1.,1.2]:
 for B in [50,100,200,400]:
  for pot in [1.,.999,.99,.95]:
   fn=lambda x:-at_wall(eq,x,T=T,B=B,potential_in_eff=pot,potential_out_eff=pot,radiation_escape_fraction=1.)['field_mass_only_net_W_kg']
   sol=minimize_scalar(fn,bounds=(12,60),method='bounded',options={'xatol':.005})
   rows.append(at_wall(eq,sol.x,T=T,B=B,potential_in_eff=pot,potential_out_eff=pot,radiation_escape_fraction=1.))
out={'scope':'Same selected dimensionless equilibrium, temperature/field/potential-recovery sensitivity. Full vacuum Larmor loss,95% heater/ion recovery,40% heat,50MeV side cost; field mass only. 200/400T are extrapolations without a coil design.','cases':rows}
(ROOT/'thermal-penning-sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
for r in rows:
 if r['B_T']==100 and r['potential_in_efficiency'] in [1.,.99]:
  print('T',r['T_MeV'],'pot',r['potential_in_efficiency'],'x',round(r['xwall'],3),'net',round(r['net_MeV_per_pn'],3),'W/kg',f"{r['field_mass_only_net_W_kg']:.3g}",'potMeV',round(r['electrostatic_work_circulated_MeV_per_pn'],2),'rad',round(r['vacuum_ion_Larmor_MeV_per_pn'],3))
