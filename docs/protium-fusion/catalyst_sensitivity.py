from catalyst_penning import *
rows=[]
for target in ['49Ti','51V','53Cr','55Mn']:
 b=json.loads((ROOT/f'catalyst-penning-{target}.json').read_text())['best'][0]
 z,a=CAT[target]['Z'],CAT[target]['A'];eq=equilibrium(z,a,b['h0'],b['c0'],b['xi'],nout=4001,rtol=1e-10)
 for heat in [.4,.6]:
  for side in [10.,50.,230.]:
   for T in np.arange(.5,2.51,.1):
    T=round(float(T),4)
    fn=lambda x:-at_wall(eq,x,target,T=T,heat=heat,side_cost=side)['field_mass_only_net_W_kg']
    op=minimize_scalar(fn,bounds=(10,60),method='bounded',options={'xatol':.005})
    rows.append(at_wall(eq,op.x,target,T=T,heat=heat,side_cost=side))
 best=max((r for r in rows if r['target']==target and r['assumptions']['heat']==.4 and r['assumptions']['net_side_restoration_MeV']==50),key=lambda r:r['field_mass_only_net_W_kg'])
 print(target,'best40pct',json.dumps({k:best[k] for k in ['T_MeV','net_MeV_per_pn','field_mass_only_net_W_kg','radius_m','pair_rotation_MeV','potential_work_MeV_per_pn','decay_inventory_per_m_kg']}),flush=True)
summary=[]
for target in ['49Ti','51V','53Cr','55Mn']:
 for heat in [.4,.6]:
  for side in [10.,50.,230.]:
   rr=[r for r in rows if r['target']==target and r['assumptions']['heat']==heat and r['assumptions']['net_side_restoration_MeV']==side]
   summary.append(max(rr,key=lambda r:r['field_mass_only_net_W_kg']))
out={'scope':'Fixed dimensionless shape selected at T=1MeV; thermal,heat recovery and side-cost sensitivity. All net values remain conditional on unmodeled engineering.100T,1MJ/kg specific strength,factor5 field allocation.','best':summary,'cases':rows}
(ROOT/'catalyst-sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print('Ti cases')
for r in summary[:6]:print('heat',r['assumptions']['heat'],'side',r['assumptions']['net_side_restoration_MeV'],'T',r['T_MeV'],'mW/kg',r['field_mass_only_net_W_kg']*1e3)
