"""Static scientific figure for the explicit catalyst-network allocation."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parent
rows=json.loads((ROOT/'network-plant.json').read_text())['cases'];sens=json.loads((ROOT/'network-plant-sensitivity.json').read_text())
ref=next(r for r in rows if r['assumptions']['B']==140 and r['assumptions']['heat']==.7 and r['assumptions']['end_in']==.999)
plt.rcParams.update({'svg.hashsalt':'protium-network2026','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
ax=axes[0,0]
for heat,color in [(.6,'#b45f06'),(.7,'#176b91')]:
 for eta,style in [(.999,'-'),(.9999,'--')]:
  ss=[r for r in rows if r['assumptions']['heat']==heat and r['assumptions']['end_in']==eta]
  ax.plot([r['assumptions']['B'] for r in ss],[r['net_W_kg']*1e3 for r in ss],style,marker='o',color=color,label=f'{heat:.0%} heat, {eta:.2%} end transfer')
ax.axhline(1,color='#444',lw=1,ls=':');ax.set(xlabel='Nominal field (T)',ylabel='Net useful output (mW/kg)',title='Complete allocated plant, including replacement');ax.legend(fontsize=8)
ax=axes[0,1]
labels=[];vals=[]
for r in ref['channel_allocations']:
 if r['volume_fraction_of_channels']>.005:
  z,a=r['ZA'];labels.append({22:'Ti',23:'V',24:'Cr',25:'Mn'}[z]+'-'+str(a));vals.append(100*r['volume_fraction_of_channels'])
labels.append('Other');vals.append(100-sum(vals))
ax.barh(labels[::-1],vals[::-1],color='#368783');ax.set(xlabel='Share of active channel volume (%)',title='Auxiliary processors use 41% of the channels',xlim=(0,65))
for i,v in enumerate(vals[::-1]):ax.text(v+.7,i,f'{v:.1f}',va='center',fontsize=9)
ax=axes[1,0]
m=ref['mass_breakdown_kg'];groups=[('Axial support',m['axial_compression_support']),('Hoop support',m['circumferential_structure']),('Radiators',m['radiators']),('End systems',m['end_system_allowance']),('Conversion',m['conversion_and_potential_handling'])]
groups.append(('Other',ref['total_mass_kg']-sum(v for _,v in groups)))
ax.barh([k for k,v in groups][::-1],[v/ref['total_mass_kg']*100 for k,v in groups][::-1],color='#6b709c');ax.set(xlabel='Share of installed mass (%)',title='Reference mass: 1.47 quadrillion kg')
ax=axes[1,1]
show=[sens['cases'][i] for i in [0,2,3,4,6,7]]
labels=['Reference','90% ion heating','90% direct conversion','99.8% end transfer','0.03% replacement/year','Extra 1 MeV/primary']
ax.barh(labels[::-1],[r['net_W_kg']*1e3 for r in show][::-1],color=['#8a6a48']*5+['#176b91']);ax.axvline(1,color='#444',ls=':',lw=1);ax.set(xlabel='Net useful output (mW/kg)',title='One assumption changed at a time')
fig.suptitle('Titanium recycling preserves a conditional route above 1 mW/kg',fontsize=15)
fig.supxlabel('T = 1.6 MeV; reference 140 T, 16 GPa hoop / 8 GPa compression. Confinement, materials and transport remain unvalidated.',fontsize=9)
fig.savefig(ROOT/'network-plant-summary.png',dpi=180);fig.savefig(ROOT/'network-plant-summary.svg',metadata={'Date':None})
print('saved figure')
