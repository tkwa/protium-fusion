"""Research figure for the conditional heavy-catalyst comparison."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'protium-catalyst-v1'})
cat=json.loads((ROOT/'catalyst-catalog.json').read_text())['candidates']
sens=json.loads((ROOT/'catalyst-sensitivity.json').read_text())
names=['49Ti','51V','53Cr','55Mn'];colors=['#126782','#b66a18','#734c9e','#34754b']
fig,ax=plt.subplots(2,2,figsize=(12.2,10.0));fig.subplots_adjust(top=.87,bottom=.18,hspace=.42,wspace=.30)
fig.suptitle('Heavy catalysts for an electron-free hydrogen reactor',fontsize=18,x=.07,ha='left',y=.97)
fig.text(.07,.925,'Calculated screens from PADF-2 + NUBASE2020. Power figures count allocated field structure only.',fontsize=11)
for n,c in zip(names,colors):
 r=next(r for r in cat if r['target']==n);rr=[x for x in r['thermal_cases'] if .4<=x['ion_T_MeV']<=2.5]
 ax[0,0].plot([x['ion_T_MeV'] for x in rr],[x['pn_ground_sigma_v_m3_s'] for x in rr],color=c,label=n,lw=2)
 rr=[x for x in sens['cases'] if x['target']==n and x['assumptions']['heat']==.4 and x['assumptions']['net_side_restoration_MeV']==50]
 ax[0,1].plot([x['T_MeV'] for x in rr],[x['field_mass_only_net_W_kg']*1e3 if x['net_MeV_per_pn']>0 else np.nan for x in rr],color=c,label=n,lw=2)
ax[0,0].set(title='A  Useful proton-to-neutron reaction rate',xlabel='Common ion temperature (MeV)',ylabel=r'$\langle\sigma v\rangle$ (m³/s)',yscale='log');ax[0,0].legend(ncol=2,frameon=False)
ax[0,1].set(title='B  Four catalysts under the same assumptions',xlabel='Common ion temperature (MeV)',ylabel='Conditional field output (mW/kg)');ax[0,1].axhline(0,color='#aaa',lw=.8);ax[0,1].set_ylim(0,.22)
for heat,c in [(.4,'#126782'),(.6,'#c05a42')]:
 for side,ls in [(10.,':'),(50.,'-'),(230.,'--')]:
  rr=[x for x in sens['cases'] if x['target']=='49Ti' and x['assumptions']['heat']==heat and x['assumptions']['net_side_restoration_MeV']==side]
  ax[1,0].plot([x['T_MeV'] for x in rr],[x['field_mass_only_net_W_kg']*1e3 if x['net_MeV_per_pn']>0 else np.nan for x in rr],c=c,ls=ls,lw=2,label=f'{heat:.0%} heat, {side:g} MeV side cost')
ax[1,0].set(title='C  Titanium: recovery and restoration costs',xlabel='Common ion temperature (MeV)',ylabel='Conditional field output (mW/kg)',ylim=(0,.65));ax[1,0].axhline(0,color='#aaa',lw=.8);ax[1,0].legend(frameon=False,fontsize=8,ncol=2,loc='upper left')
for i,n in enumerate(names):
 b=next(x for x in sens['best'] if x['target']==n and x['assumptions']['heat']==.4 and x['assumptions']['net_side_restoration_MeV']==50)
 power=b['reaction_rate_per_m_s']*b['net_MeV_per_pn']*1.602176634e-13
 ax[1,1].scatter(i,b['decay_inventory_per_m_kg']/power*1e9,s=55,c=colors[i],marker='o',label='Cold decay inventory' if i==0 else None)
 ax[1,1].scatter(i,b['field_mass_per_m_kg']/power*1e9,s=55,c=colors[i],marker='s',label='Allocated field structure' if i==0 else None)
ax[1,1].set(title='D  Inventory versus field mass at best temperature',ylabel='Mass per GW of conditional output (kg)',yscale='log',xticks=range(4),xticklabels=names);ax[1,1].legend(frameon=False,fontsize=9,loc='center left')
for a in ax.flat:a.grid(alpha=.16)
fig.text(.07,.035,'100 T; 1 MJ/kg structural specific strength; field mass factor 5; 70% packing; 95% ion heating and recovery;\n99.9% potential-energy recovery in each direction. Finite-length confinement, coils, cooling, neutron collection and\nclosed side-reaction processing are not engineered. The whole-plant target is 1 mW/kg.\nCurves show positive model output; negative cases are omitted.',fontsize=9,color='#444',linespacing=1.5)
fig.savefig(ROOT/'catalyst-comparison.png',dpi=170)
fig.savefig(ROOT/'catalyst-comparison.svg',metadata={'Date':None})

# Normalize generated SVG whitespace so repository checks remain reproducible.
svg=ROOT/'catalyst-comparison.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
