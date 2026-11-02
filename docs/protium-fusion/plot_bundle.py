"""Scientific summary of conditional shared-solenoid allocations."""
from pathlib import Path
import json,gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parent
rows=json.loads(gzip.decompress((ROOT/'bundle-scan.json.gz').read_bytes()))['cases']
sensitivity=json.loads((ROOT/'bundle-sensitivity.json').read_text());ref=sensitivity['reference']
endrows=json.loads((ROOT/'end-transfer-sensitivity.json').read_text())['cases']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(13.8,9.4))
fig.subplots_adjust(left=.09,right=.97,bottom=.17,top=.87,hspace=.42,wspace=.34)
fig.suptitle('Shared-magnet power depends on materials, recovery and end transfer',x=.08,y=.97,ha='left',fontsize=16,fontweight='bold')
fig.text(.08,.925,'Titanium-49 cycle · 1.6 MeV ions · future material and recovery assumptions · no engineered reactor demonstrated',fontsize=11,color='#414d5a')
colors=['#9b653a','#447a9b','#79529a']
ax=axs[0,0]
for (sh,sc),col in zip([(4e9,2e9),(8e9,4e9),(16e9,8e9)],colors):
 data=sorted([r for r in rows if r['assumptions']['sigma_hoop_Pa']==sh and r['assumptions']['sigma_comp_Pa']==sc and r['assumptions']['heat_efficiency']==.6 and r['assumptions']['T_MeV']==1.6],key=lambda r:r['assumptions']['nominal_B_T'])
 x=[r['assumptions']['nominal_B_T'] for r in data]
 ax.plot(x,[1000*r['net_W_kg'] for r in data],'-o',ms=3,color=col,label=f'{sh/1e9:g} / {sc/1e9:g} GPa hoop / compression')
er=sorted([r for r in endrows if r['assumptions']['heat_efficiency']==.6 and r['assumptions']['end_transfer_in']==.999],key=lambda r:r['assumptions']['nominal_B_T'])
ax.plot([r['assumptions']['nominal_B_T'] for r in er],[r['net_W_kg']*1000 for r in er],'--',lw=2,color='#333',label='16 / 8 GPa + explicit end transfer')
ax.axhline(1,color='#333',lw=1,ls='--');ax.axhline(0,color='#aaa',lw=.7)
ax.set(title='A  Base allocation and added end transfer',xlabel='Nominal central field (T)',ylabel='Net useful power / total mass (mW/kg)',ylim=(-.4,1.85))
ax.legend(frameon=False,fontsize=8.6,loc='upper left');ax.grid(alpha=.15)
ax=axs[0,1]
mb=ref['mass_breakdown_kg']
groups={'Axial compression':mb['axial_compression_support'],'Hoop structure':mb['circumferential_structure'],'End allowance':mb['end_system_allowance'],'Radiators':mb['radiators'],'All other systems':sum(v for k,v in mb.items() if k not in ['axial_compression_support','circumferential_structure','end_system_allowance','radiators'])}
labels=list(groups);values=np.array(list(groups.values()))/ref['total_mass_kg']*100
ax.barh(labels[::-1],values[::-1],color=['#8eabb5','#5e8198','#a08bb5','#79529a','#bdcbd1'][::-1]);ax.set_xlim(0,65)
for i,v in enumerate(values[::-1]):ax.text(v+.7,i,f'{v:.1f}%',va='center',fontsize=9)
ax.set(title='B  Mass allocation in the 120 T base case',xlabel='Share of complete allocated mass (%)')
ax=axs[1,0]
for b,col in zip([100,120,180],['#447a9b','#79529a','#b76c3c']):
 data=sorted([r for r in endrows if r['assumptions']['nominal_B_T']==b and r['assumptions']['end_transfer_in']==.999],key=lambda r:r['assumptions']['heat_efficiency'])
 ax.plot([r['assumptions']['heat_efficiency']*100 for r in data],[r['net_W_kg']*1000 for r in data],'-o',color=col,label=f'{b} T')
ax.axhline(1,color='#333',lw=1,ls='--');ax.set(title='C  Explicit end transfer at 99.9% each way',xlabel='Thermal-to-electric conversion (%)',ylabel='After replacement manufacturing (mW/kg)',ylim=(0,1.85));ax.legend(frameon=False);ax.grid(alpha=.15)
ax=axs[1,1]
keys=['side_cost','heater','potential_in','extraction_MeV_per_pn','annual_replacement_fraction']
vals=[230,.9,.99,1.,.001]
labels=['230 MeV side restoration','90% ion heating','99% potential injection','1 MeV extraction loss','0.1% replacement / year']
out=[]
for k,v in zip(keys,vals):
 rr=next(r['result'] for r in sensitivity['cases'] if r['parameter']==k and r['value']==v);out.append(rr['net_W_kg']*1000)
ax.barh(labels[::-1],out[::-1],color=['#9e5551' if v<0 else '#77899b' for v in out[::-1]])
ax.axvline(ref['net_W_kg']*1000,color='#79529a',ls='--',label=f"Reference {ref['net_W_kg']*1000:.2f}")
ax.axvline(0,color='#aaa',lw=.7);ax.set(title='D  One changed assumption at a time',xlabel='After replacement manufacturing (mW/kg)',xlim=(-2.1,1.35));ax.legend(frameon=False,loc='upper left',fontsize=9)
fig.text(.08,.095,'Solid traces in A include a 0.2 MeV extraction allowance; the dashed trace and panel C add a 25 kT end barrier and its conversion equipment.',fontsize=9)
fig.text(.08,.071,'Default: 95% ion heating/direct conversion; 99.9% potential recovery each way; 50 MeV per side event; 30 K radiators; 100 MJ/kg embodied energy;',fontsize=9)
fig.text(.08,.047,'0.01% annual replacement. All traces include complete mass allowances. Recycling, extraction hardware, shielding and conductor life remain unvalidated.',fontsize=9)
for fmt in ['png','svg']:
 p=ROOT/f'bundle-solenoid-summary.{fmt}';fig.savefig(p,dpi=180)
 if fmt=='svg':p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
print(ROOT/'bundle-solenoid-summary.png')
