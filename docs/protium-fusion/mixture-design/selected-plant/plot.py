"""Standalone scientific figure from the selected source-family plant accounts."""
from pathlib import Path
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
trade=json.loads((HERE/'receiver-tradeoff.json').read_text());service=json.loads((HERE/'services-sensitivity.json').read_text());point=json.loads((HERE/'point.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
colors=['#176789','#b55b28','#6d7786']
series=[('EM20',0,'Coupled electromagnetic transport'),('all_to_cold',0,'All saved prompt heat sent cold'),('EM20',.005,'Additional 0.5% sent cold')]
rows=[]
for (policy,extra,label),color in zip(series,colors):
 data=sorted([r for r in trade['rows'] if r['prompt_policy']==policy and r['missed_hot_routing']=='cold_equipment' and r['additional_missed_hot_fraction_HplusUprecursor']==extra],key=lambda r:r['receiver_sigma_kg_m2'])
 axes[0].plot([r['receiver_sigma_kg_m2'] for r in data],[r['mW_kg'] for r in data],'-o',ms=4,lw=1.5,color=color,label=label)
 rows += [{'series':label,'receiver_kg_m2':r['receiver_sigma_kg_m2'],'service_kg_m2':1000,'net_mW_kg':r['mW_kg']} for r in data]
axes[0].axvspan(250,350,color='#c5dce5',alpha=.35,zorder=0)
axes[0].set_xscale('log');axes[0].set_xticks([200,350,500,1000,2000],labels=['200','350','500','1,000','2,000'])
axes[0].set(xlabel='Hot graphite receiver (kg/m²)',ylabel='Net specific power (mW/kg)',title='Receiver mass and retained heat')
axes[0].legend(fontsize=8,loc='lower left',frameon=False)
data=service['rows'];axes[1].plot([r['service_kg_m2'] for r in data],[r['mW_kg'] for r in data],'-o',color=colors[0],ms=5)
axes[1].scatter([1000],[1000*point['plant']['net_W_kg']],s=65,facecolors='white',edgecolors=colors[0],zorder=3)
axes[1].set(xlabel='General service allocation (kg/m²)',ylabel='Net specific power (mW/kg)',title='Equipment mass remains consequential')
axes[1].text(.03,.06,'350 kg/m² hot receiver\nLower allocations need a completed equipment bill.',transform=axes[1].transAxes,fontsize=8,ha='left')
rows += [{'series':'Service allocation sensitivity','receiver_kg_m2':350,'service_kg_m2':r['service_kg_m2'],'net_mW_kg':r['mW_kg']} for r in data]
for ax in axes:
 ax.grid(axis='y',alpha=.18);ax.spines[['top','right']].set_visible(False)
fig.suptitle('Staged ordinary-hydrogen fusion: conditional plant comparison',fontsize=13)
for ext in ['png','pdf','svg']:fig.savefig(HERE/f'plant-tradeoffs.{ext}',dpi=190)
with (HERE/'plant-tradeoffs.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(HERE/'plot-provenance.json').write_text(json.dumps({'inputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'receiver-tradeoff.json',HERE/'services-sensitivity.json',HERE/'point.json',Path(__file__)]},'scope':'Lines connect calculated scenarios; they are not an interpolated engineering qualification. Fixed radius1.1e9m,200K radiators,60% driver/heat efficiencies. All missed hot energy is charged to cold equipment. Source assembly and hardware assumptions remain conditional.','rows':len(rows)},indent=2)+'\n')
print('Saved PNG/PDF/SVG, CSV and provenance.')
