"""Standalone scientific figure for the current catcher follow-up."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
d=read(P/'stopping/surface-stopping.json');thermal=read(P/'geometry/temperature-histories.json');ref=read(P/'refurbishment/results.json')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,3,figsize=(16,4.7),layout='constrained')
colors=['#377eb8','#e41a1c','#4daf4a','#984ea3']
for stage,c in zip(['A','B1','B2','C13'],colors):
 s=d['stages'][stage]['species']['p'];z=np.array(s['packet_mean_projected_depth_um']);w=np.array(s['packet_nuclei_m2']);k=np.argsort(z)
 axes[0].plot(z[k],np.cumsum(w[k])/w.sum(),label=f'{stage} H',color=c,lw=1.7)
s=d['stages']['B1']['species']['he4'];z=np.array(s['packet_mean_projected_depth_um']);w=np.array(s['packet_nuclei_m2']);k=np.argsort(z)
axes[0].plot(z[k],np.cumsum(w[k])/w.sum(),label='B1 He',color='#222222',lw=2,ls='--')
axes[0].axvline(150,color='#999999',ls=':',label='150 µm active layer')
axes[0].set(xscale='log',xlim=(.03,200),ylim=(0,1.03),xlabel='Mean projected stopping depth (µm)',ylabel='Cumulative nucleus fraction',title='Current energy spectra mapped to stopping depth')
axes[0].legend(fontsize=8,loc='lower right');axes[0].grid(alpha=.15)
p=max(thermal['profiles'],key=lambda x:x['maximum_skin_front_K']);r=list(csv.DictReader(Path(p['csv']).open()))
t=np.array([float(x['relative_to_prompt_arrival_s']) for x in r]);T=np.array([float(x['front_bulk_k50_K']) for x in r]);q=np.array([float(x['material_W_m2'])+float(x['delayed_decay_W_m2']) for x in r])
for k,c in [(4.,'#377eb8'),(1.,'#ff7f00'),(.5,'#e41a1c')]:
 curve=T+q*(150e-6/k+50e-6/(50*.25))
 axes[1].plot(t,curve,label=f'Active-layer k = {k:g} W/(m K)',color=c,lw=1.7)
axes[1].set(xlim=(-100,2200),ylim=(1000,1750),xlabel='Time from prompt arrival (s)',ylabel='Front temperature (K)',title='Hottest scheduled C13 material-arrival window')
axes[1].legend(fontsize=8);axes[1].grid(alpha=.15)
labels=['100 nm face each pulse\n2 MJ/kg input','All NRT proxy events\n10 eV/event input','150 µm monthly\n2 MJ/kg input','150 µm each pulse\n10 kJ/kg input','150 µm each pulse\n2 MJ/kg input']
vals=[.052893986668713984,.18940849643963228,.22815423226100842,.3967049000153549,79.34098000307098]
y=np.arange(len(vals));axes[2].barh(y,vals,color=['#4daf4a']*3+['#377eb8','#e41a1c']);axes[2].set_yticks(y,labels,fontsize=8);axes[2].invert_yaxis();axes[2].set(xscale='log',xlim=(.01,250),xlabel='Electrical input (W per receiver m²)',title='Distinct refurbishment operating scenarios')
for yy,v in zip(y,vals):axes[2].text(v*1.08,yy,f'{v:.3g}',va='center',fontsize=8)
axes[2].grid(axis='x',alpha=.15)
fig.savefig(P/'catcher-summary.png',dpi=180)
fig.savefig(P/'catcher-summary.svg')
(P/'plot-provenance.json').write_text(json.dumps({'inputs_sha256':{str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [P/'stopping/surface-stopping.json',P/'geometry/temperature-histories.json',P/'refurbishment/results.json',Path(p['csv'])]},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'caption':'Stopping curves use per-energy mean projected ranges; straggling/reflection/retention are separate. Thermal curves exclude microsecond prompt peaks, which have separately archived upper histories. Refurbishment bars are explicitly different task/cost scenarios, not measured machines or probabilities.'},indent=2)+'\n')
