"""Plot the diagnosed cold-tail failure and completed grid refinements."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser()
ap.add_argument('--raw',type=Path,required=True)
ap.add_argument('--old-trace',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args()
fig,ax=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
for path,label,color in [(a.old_trace,'Original viscosity: failed','#a33b33'),
 (a.raw/'B2-64-compression-cap.npz','Corrected viscosity: completed','#167a69')]:
 with np.load(path) as data:
  nn=len(data['species']);nv=int(data['variables_per_zone']);nz=int(data['zones'])
  temperatures=data['state'].reshape(nz,nv,-1)[:,nn,:]*300
  ax[0].semilogy(data['time_s']*1e6,temperatures[-1],label=label,color=color,lw=2)
ax[0].set(xlabel='Time (microseconds)',ylabel='Outer-cell temperature (MK)',
 title='The artificial cold tail disappears')
ax[0].legend(frameon=False,fontsize=9,loc='lower left')
zones=[];gains=[]
for zone in [64,128,256]:
 name=f'B2-{zone}-compression-cap'+('' if zone==64 else '-energy-matched')
 path=a.raw/(name+'.json')
 if path.exists():
  r=json.loads(path.read_text());zones.append(zone);gains.append(r['deposited_nuclear_over_initial_energy'])
ax[1].plot(zones,gains,'o-',color='#375a8c',lw=2)
for zone,gain in zip(zones,gains):ax[1].annotate(f'{gain:.4f}',(zone,gain),xytext=(0,9),textcoords='offset points',ha='center')
ax[1].set(xlabel='Number of radial cells',ylabel='Prompt deposited nuclear gain',
 title='Identical preparation, matched input energy')
ax[1].set_xticks(zones)
ax[1].margins(x=.15,y=.25)
for axis in ax:
 axis.grid(alpha=.2)
 axis.spines[['top','right']].set_visible(False)
fig.savefig(a.output.with_suffix('.png'),dpi=180)
fig.savefig(a.output.with_suffix('.svg'))
