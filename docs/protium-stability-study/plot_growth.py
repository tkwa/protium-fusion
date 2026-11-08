"""Finite-time growth of controlled low-order disturbances during burning."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent;base=Path('/var/tmp/tkwa-chatgpt/protium-stability/runs')
def read(name):
    p=base/name/'diagnostics.json'
    if not p.exists():p=root/'results/cases'/(name+'.json')
    return json.loads(p.read_text())
c=read('A-384x64-smooth');h=c['history'];t=np.array(h['time'])*30
q=np.array(h['product_mass'])/h['product_mass'][-1]
t10,t90=np.interp([.1,.9],q,t)
fig,ax=plt.subplots(1,2,figsize=(9,3.8),layout='constrained')
for name in ['A-384x64-l2-a0001','A-384x64-l2-a001','A-384x64-l2-a003','A-384x64-l2-a01']:
    d=read(name);a=d['config']['amplitude'];rows=d['snapshots']
    tt=[r['time_us'] for r in rows]
    ax[0].plot(tt,[r['density_rms']/a for r in rows],label=f'{100*a:g}% velocity seed')
    ax[1].plot(tt,[r['radii']['0.5']['rms_fraction']/a for r in rows])
for panel in ax:
    panel.axvspan(t10,t90,color='grey',alpha=.15)
    panel.set(xlim=(0,42),xlabel='Time from prepared initial state (microseconds)')
    panel.grid(alpha=.2);panel.spines[['top','right']].set_visible(False)
ax[0].set(ylabel='Angular density RMS / seed amplitude',title='Relative density structure')
ax[1].set(ylabel='RMS distortion / seed amplitude',title='Half-mass surface distortion')
ax[0].legend(frameon=False,fontsize=8)
fig.suptitle('Finite-time amplification; shaded interval produces 10–90% of smooth output')
out=root/'figures';out.mkdir(exist_ok=True)
fig.savefig(out/'disturbance-growth.png',dpi=180);fig.savefig(out/'disturbance-growth.svg')
