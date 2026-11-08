"""Plot subjective threshold intervals and their dependence on assembly judgment."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
d=json.loads((root/'results/subjective-assessment.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')

def interval(ax,y,s,label,color):
    lo,hi=s['interval80'];m=s['median']
    ax.errorbar(m,y,xerr=[[m-lo],[hi-m]],fmt='o',color=color,capsize=5,lw=3)
    ax.text(m,y+.14,label,ha='center',va='bottom',fontsize=9)

a=axes[0]
interval(a,0,d['arrival_time_us'],'Subjective median and 80% interval','#2874a6')
a.axvline(d['anchors']['timing_us'],ls='--',color='#555555',label='Hot-flow numerical anchor')
a.set(xscale='log',xlabel='Maximum RMS arrival error (microseconds)',yticks=[],ylim=(-.5,.65),title='Timing tolerance')
a.legend(loc='lower left',frameon=False,fontsize=8)

a=axes[1]
for y,(key,s) in enumerate(d['packet_count_prior_sensitivity'].items()):
    interval(a,y,s,f'Assumed median compression factor {float(key):g}', '#2874a6' if float(key)==4 else '#888888')
a.axvline(d['anchors']['hot_packet_count'],ls='--',color='#555555',label='Hot-flow numerical anchor')
a.set(xscale='log',xlabel='Minimum shaped-packet count',yticks=[],ylim=(-.5,2.65),title='Packet count: assembly-prior sensitivity')
a.legend(loc='lower left',frameon=False,fontsize=8)
for a in axes:
    a.grid(axis='x',alpha=.2)
    for side in ['top','right','left']:a.spines[side].set_visible(False)
fig.suptitle('At least 90% of smooth-pulse output; isolated limits, conditional on feasible preparation',fontsize=11)
fig.savefig(root/'figures/subjective-intervals.png',dpi=180)
fig.savefig(root/'figures/subjective-intervals.svg')
