#!/usr/bin/env python3
"""Static research figure; bars are conditional accounts, without probability bands."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent

def main():
    d=json.loads((P/'plant-results.json').read_text());cases={x['name']:x for x in d['cases']}
    names=['Selected reference','Large N15 thermal','Thermal, services 700','Hybrid, 1 kW/kg + 0.1 MJ/kg','Hybrid, 10 kW/kg + 1 MJ/kg']
    vals=[d['reference'],cases['spherical_large_N15_16']['net_W_kg'],cases['spherical_large_N15_16_services700']['net_W_kg'],cases['hybrid16']['net_W_kg'],cases['hybrid16_both_strong_ratings']['net_W_kg']]
    fig,ax=plt.subplots(1,2,figsize=(13,4.8),gridspec_kw={'width_ratios':[1.7,1]},layout='constrained')
    colors=['#738095','#247f8a','#60a6aa','#ba815b','#8b6aa4']
    ax[0].barh(names[::-1],vals[::-1],color=colors[::-1]);ax[0].axvline(.0208*16,color='#b95050',ls='--',lw=1.5,label='16 × 20.8 mW/kg')
    for y,v in enumerate(vals[::-1]):ax[0].text(v+.005,y,f'{v:.3f}',va='center',fontsize=9)
    ax[0].set(xlabel='Operating net electricity / installed mass (W/kg)',xlim=(0,.415),title='Conditional plant accounts')
    ax[0].legend(loc='lower right',fontsize=8);ax[0].grid(axis='x',alpha=.2);ax[0].set_axisbelow(True)
    z=[32,64,128,256];q=[]
    for n in z:q.append(json.loads((P/f'pulse-n15-m100-v1-z{n}-t1.json').read_text())['desired_product_per_initial_catalyst'])
    ax[1].plot(z,q,'o-',color='#247f8a',label='100× N15 pulse mass')
    small=json.loads((P/'pulse-n15-m1-v1-z128-t1.json').read_text())['desired_product_per_initial_catalyst']
    ax[1].scatter([128],[small],marker='s',color='#738095',label='Reproduced small control',zorder=5)
    ax[1].set(xlabel='Radial cells',ylabel='N15 converted / initial N15',title='New spherical pulse results',ylim=(.64,.83),xticks=z)
    ax[1].grid(alpha=.2);ax[1].legend(fontsize=8)
    fig.suptitle('Larger N15 shots relieve dispatch; receiver and source qualification remain open',fontsize=12)
    fig.savefig(P/'comparison.png',dpi=170);fig.savefig(P/'comparison.svg')

if __name__=='__main__':main()
