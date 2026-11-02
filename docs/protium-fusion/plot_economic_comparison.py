"""Scenario comparisons for plant energy capital and local solar collection."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(data,out):
    d=json.loads((data/'economic-comparison.json').read_text());s=d['nonneutral_reference']['net_W_kg']
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'protium-economics-v1'})
    fig,ax=plt.subplots(1,3,figsize=(13,4.7),layout='constrained')
    e=np.geomspace(.1,1000,200)
    for spec,name,color in [(.001,'1 mW/kg target','#765da6'),(s,'Conditional ion plant','#126e82')]:
        ax[0].loglog(e,e*1e6/(spec*d['constants']['seconds_per_year']),color=color,lw=2,label=name)
    ax[0].set(xlabel='Construction energy (MJ/kg)',ylabel='Simple energy payback (years)',title='Slow repayment at low specific power',xlim=(.1,1000))
    ax[0].legend(frameon=False,fontsize=8);ax[0].grid(alpha=.2,which='major')
    r=np.geomspace(10,10000,200)
    for sigma,color in [(.1,'#298c55'),(1.,'#126e82'),(10.,'#b34327')]:
        ax[1].loglog(r,.3*1361/(sigma*r*r)*1e3,lw=2,color=color,label=f'{sigma:g} kg/m² complete collector')
    ax[1].axhline(1,color='#765da6',ls='--',label='1 mW/kg')
    ax[1].set(xlabel='Distance from a Sun-like star (AU)',ylabel='Net solar specific power (mW/kg)',title='Solar comparison · 30% net conversion',xlim=(10,10000),ylim=(.001,1e5))
    ax[1].legend(frameon=False,fontsize=8,loc='upper right');ax[1].grid(alpha=.2,which='major')
    rows=d['opaque_replacement_cases'];x=np.array([v['replacement_fraction_year'] for v in rows])*100
    ax[2].plot(x,[v['throughput_over_original'] for v in rows],'-o',lw=2,color='#b34327',ms=5)
    for i,row in enumerate(rows):
        if row['radiator_area_exceeds_collector_area']:
            ax[2].plot(x[i],row['throughput_over_original'],'o',mfc='white',mec='#b34327',ms=7)
            ax[2].annotate('Radiator exceeds\ncollector surface',xy=(x[i],row['throughput_over_original']),xytext=(.047,3.9),fontsize=8,arrowprops={'arrowstyle':'-','color':'#777777'})
    ax[2].set(xlabel='Installed mass replaced per year (%)',ylabel='Pulse throughput / original allocation',title='Opaque plant · retain 1 mW/kg net',xlim=(0,.1),ylim=(.9,4.5))
    ax[2].grid(alpha=.2)
    fig.suptitle('Economic requirements of low-specific-power hydrogen fusion',fontsize=14)
    fig.supxlabel('Scenario inputs, not demonstrated component performance. Opaque replacement uses 100 MJ/kg manufacturing,\n1 kg/W factory capacity and ten years of spares; higher pulse cadence still requires engineering qualification.',fontsize=9)
    fig.savefig(out.with_suffix('.png'),dpi=180);fig.savefig(out.with_suffix('.svg'),metadata={'Date':None});plt.close(fig)
    svg=out.with_suffix('.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path(__file__).parent);ap.add_argument('--out',type=Path);a=ap.parse_args();plot(a.data,a.out or a.data/'economic-comparison')
