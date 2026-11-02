"""Plot the finite-axial-energy distant-encounter diagnostic."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args()
    data=json.loads((args.outdir/'guiding-center-response.json').read_text())
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(12.6,5.4))
    colors=['#4a759c','#ad4937']
    for case,color in zip(data['cases'],colors):
        rows=case['states'][1:];E=np.array([r['axial_energy_eV'] for r in rows])
        label='Ion relative flow included' if case['include_relative_flow'] else 'No relative flow'
        ax[0].plot(E,[r['axial_heating_eV_s']/1e6 for r in rows],color=color,label=label)
        ax[1].plot(E,[r['electron_radial_drift_m_s'] for r in rows],color=color,label=label)
    ax[0].set(yscale='symlog',ylim=(-250,2500),ylabel='Mean axial energy change (MeV/s)')
    ax[1].set_yscale('symlog',linthresh=.001)
    ax[1].set(ylabel='Electron guiding-center radial drift (m/s)')
    for a in ax:
        a.set_xscale('log');a.set_xlabel('Initial axial kinetic energy (eV)')
        a.axhline(0,color='0.5',lw=.7);a.grid(alpha=.15)
    ax[0].legend(fontsize=9,loc='upper right')
    fig.suptitle('Distant encounters: consistent drift changes the heating result',fontsize=16,y=.99)
    fig.text(.5,.025,'Median titanium radius · 137.15 T · prescribed log interval 4 · fixed magnetic moment\nLocal weak-encounter limit; this is not a stationary electron population or a revised plant load.',ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.11,1,.94])
    fig.savefig(args.outdir/'guiding-center-summary.png',dpi=160)
    fig.savefig(args.outdir/'guiding-center-summary.svg')

if __name__=='__main__':main()
