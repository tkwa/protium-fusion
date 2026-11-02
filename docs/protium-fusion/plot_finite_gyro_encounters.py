from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args();root=args.outdir
    v=json.loads((root/'finite-gyro-encounter-validation.json').read_text())
    errors=np.array(next(c['metric'] for c in v['checks'] if c['name']=='Impact and gyroangle refinement reduces both errors'))
    a=json.loads((root/'finite-gyro-pilot-analysis.json').read_text());rows=a['source_cases']
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.6))
    for i,label,color in [(0,'Impulse / energy covariance','#3979a8'),(1,'Mean drift','#b64b39')]:
        axes[0].semilogy(range(3),errors[:,i],'-o',label=label,color=color)
    axes[0].set_xticks(range(3),['40 × 8 × 8','64 × 16 × 16','96 × 24 × 24'])
    axes[0].set(xlabel='Impact-radius × impact-angle × gyroangle nodes',ylabel='Relative error against spectral calculation',title='Independent orbit calculation converges')
    axes[0].legend(fontsize=9);axes[0].grid(alpha=.15)
    warm=sorted([r for r in rows if r['vperp_m_s']>0 and r['quadrature']['nk']==96],key=lambda r:r['cut_factor'])
    for key,label,color in [('perpendicular_heating_eV_s','Perpendicular energy','#3979a8'),('parallel_heating_eV_s','Parallel energy','#b64b39'),('total_heating_eV_s','Total energy','#6b6a42')]:
        axes[1].plot([r['cut_factor'] for r in warm],[r['matched'][key]/1e6 for r in warm],'-o',label=label,color=color)
        fine=next(r for r in rows if r['quadrature']['nk']==160)
        axes[1].plot(1,fine['matched'][key]/1e6,'D',color=color,markerfacecolor='white',markersize=7)
    axes[1].set_xscale('log',base=2);axes[1].set_xticks([1,2,4],['1','2','4'])
    axes[1].set(xlabel='Filter wave number / (electron gyrofrequency / Ti thermal speed)',ylabel='Mean energy change (MeV/s)',title='Warm source probe still needs refinement')
    axes[1].legend(loc='center right',fontsize=9);axes[1].axhline(0,color='0.5',lw=.7);axes[1].grid(alpha=.15)
    fig.suptitle('Finite electron gyroradius: independent verification and source pilot',fontsize=15,y=.99)
    fig.text(.5,.03,'Left: smooth screened potential, full Lorentz response, no ion recoil. Right: median Ti radius, prescribed total log 15.\nOpen diamonds use finer quadrature; each point describes one velocity state, not a stationary population.',ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.12,1,.93])
    fig.savefig(root/'finite-gyro-encounter-summary.png',dpi=160);fig.savefig(root/'finite-gyro-encounter-summary.svg')
    svg=root/'finite-gyro-encounter-summary.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')

if __name__=='__main__':main()
