"""Scientific diagnostic figures for the conditional lepton-birth ledger."""
from pathlib import Path
import argparse,gzip,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot(data,out):
    c=json.loads(gzip.decompress((data/'bricc-pair-coefficients-refined.json.gz').read_bytes()))
    p=json.loads((data/'prompt-pair-supply-refined.json').read_text())
    d=json.loads((data/'residual-lepton-births.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans','svg.hashsalt':'protium-prompt-pairs-v1'})
    fig,ax=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    colors=['#126e82','#b34327','#298c55','#a57414','#765da6','#797979']
    for i,(m,color) in enumerate(zip(c['mult'],colors)):
        ax[0].semilogy(c['energy_MeV'],np.asarray(c['alpha_pair_over_gamma']['23'])[:,i],label=m,color=color,lw=1.7)
    ax[0].set(xlabel='Transition energy (MeV)',ylabel='Pair / gamma partial-rate ratio',title='Evaluated ordinary pair conversion · Z = 23',xlim=(1.02,8))
    ax[0].legend(ncol=3,frameon=False,loc='lower right');ax[0].grid(alpha=.18,which='major')
    for pop,name,style in [('ground_state','Ground-state neutral-rate proxy','-'),('fastest_known_weak_state','Fastest known weak-state proxy','--')]:
        s=[s for s in d['scenarios'] if s['population']==pop]
        ax[1].loglog([x['residence_s'] for x in s],[x['total_neutral_rate_proxy_per_primary'] for x in s],style,color='#765da6',lw=1.8,label=name)
    for m,color in [('E1','#126e82'),('M3','#797979')]:
        ax[1].axhline(p['per_primary']['pair_branch_'+m],color=color,lw=1.5,label=f'Photon cascade assigned {m}')
    ax[1].set(xlabel='Assigned residual residence (s)',ylabel='Birth / decay proxy per primary Ti(p,n)',title='Finite-residence decay versus prompt pairs',xlim=(1e-6,100),ylim=(1e-12,.01))
    ax[1].grid(alpha=.18,which='major');ax[1].legend(frameon=False,fontsize=8,loc='lower right')
    fig.suptitle('Conditional lepton production in the 1.6 MeV catalyst network',fontsize=14)
    fig.supxlabel('Pair curves cover 1.1–8 MeV; E0 and stripped-ion cascade changes are unresolved.\nDecay curves deliberately count EC as a lepton proxy and prescribe the residence.',fontsize=9)
    fig.savefig(out.with_suffix('.png'),dpi=180);fig.savefig(out.with_suffix('.svg'),metadata={'Date':None});plt.close(fig)
    svg=out.with_suffix('.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path(__file__).parent);ap.add_argument('--out',type=Path);args=ap.parse_args();plot(args.data,args.out or args.data/'prompt-pair-supply')
