"""Plot the prescribed gyro-averaged electron escape benchmark."""
from pathlib import Path
import argparse,json,gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).parent
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    main=json.loads((args.outdir/'lepton-escape-main.json').read_text())
    analysis=json.loads((args.outdir/'lepton-escape-analysis.json').read_text())
    prior=json.loads(gzip.decompress((args.root/'lepton-collision-screen.json.gz').read_bytes()))
    iso=next(r for r in prior['cases'] if r['reaction_fraction_inside']==.5 and r['log_lambda']==15 and r['include_relative_flow'])
    trajectories=json.loads(gzip.decompress((args.outdir/'lepton-escape-main-trajectories.json.gz').read_bytes()))
    times=np.sort(np.array(trajectories['rows'])[:,0])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    edges=np.array(main['energy_histogram_edges_eV'])/1000
    density=np.array(main['occupation_probability_per_energy_bin'])/np.diff(edges)
    axes[0,0].stairs(density,edges,color='#266c8b',label='Gyro-averaged escape trajectories')
    axes[0,0].plot(np.array(prior['energy_eV'])/1000,np.array(iso['density_per_eV'])*1000,
        color='#ad604b',ls='--',label='Earlier isotropic energy closure')
    axes[0,0].set(xlim=(0,12),xlabel='Energy in local E × B frame (keV)',ylabel='Occupation density per keV',title='The energy distribution changes little')
    axes[0,0].legend(fontsize=8)
    pitch=np.array(main['occupation_probability_per_pitch_bin'])/.05
    axes[0,1].stairs(pitch,main['pitch_bin_edges'],color='#266c8b')
    axes[0,1].axhline(.5,color='#ad604b',ls='--',label='Isotropic distribution')
    axes[0,1].set(xlabel='Parallel velocity / total speed',ylabel='Occupation density',ylim=(.35,.7),title='Parallel and transverse energies differ')
    axes[0,1].legend(fontsize=8)
    axes[1,0].plot(times*1000,1-(np.arange(len(times))+.5)/len(times),color='#266c8b',label='4,096 electron trajectories')
    t=np.linspace(0,.04,401);odd=np.arange(1,602,2)
    L=main['active_length_m']
    for D,label,color in [(iso['isotropized_green_kubo_diffusion_m2_s'],'Diffusion: earlier isotropic closure','#ad604b'),
        (analysis['finite_length']['free_space_pooled_diffusion_m2_s']['mean'],'Diffusion: free-space control','#4a8067')]:
        survival=8/np.pi**2*np.sum(np.exp(-np.pi**2*D*t[:,None]*odd[None,:]**2/L**2)/odd[None,:]**2,axis=1)
        axes[1,0].plot(t*1000,survival,ls='--',color=color,label=label)
    axes[1,0].set(xlim=(0,30),ylim=(0,1),xlabel='Time since birth (ms)',ylabel='Fraction remaining',title='Absorbing ends require a kinetic correction')
    axes[1,0].legend(fontsize=8)
    lengths=[];means=[];errors=[]
    for name in ['half','main','double']:
        row=json.loads((args.outdir/f'lepton-escape-{name}.json').read_text())
        lengths.append(row['active_length_m']/1000)
        means.append(row['residence_s']['mean']*1000)
        errors.append(row['residence_s']['standard_error']*1000)
    xx=np.linspace(0,2.05,200);a,b,c=[analysis['finite_length']['coefficients_s'][k] for k in ['a','b','c']]
    axes[1,1].errorbar(lengths,means,yerr=errors,fmt='o',color='#266c8b',capsize=3,label='Trajectory mean ± sampling SE')
    axes[1,1].plot(xx*L/1000,1000*(a*xx*xx+b*xx+c),color='#266c8b',label='Three-point quadratic interpolation')
    D=analysis['finite_length']['free_space_pooled_diffusion_m2_s']['mean']
    axes[1,1].plot(xx*L/1000,1000*(xx*L)**2/(12*D),color='#4a8067',ls='--',label='Bulk diffusion term alone')
    axes[1,1].set(xlim=(0,42),ylim=(0,18),xlabel='Absorbing interval length (km)',ylabel='Mean residence (ms)',title='Length dependence agrees with the free-space control')
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Electron escape at the median titanium reaction radius',fontsize=15)
    fig.supxlabel('Prescribed unmagnetized encounter coefficients, fixed radius, 1 eV births and ideal absorbing ends. Strong-field kinetics remain open.',fontsize=9)
    for ext in ['png','svg']:fig.savefig(args.outdir/f'lepton-escape-summary.{ext}',dpi=170)
    path=args.outdir/'lepton-escape-summary.svg'
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
