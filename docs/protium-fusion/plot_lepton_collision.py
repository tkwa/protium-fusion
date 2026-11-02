"""Plot the diagnostic collision/radiation closure with its scope visible."""
from pathlib import Path
import argparse,json,gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).parent
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=ROOT)
    ap.add_argument('--root',type=Path,default=ROOT);args=ap.parse_args()
    data=json.loads(gzip.decompress((args.outdir/'lepton-collision-screen.json.gz').read_bytes()))
    E=np.array(data['energy_eV']);rows=data['cases']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for logL,color in [(5.,'#ad604b'),(15.,'#266c8b'),(30.,'#4a8067')]:
        selected=[r for r in rows if r['log_lambda']==logL and r['include_relative_flow']]
        mid=next(r for r in selected if r['reaction_fraction_inside']==.5)
        radii=np.array([r['local_state']['radius_m'] for r in selected])*1000
        label=f'Prescribed ln Λ = {logL:g}'
        axes[0,0].plot(E/1000,np.array(mid['density_per_eV'])*1000,label=label,color=color)
        axes[0,1].plot(radii,[r['radiation_eV_s']/1e6 for r in selected],'o-',label=label,color=color)
        axes[1,0].plot(radii,[r['electron_radial_force_drift_m_s']*1000 for r in selected],'o-',color=color)
        axes[1,1].plot(radii,[r['isotropized_green_kubo_diffusion_m2_s']/1e9 for r in selected],'o-',color=color)
    axes[0,0].set(xlabel='Kinetic energy in local E × B frame (keV)',ylabel='Probability density per keV',xlim=(0,16),title='Energy distribution at the median reaction radius')
    axes[0,0].legend(fontsize=8)
    axes[0,1].set(xlabel='Radius (mm)',ylabel='Radiation per lepton (MeV/s)',title='Ion collisions replenish the radiated energy')
    axes[0,1].legend(fontsize=8)
    axes[1,0].set(xlabel='Radius (mm)',ylabel='Electron outward force drift (mm/s)',title='Local angular-momentum transfer drives radial motion')
    axes[1,1].set(xlabel='Radius (mm)',ylabel='Homogeneous axial diffusivity (10⁹ m²/s)',title='Energy-dependent velocity correlations affect transport')
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Trace-lepton kinetic benchmark for the titanium source',fontsize=15)
    fig.supxlabel('Diagnostic unmagnetized, isotropic closure. Strong-field collisions, birth rates and finite collectors remain unresolved.',fontsize=10)
    for ext in ['png','svg']:fig.savefig(args.outdir/f'lepton-collision-summary.{ext}',dpi=170)
    svg=args.outdir/'lepton-collision-summary.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
