"""Show encounter-scale ordering and the limited cold-pulse response."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cold_magnetized_response import pulse_fractions

ROOT=Path(__file__).parent
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
 scales=json.loads((args.outdir/'magnetized-lepton-scales.json').read_text())['cases']
 rates=json.loads((args.outdir/'cold-magnetized-response.json').read_text())['cases']
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 radii=np.array([r['radius_m'] for r in scales])*1000
 axes[0,0].semilogy(radii,[r['electron_rms_perpendicular_gyroradius_m']/r['local_state']['ion_static_screening_m'] for r in scales],'o-',label='Electron (earlier energy benchmark)',color='#266c8b')
 for name,color in [('H','#ad604b'),('Ti49','#4a8067')]:
  axes[0,0].semilogy(radii,[next(s for s in r['species'] if s['name']==name)['ion_gyroradius_over_static_screening'] for r in scales],'o-',label=name,color=color)
 axes[0,0].axhline(1,color='.5',lw=1,ls=':')
 axes[0,0].set(xlabel='Source radius (mm)',ylabel='Rms gyroradius / ion static screening length',title='The three species have different encounter regimes')
 axes[0,0].legend(fontsize=8)
 aa=np.geomspace(.001,20,400);par,perp=pulse_fractions(aa,1/3)
 axes[0,1].semilogx(aa,par,label='Parallel energy',color='#266c8b')
 axes[0,1].semilogx(aa,perp,label='Transverse energy',color='#ad604b')
 axes[0,1].semilogx(aa,par+perp,label='Total energy',color='#4a8067')
 axes[0,1].set(xlabel='Gyrofrequency × encounter duration  |ωₑ|b/u',ylabel='Fraction of zero-field pulse energy',ylim=(0,1.08),title='Cold electron, isotropically oriented weak pulses')
 axes[0,1].legend(fontsize=8)
 for name,color in [('H','#ad604b'),('Ti49','#266c8b')]:
  chosen=[r for r in rates if r['species']==name and r['outer_fraction_of_static_screening']==.1 and r['inner_scale_multiplier']==10.]
  axes[1,0].semilogy(radii,[r['total_heating_eV_s']/1e6 for r in chosen],'o-',color=color,label=name+', with magnetic field')
  axes[1,0].semilogy(radii,[r['unmagnetized_heating_eV_s']/1e6 for r in chosen],'--',color=color,label=name+', zero field')
 axes[1,0].set(xlabel='Source radius (mm)',ylabel='Initial electron heating (MeV/s)',title='Cutoff screen: inner scale ×10, outer scale ×0.1')
 axes[1,0].legend(fontsize=8)
 for inner,color in [(1.,'#266c8b'),(10.,'#ad604b')]:
  for outer,style in [(.1,'-'),(1.,'--')]:
   chosen=[r for r in rates if r['species']=='Ti49' and r['outer_fraction_of_static_screening']==outer and r['inner_scale_multiplier']==inner]
   axes[1,1].plot(radii,[r['magnetized_over_unmagnetized'] for r in chosen],'o',ls=style,color=color,label=f'Inner ×{inner:g}, outer ×{outer:g}')
 axes[1,1].set(xlabel='Source radius (mm)',ylabel='Initial heating / same-cutoff zero-field result',ylim=(.4,.9),title='Cutoff choices retain substantial axial heating')
 axes[1,1].legend(fontsize=8)
 for ax in axes.flat:ax.grid(alpha=.2)
 fig.suptitle('Magnetized encounter scales and a cold-electron response test',fontsize=15)
 fig.supxlabel('Prescribed straight-ion force pulses. The cutoff rate screen does not determine the warm electron distribution, residence or plant cost.',fontsize=9)
 for ext in ['png','svg']:fig.savefig(args.outdir/f'cold-magnetized-summary.{ext}',dpi=170)
 path=args.outdir/'cold-magnetized-summary.svg';path.write_text('\n'.join(s.rstrip() for s in path.read_text().splitlines())+'\n')
