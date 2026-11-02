"""Plot ideal wall insulation and the self-consistent static axial response."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from electron_wall_orbits import Cylinder
from axial_bias_screen import response
ROOT=Path(__file__).parent

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
 plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
 wall=json.loads((args.outdir/'electron-wall-orbits.json').read_text());c=Cylinder(plant['selected_source_profiles']['22049']);eq=c.eq
 x=np.linspace(eq['x'][0],c.src['xwall'],1001);r=x*c.scale
 h=np.interp(x,eq['x'],eq['h']);n=np.interp(x,eq['x'],eq['c']);birth=h*n*x;birth/=np.trapezoid(birth,r*1000)
 resp=response(c);transfer=np.array([resp['field_transfer'](rr) for rr in r])
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(12,8.4),layout='constrained')
 ax=axes[0,0];ax.plot(r*1000,birth,color='#266c8b',lw=2,label='Radial reaction probability')
 depth=next(q['maximum_access_depth_m'] for q in wall['all_direction_access_bounds'] if q['charge_sign']==-1 and q['initial_MeV']==30)
 ax.axvspan((c.R-depth)*1000,c.R*1000,color='#bb6557',alpha=.3,label='Accessible from wall at 30 MeV')
 ax.axvline(wall['profiles'][0]['reaction_quantiles'][-2]['radius_m']*1000,color='#666',ls=':',label='99% of reactions inside')
 ax.set(xlabel='Radius (mm)',ylabel='Probability density per mm',title='The hot source sits well inside the channel wall',xlim=(0,c.R*1000))
 ax.legend(fontsize=8)
 ax=axes[0,1];rr=np.linspace(.0005,c.R-.000001,160)
 for sign,label,color in [(-1,'Electron','#266c8b'),(1,'Positron','#a55446')]:
  ax.plot(rr*1000,[c.minimum_wall_energy(v,sign) for v in rr],label=label,color=color)
 ax.set(xlabel='Radius to be reached (mm)',ylabel='Necessary initial wall energy (MeV)',title='All-direction access bound from conserved momentum',ylim=(0,450),xlim=(0,c.R*1000))
 ax.legend();ax.grid(alpha=.2)
 ax=axes[1,0];ax.semilogy(r*1000,transfer,color='#447f69',lw=2)
 ax.axvline(wall['profiles'][0]['reaction_quantiles'][1]['radius_m']*1000,color='#666',ls=':',label='Median reaction radius')
 ax.set(xlabel='Radius (mm)',ylabel='Internal potential change / wall voltage change',title='Slow axial wall bias is screened in the reacting core',xlim=(0,c.R*1000),ylim=(1e-5,2))
 ax.legend(fontsize=8);ax.grid(alpha=.2)
 ax=axes[1,1];orbits=wall['titanium_wall_orbits'];energies=sorted(set(o['initial_MeV'] for o in orbits))
 for sign,label,color in [(-1,'Electron','#266c8b'),(1,'Positron','#a55446')]:
  vals=[max(o['return_time_s'] for o in orbits if o['charge_sign']==sign and o['initial_MeV']==K)*1e12 for K in energies]
  ax.loglog(np.array(energies)*1e6,vals,'o-',label=label,color=color)
 ax.set(xlabel='Initial kinetic energy at wall (eV)',ylabel='Largest sampled return time (ps)',title='Conservative wall orbits return on picosecond scales')
 ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle('Electron control in the conditional titanium reactor',fontsize=15)
 fig.supxlabel('Static axisymmetric fields; collisions, radiation reaction, ends and driven collective response require separate models.',fontsize=10)
 for ext in ['png','svg']:fig.savefig(args.outdir/f'electron-control-summary.{ext}',dpi=170)
 svg=args.outdir/'electron-control-summary.svg'
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
