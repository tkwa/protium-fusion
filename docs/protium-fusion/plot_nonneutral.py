"""Plot conditional conversion efficiencies, without implying plant feasibility."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.colors import TwoSlopeNorm
from nonneutral_screen import nuclear_and_drag,angular_recycle
ROOT=Path(__file__).parent
plt.rcParams['svg.hashsalt']='protium-nonneutral'
d=nuclear_and_drag(4.5,.1)
a=np.linspace(.8,1,241);c=np.linspace(.8,1,241);A,C=np.meshgrid(a,c)
th=.4;g=C+(1-C)*th
base=g*d['directly_recoverable_ion_energy_MeV_per_pn']+th*d['thermal_energy_MeV_per_pn']-d['beam_preparation_MeV_per_pn']/A
side=50*(d['B10_residual_per_pn']+d['N14_capture_per_pn'])
base-=side
extra=angular_recycle(d,np.pi/2)['additional_recirculated_beam_energy_MeV_per_pn']
fig,axes=plt.subplots(1,2,figsize=(11.5,6.2),sharex=True,sharey=True)
fig.subplots_adjust(left=.08,right=.88,bottom=.2,top=.88,wspace=.08)
for ax,z,title in zip(axes,[base,base-extra*(1/A-g)],['Energy-transfer screen','With 90° acceptance and ion recycling']):
 im=ax.pcolormesh(A,C,z,cmap='RdBu',norm=TwoSlopeNorm(vmin=-10,vcenter=0,vmax=5),shading='auto',rasterized=True)
 cs=ax.contour(A,C,z,levels=[0,1,3],colors=['black','#444444','#666666'],linewidths=[1.5,.8,.8])
 ax.clabel(cs,fmt=lambda x:f'{x:g} MeV',fontsize=9)
 ax.scatter([.95],[.95],s=35,facecolor='white',edgecolor='black',zorder=4)
 ax.set_title(title,fontsize=12)
 ax.set_xlabel('Electricity-to-ion acceleration efficiency',labelpad=9)
 ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0))
 ax.set_xticks([.8,.85,.9,.95,1.]);ax.set_yticks([.8,.85,.9,.95,1.])
 ax.set_aspect('equal')
axes[0].set_ylabel('Ion energy converted directly to electricity')
cb=fig.colorbar(im,cax=fig.add_axes([.91,.23,.018,.6]),extend='min');cb.set_label('Net MeV per useful proton-to-neutron event')
fig.suptitle('Carbon-13 beam cycle: recovery requirements before remaining plant costs',fontsize=13,y=.98)
fig.text(.08,.055,'4.5 MeV protons; 0.1 MeV RMS width; carbon Coulomb log 25; 40% heat conversion.\nClosed deuterium cycle assumed. White dot: 95% / 95%. Neither panel establishes a reactor design.',fontsize=9)
fig.savefig(ROOT/'nonneutral-recovery.png',dpi=160,bbox_inches='tight')
fig.savefig(ROOT/'nonneutral-recovery.svg',dpi=160,bbox_inches='tight',metadata={'Date':None})
p=ROOT/'nonneutral-recovery.svg'
p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
