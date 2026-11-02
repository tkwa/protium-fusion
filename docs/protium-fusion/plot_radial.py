"""Create a scientific figure from saved spatial-model outputs."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).parent
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11.7,4.8))
colors=['#176b93','#c76b1b','#6654a3']
names=['radial-A-100m-64','radial-A-300m-1024-reduced','radial-B2-300m-64-rad']
labels=['N14 stage · 100 m','N14 stage · 300 m','C12 stage · 300 m']
for name,label,col in zip(names,labels,colors):
    r=json.loads((ROOT/(name+'.json')).read_text())
    mass=np.array(r['profile']['mass_fraction_of_total'])
    m=np.cumsum(mass)-mass/2
    axes[1].plot(m,100*np.array(r['profile']['desired_product_per_initial_catalyst']),
                 color=col,label=label,lw=2)
axes[1].set(xlabel='Enclosed initial mass / total mass',ylabel='Desired isotope conversion (%)',
            ylim=(-2,104),xlim=(0,1),title='Surface expansion leaves an unburned layer')
axes[1].legend(frameon=False,fontsize=9,loc='lower left')

data=np.load(ROOT/'radial-A-300m-1024-reduced.npz')
n=int(data['zones']);nv=int(data['variables_per_zone']);nn=len(data['species'])
t=data['time_s'];z=data['state'].reshape((n,nv,-1))
x=np.linspace(0,1,n+1);mass=np.diff(x**3);m=np.cumsum(mass)-mass/2
for target,col in zip([0,20,50,100,340],['#999999','#e1a44c','#cb6337','#377ba8','#142f4a']):
    k=int(np.argmin(abs(t*1e6-target)))
    axes[0].plot(m,z[:,nn,k]*300,color=col,label=f'{t[k]*1e6:.0f} µs',lw=1.8)
axes[0].set(xlabel='Enclosed initial mass / total mass',ylabel='Temperature (million K)',
            xlim=(0,1),ylim=(0,1550),title='Nitrogen stage · initial radius 300 m')
axes[0].legend(frameon=False,fontsize=9,ncol=3)
for ax in axes: ax.grid(axis='y',alpha=.18)
fig.suptitle('Opaque CNO pulses: first spatial checks',fontsize=17,x=.06,ha='left',y=.98)
fig.text(.06,.04,'Initially compressed fuel: ρ = 10⁸ kg/m³, T = 300 million K. Spherical model; no implosion or driver.\n'
         'Radiation runs use constant opacity 0.001 m²/kg. Screening and Coulomb EOS corrections remain to be tested.',
         fontsize=9,color='#555555')
fig.subplots_adjust(left=.07,right=.98,top=.80,bottom=.23,wspace=.28)
fig.savefig(ROOT/'radial-pulse-check.png',dpi=170)
fig.savefig(ROOT/'radial-pulse-check.svg')
svg=ROOT/'radial-pulse-check.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
print(ROOT/'radial-pulse-check.png')
