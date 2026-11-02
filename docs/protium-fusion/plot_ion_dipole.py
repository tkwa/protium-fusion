"""Frozen-orbit spectra and rigid-centroid frequency crossings."""
from pathlib import Path
import argparse,json,gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from electron_wall_orbits import Cylinder


def plot(root,data,out):
    s=json.load(gzip.open(data/'ion-dipole-spectrum-refined.json.gz','rt'))
    other=json.load(gzip.open(data/'ion-dipole-spectrum-refined-independent.json.gz','rt'))
    d=json.loads((data/'rigid-ion-dipole.json').read_text());c=Cylinder(d['source'])
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'protium-ion-dipole-v1'})
    fig,ax=plt.subplots(1,2,figsize=(11,4.7),layout='constrained')
    for i,(name,color) in enumerate([('Protons','#126e82'),('Titanium-49','#b34327')]):
        for j,src in enumerate([s,other]):
            row=src['rows'][i];freq=-(np.array(row['rotating_frequency_rad_s'])+src['ion_rotation_rad_s'])/(2*np.pi*1e6)
            p=np.array(row['bare_rotating_heat_MeV_per_primary_per_E2']);mask=(freq>20)&(freq<255);order=np.argsort(freq[mask])
            ax[0].semilogy(freq[mask][order],p[mask][order],color=color,alpha=1 if j==0 else .45,lw=1.5 if j==0 else .9,label=name if j==0 else None)
    ax[0].set(xlabel='Clockwise laboratory drive frequency (MHz)',ylabel='Rotating heat / primary / field² (MeV per (V/m)²)',title='Frozen-ion response to a uniform circular field',xlim=(20,255),ylim=(1e-8,1e-3))
    ax[0].legend(frameon=False);ax[0].grid(alpha=.2,which='major')
    r=np.linspace(c.R*1e-5,c.R,1001)
    f=-c.T*1e6/c.scale*c.spline(r/c.scale,1)/(c.B*r)/(2*np.pi*1e6)
    ax[1].plot(r*1e3,f,color='#765da6',lw=2,label='Lepton guiding-center rotation')
    for row in d['modes']:
        if not row['lepton_ExB_matching_radii_m']:continue
        frequency=-row['omega_real_rad_s']/(2*np.pi*1e6);radius=row['lepton_ExB_matching_radii_m'][0]*1e3
        ax[1].axhline(frequency,color='#777777',ls='--',lw=1)
        ax[1].plot(radius,frequency,'o',color='#333333',ms=5)
        ax[1].annotate(f'{radius:.3f} mm',xy=(radius,frequency),xytext=(radius-3,frequency-25 if radius<10 else frequency+25),fontsize=9,arrowprops=dict(arrowstyle='-',color='#777777'))
    ax[1].set(xlabel='Radius (mm)',ylabel='Clockwise laboratory frequency (MHz)',title='Two rigid-centroid poles cross the clearing sweep',xlim=(0,c.R*1e3),ylim=(0,275))
    ax[1].legend(frameon=False,fontsize=9,loc='upper right');ax[1].grid(alpha=.2)
    fig.suptitle('Transverse clearing in the 1.6 MeV nonneutral source',fontsize=14)
    fig.supxlabel('Two 16,384-particle orbit ensembles; dashed lines are rigid-centroid model frequencies.\nInduced fields, kinetic damping, lepton capture and complete electrical cost require a coupled model.',fontsize=9)
    fig.savefig(out.with_suffix('.png'),dpi=180);fig.savefig(out.with_suffix('.svg'),metadata={'Date':None});plt.close(fig)
    svg=out.with_suffix('.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--data',type=Path,default=Path(__file__).parent);ap.add_argument('--out',type=Path);a=ap.parse_args();plot(a.root,a.data,a.out or a.data/'ion-dipole')
