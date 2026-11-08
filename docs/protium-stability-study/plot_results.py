"""Render scientific figures from retained diagnostics and selected field dumps."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from analyze import snapshot
from packet_geometry import meridional

root=Path(__file__).resolve().parent
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
out=root/'figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})

def read(name):
    p=base/'runs'/name/'diagnostics.json'
    if not p.exists():p=root/'results/cases'/(name+'.json')
    if not p.exists():return None
    x=json.loads(p.read_text());e=x['execution']
    return x if e and e.get('reached_time_limit') and not e.get('fatal_error') else None

def output(d):
    return d.get('nuclear_output_J',d['gain']*d['initial_energy_J'])

fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
control=read('A-384x64-ref-cos2-s0us-hlle')
if control:
    for ell,color in [(1,'#2874a6'),(2,'#d35400'),(8,'#239b56'),(16,'#7d3c98')]:
        rows=[(0,100)]
        for case in (base/'runs').glob('A-384*ref-cos*-s*us*'):
            p=case/'diagnostics.json'
            if not p.exists():continue
            d=read(case.name)
            if not d or d['config']['ell']!=ell or not d['config'].get('binary') or 'hlle' not in d['config']['binary']:continue
            if d['config']['nt']!=(128 if ell==16 else 64):continue
            if d['config'].get('phase_taper_m',0)>0:continue
            sigma=d['config'].get('timing_us',0)
            if sigma:rows.append((sigma,100*output(d)/output(control)))
        rows=sorted(set(rows))
        if len(rows)>1:
            x,y=np.array(rows).T;ax[0].plot(x,y,'o-',color=color,label=f'Angular order {ell}')
early=read('A-384x64-ref0-cos8-s0us')
if early:
    rows=[(0,100)]
    for sigma in [3,4,4.5,5,5.5,6]:
        d=read(f'A-384x64-ref0-cos8-s{sigma}us')
        if d:rows.append((sigma,100*output(d)/output(early)))
    if len(rows)>1:
        x,y=np.array(rows).T;ax[0].plot(x,y,'^--',color='#17202a',label='Order 8, earlier seed')
for grid,style in [(32,'s--'),(48,'s-'),(64,'D-')]:
    c=read(f'A-cart{grid}-wide-ref0-gaussian-s0us')
    if not c:continue
    rows=[(0,100)]
    for s in [1,2,3,4,4.5,5,6]:
        d=read(f'A-cart{grid}-wide-ref0-gaussian-s{s}us')
        if d:rows.append((s,100*output(d)/output(c)))
    if len(rows)>1:
        x,y=np.array(rows).T;ax[0].plot(x,y,style,label=f'3D patch jitter, base {grid}',color='#b03a2e',alpha=.7 if grid==32 else 1)
ax[0].set(xlabel='RMS arrival-time error (microseconds)',ylabel='Nuclear output retained (%)',title='Timing errors in the evolved pulse')
if len(ax[0].lines):ax[0].legend(frameon=False,fontsize=8)
smooth=read('A-384x128-smooth-hllc') or read('A-384x64-smooth')
def packet_plot(x,y,style,color,label,alpha=1):
    # Use the same interpolation coordinate as infer_thresholds.py.
    if len(x)>1:
        dense=np.geomspace(min(x),max(x),200)
        yy=np.interp(1/dense,(1/x)[::-1],y[::-1])
        ax[1].semilogx(dense,yy,'--' if '--' in style else '-',color=color,alpha=alpha)
    ax[1].semilogx(x,y,style[0],color=color,label=label,alpha=alpha)
if smooth:
    rows=[]
    for k in [4,5,6,8,16,32]:
        d=read(f'A-384x128-packets-k{k}')
        if d:rows.append((meridional(k)['circular_cap_equivalent_N'],100*output(d)/output(smooth)))
    if rows:
        x,y=np.array(rows).T;packet_plot(x,y,'o-','#2874a6','2D sector equivalents')
for grid,style in [(32,'s--'),(48,'s-'),(64,'D-')]:
    c=read(f'A-cart{grid}-smooth')
    if not c:continue
    rows=[]
    for n in [32,64,96,128,256]:
        d=read(f'A-cart{grid}-packets-n{n}')
        if d:rows.append((n,100*output(d)/output(c)))
    if rows:
        x,y=np.array(rows).T;packet_plot(x,y,style,'#b03a2e',f'3D packets, base {grid}',alpha=.7 if grid==32 else 1)
ax[1].set(xlabel='Packet count (2D points use geometric equivalents)',ylabel='Nuclear output retained (%)',title='One bulk direction per shaped packet')
if len(ax[1].lines):ax[1].legend(frameon=False,fontsize=8)
for a in ax:
    a.axhline(90,color='#a93226',linestyle='--',linewidth=1)
    a.grid(alpha=.2)
    minimum=min([float(np.min(line.get_ydata())) for line in a.lines if len(line.get_ydata())]+[60.])
    a.set_ylim(max(0,min(60,10*np.floor(minimum/10))),103)
fig.suptitle('Prepared 30 m-reference nitrogen/proton pulse',fontsize=13)
fig.savefig(out/'tolerance-curves.png',dpi=180);fig.savefig(out/'tolerance-curves.svg');plt.close(fig)

c=read('A-768-1d-preheat-f0')
if c:
    fig,a=plt.subplots(figsize=(6.5,4.4),layout='constrained')
    for mode,label,color in [('preheat','Added merger heat','#d35400'),('preheatfloor','Existing preheat credited','#239b56')]:
        rows=[(0,100)]
        for f in [.01,.02,.025,.03,.05,.075,.1]:
            d=read(f'A-768-1d-{mode}-f{f}')
            if d:rows.append((100*f,100*output(d)/output(c)))
        if len(rows)>1:
            x,y=np.array(rows).T;a.plot(x,y,'o-',color=color,label=label)
    if smooth:
        for k in [4,5,6,8,16,32]:
            d=read(f'A-384x128-packets-k{k}')
            if d:
                a.plot(100*meridional(k)['direction_sine_rms']**2,100*output(d)/output(smooth),'s',color='#2874a6',label='2D packet merger' if k==4 else None)
    a.axhline(90,color='#a93226',ls='--',lw=1);a.grid(alpha=.2)
    a.set(xlabel='Kinetic-energy fraction assigned to early heat (%)',ylabel='Nuclear output retained (%)',title='Energy conserved; compression energy becomes heat',xlim=(0,10),ylim=(65,102))
    a.legend(frameon=False,fontsize=9)
    fig.savefig(out/'early-heat-sensitivity.png',dpi=180);fig.savefig(out/'early-heat-sensitivity.svg');plt.close(fig)

names=['A-384x128-smooth-hllc','A-384x128-packets-k4','A-384x128-packets-k8']
titles=['Smooth pulse','About 40 packets (sector equivalent)','About 160 packets (sector equivalent)']
available=[]
for name,title in zip(names,titles):
    d=base/'runs'/name
    groups={}
    for f in d.glob('*.vtk'):groups.setdefault(f.name.split('.')[-2],[]).append(f)
    if not groups:continue
    candidates=[]
    for files in groups.values():
        # Header-only scan; do not repeatedly load whole fields to select a time.
        import re
        head=files[0].read_bytes()[:250].decode('ascii',errors='replace')
        t=float(re.search(r'time=\s*([\deE+.-]+)',head).group(1))
        candidates.append((abs(t*30-27),files))
    _,files=min(candidates,key=lambda x:x[0]);available.append((title,snapshot(files)))
if available:
    fig,axes=plt.subplots(1,len(available),figsize=(4*len(available),4.4),layout='constrained',squeeze=False)
    for a,(title,(t,r,theta,d)) in zip(axes[0],available):
        rr,tt=np.meshgrid(r*30,theta);x=rr*np.sin(tt);z=rr*np.cos(tt)
        rho=d['dens']*1e8
        im=a.pcolormesh(x,z,rho,norm=LogNorm(1e6,3e8),cmap='viridis',shading='flat',rasterized=True)
        a.pcolormesh(-x,z,rho,norm=LogNorm(1e6,3e8),cmap='viridis',shading='flat',rasterized=True)
        a.set(xlim=(-65,65),ylim=(-65,65),aspect='equal',xlabel='Radius across axis (m)',ylabel='Position along axis (m)',title=f'{title}\n{t*30:.1f} microseconds')
    fig.colorbar(im,ax=list(axes[0]),label='Density (kg/m³)',shrink=.8)
    fig.savefig(out/'merger-density.png',dpi=180);fig.savefig(out/'merger-density.svg');plt.close(fig)
print('Figures saved to',out)
