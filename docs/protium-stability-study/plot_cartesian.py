"""Central slices of the 3D control and packet calculations on their native AMR mesh."""
import json,re,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm

base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
sys.path.insert(0,str(base/'athena/vis/python'));import athena_read
out=Path(__file__).resolve().parent/'figures';out.mkdir(exist_ok=True)
panels=[]
for name,label in [('A-cart48-smooth','Smooth pulse'),('A-cart48-packets-n64','64 shaped packets'),('A-cart48-packets-n128','128 shaped packets')]:
    d=base/'runs'/name
    if not (d/'execution.json').exists():continue
    e=json.loads((d/'execution.json').read_text())
    if not e.get('reached_time_limit') or e.get('fatal_error'):continue
    groups={}
    for f in d.glob('*.vtk'):groups.setdefault(f.name.split('.')[-2],[]).append(f)
    choices=[]
    for files in groups.values():
        with files[0].open('rb') as f:head=f.read(300).decode('ascii',errors='replace')
        t=float(re.search(r'time=\s*([\deE+.-]+)',head).group(1))
        choices.append((abs(t*30-31.5),t,files))
    if not choices:continue
    _,t,files=min(choices,key=lambda x:x[0]);polygons=[];values=[]
    for f in files:
        x,y,z,fields=athena_read.vtk(f)
        if not z[0]<=0<z[-1]:continue
        k=np.searchsorted(z,0,side='right')-1
        x=x*30;y=y*30
        for j in range(len(y)-1):
            if y[j]>65 or y[j+1]<-65:continue
            for i in range(len(x)-1):
                if x[i]>65 or x[i+1]<-65:continue
                polygons.append([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])
                values.append(fields['dens'][k,j,i]*1e8)
    panels.append((label,t,polygons,values))
if panels:
    fig,axes=plt.subplots(1,len(panels),figsize=(4*len(panels),4.3),layout='constrained',squeeze=False)
    for a,(label,t,polygons,values) in zip(axes[0],panels):
        pc=PolyCollection(polygons,array=np.asarray(values),norm=LogNorm(1e6,3e8),cmap='viridis',edgecolor='none',rasterized=True)
        a.add_collection(pc);a.set(xlim=(-65,65),ylim=(-65,65),aspect='equal',xlabel='x (m)',ylabel='y (m)',title=f'{label}\n{30*t:.1f} microseconds')
    fig.colorbar(pc,ax=list(axes[0]),label='Density in central slice (kg/m³)',shrink=.8)
    fig.savefig(out/'three-dimensional-merger.png',dpi=180)
    fig.savefig(out/'three-dimensional-merger.svg')
    print(out/'three-dimensional-merger.png')
