"""Static scientific figure from the saved finite open-loop station scan.

No model is evaluated or changed. Run with the existing matplotlib runtime:
PYTHONPATH=/tmp/protium-python python plot_station_sensitivity.py
Inputs default to station.json and station-inputs.json beside this script.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os

os.environ.setdefault('MPLCONFIGDIR','/var/tmp/tkwa-chatgpt/protium-station-figure/mpl')
Path(os.environ['MPLCONFIGDIR']).mkdir(parents=True,exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--input',type=Path,default=HERE/'station.json')
p.add_argument('--assumptions',type=Path,default=HERE/'station-inputs.json')
p.add_argument('--output-directory',type=Path,default=HERE)
a=p.parse_args()
data=json.loads(a.input.read_text());assumptions=json.loads(a.assumptions.read_text())
input_hash=hashlib.sha256(a.input.read_bytes()).hexdigest()
assumption_hash=hashlib.sha256(a.assumptions.read_bytes()).hexdigest()
assert data['inputs_sha256']==assumption_hash,'Station scan does not match current assumptions'
cycle=assumptions['cycle']
assert math.isclose(cycle['total_ongoing_C_H_mass_budget'],.007,abs_tol=1e-12)
assert math.isclose(cycle['fresh_C_H_mass_ratio'],.0056,abs_tol=1e-12)
assert math.isclose(cycle['eta_driver'],.6) and math.isclose(cycle['eta_heat'],.4)
assert math.isclose(cycle['energy_capture_fraction'],.99)
rows=[r for r in data['family'] if math.isclose(r['absorber_sigma_kg_m2'],5000)]
services=sorted({r['service_sigma_kg_m2'] for r in rows})
ratings=sorted({r['peak_driver_rating_W_kg'] for r in rows})
timings=['mass999','full_resolved']
assert len(services)==3
styles=[('#087E8B','o'),('#CC6414','s'),('#7254A3','^')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
 'axes.titlesize':11,'axes.titleweight':'bold','axes.spines.top':False,
 'axes.spines.right':False,'axes.edgecolor':'#87939D','xtick.color':'#36434D',
 'ytick.color':'#36434D','text.color':'#26333E','axes.labelcolor':'#26333E',
 'svg.fonttype':'none','svg.hashsalt':'protium-open-loop-station-sensitivity-v1',
 'savefig.facecolor':'white'})
fig,axes=plt.subplots(1,2,figsize=(12.3,7.5),sharex=True,sharey=True)
fig.subplots_adjust(left=.082,right=.975,top=.735,bottom=.285,wspace=.16)
fig.suptitle('Open-loop station: source and service mass sensitivity',x=.5,y=.965,fontsize=16,weight='bold')
fig.text(.5,.91,'5,000 kg/m² absorber · two separate pulse-clearance assumptions',ha='center',fontsize=11,color='#4E5C68')
handles=[Line2D([0],[0],color=color,marker=marker,lw=2,markersize=6,
 label=f'{service:,.0f} kg/m²') for service,(color,marker) in zip(services,styles)]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.872),ncol=3,
 frameon=False,title='Processing, structure and thermal service mass per area',title_fontsize=10,
 columnspacing=2.8,handlelength=2.5)
plot_points=[]
all_y=[r['best']['net_W_kg']*1000 for r in rows if r['best'] is not None]
ymax=max(4.,math.ceil(max(all_y)*2)/2+.25)
for ax,timing,title in zip(axes,timings,['(a) 99.9% mass window','(b) All resolved shells']):
 ax.set_title(title,pad=14,loc='left')
 ax.set_xscale('log');ax.set_xlim(min(ratings)*.78,max(ratings)*1.70)
 ax.set_ylim(0,ymax)
 ax.xaxis.set_major_locator(FixedLocator(ratings))
 ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'))
 ax.xaxis.set_minor_locator(NullLocator())
 ax.yaxis.set_major_locator(FixedLocator([x*.5 for x in range(int(ymax*2)+1)]))
 ax.grid(axis='y',color='#E2E7EB',linewidth=.7,zorder=0)
 ax.set_axisbelow(True);ax.tick_params(length=3)
 ax.set_xlabel('Peak source rating (W/kg of source hardware)',labelpad=12,fontsize=10)
 ax.axhline(1,color='#7C8891',lw=1,ls=(0,(4,4)),zorder=1)
 feasible_count=0
 for service,(color,marker) in zip(services,styles):
  chosen=sorted([r for r in rows if r['timing_basis']==timing and r['service_sigma_kg_m2']==service],key=lambda r:r['peak_driver_rating_W_kg'])
  assert len(chosen)==len(ratings)
  xs=[];ys=[]
  for row in chosen:
   best=row['best'];rating=row['peak_driver_rating_W_kg']
   value=best['net_W_kg']*1000 if best else None
   if best:
    assert best['feasible'] and best['timing_basis']==timing
    assert math.isclose(best['absorber_sigma_kg_m2'],5000)
    feasible_count+=1
   xs.append(rating);ys.append(math.nan if value is None else value)
   plot_points.append(dict(timing_basis=timing,service_sigma_kg_m2=service,peak_source_rating_W_kg=rating,net_mW_kg=value,feasible=best is not None))
  ax.plot(xs,ys,color=color,marker=marker,lw=2,markersize=6.5,zorder=3)
  valid=[(x,y) for x,y in zip(xs,ys) if math.isfinite(y)]
  if valid:
   x,y=valid[-1]
   ax.annotate(f'{y:.3f}',(x,y),xytext=(7,8),textcoords='offset points',color=color,
               fontsize=10,fontweight='bold',ha='left',va='bottom')
 if timing=='mass999':
  ax.text(14,1.055,'1 mW/kg reference',color='#6B7780',fontsize=9)
  ax.text(.03,.93,'Tail handling required beyond this window',transform=ax.transAxes,
          color='#687681',fontsize=9,va='top')
 else:
  total=sum(r['timing_basis']==timing for r in rows)
  missing=total-feasible_count
  ax.text(.48,.63,f'{missing} of {total} allocations have no feasible point\nin the tested family',
          transform=ax.transAxes,ha='center',va='center',fontsize=11,color='#596976',
          linespacing=1.5,bbox=dict(boxstyle='round,pad=.75',facecolor='#F3F6F8',edgecolor='none'))
  ax.text(.48,.39,'Only feasible points are shown.',transform=ax.transAxes,
          ha='center',va='center',fontsize=9,color='#687681')
axes[0].set_ylabel('Net electrical specific power (mW/kg)',labelpad=10)
fig.text(.082,.185,'Conditional engineering allocations — not measured designs.',fontsize=11,weight='bold')
fig.text(.082,.139,'Fresh C:H mass budget: 0.007 kg/kg total = 0.0056 fuel-process makeup + 0.0014 net material makeup.',fontsize=10)
fig.text(.082,.099,'60% preparation efficiency · 40% heat conversion · 99% energy capture · Supplied feed; extraction and delivery excluded.',fontsize=9.3,color='#52616C')
fig.text(.082,.055,'Markers show the best feasible points in a finite scan; connecting lines guide the eye. Corrected 64-zone prescribed pulses.',fontsize=9,color='#52616C')
a.output_directory.mkdir(parents=True,exist_ok=True)
svg=a.output_directory/'station-sensitivity.svg';png=a.output_directory/'station-sensitivity.png'
fig.savefig(svg,metadata={'Date':None,'Description':'Conditional engineering allocations from station.json; no measured station performance.'})
fig.savefig(png,dpi=180,metadata={'Software':'Matplotlib; plot_station_sensitivity.py'})
plt.close(fig)
metadata=dict(title='Open-loop station: source and service mass sensitivity',scope='Static display of saved finite conditional allocations; no model evaluation.',station_sha256=input_hash,assumptions_sha256=assumption_hash,plot_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),absorber_sigma_kg_m2=5000,total_C_H_mass_budget=.007,fuel_process_C_H_mass_budget=.0056,omitted_infeasible_points=sum(not r['feasible'] for r in plot_points),points=plot_points,output_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [svg,png]})
(a.output_directory/'station-sensitivity.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(outputs=[str(svg),str(png)],feasible_markers=sum(p['feasible'] for p in plot_points),omitted_points=metadata['omitted_infeasible_points']),indent=2))
