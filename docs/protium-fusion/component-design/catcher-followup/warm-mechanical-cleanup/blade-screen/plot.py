from pathlib import Path
import json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
rows=json.loads((p/'records/grid.json').read_text())['cases']
colors={0.:'#1768ac',.5:'#d97706',1.:'#8b4ca0'}
fig,axs=plt.subplots(2,2,figsize=(10.5,7.4),sharex=True)
for i,gas in enumerate(['He','H2']):
 for a in [0.,.5,1.]:
  for L,style in [(1.,'-'),(3.,'--')]:
   rr=sorted([r for r in rows if r['parameters']['gas']==gas and r['parameters']['alpha']==a and r['parameters']['L_pitch']==L],key=lambda r:r['parameters']['theta'])
   x=[r['parameters']['theta'] for r in rr]
   for j in [0,1]:
    y=[r['K0'] if j==0 else r['net_long_cascade_capture'] for r in rr]
    if j==0:
     e=[r['K0']*1.95996398454*math.sqrt((1-r['forward']['T'])/(100000*r['forward']['T'])+(1-r['reverse']['T'])/(100000*r['reverse']['T'])) for r in rr]
    else:e=[(r['net_capture_normal95'][1]-r['net_capture_normal95'][0])/2 for r in rr]
    axs[i,j].errorbar(x,y,yerr=e,color=colors[a],linestyle=style,marker='o' if L==1 else None,markersize=3.5,capsize=2,linewidth=1.5,alpha=1 if L==1 else .6,label=f'diffuse probability {a:g}' if L==1 else None)
 for j in [0,1]:
  axs[i,j].grid(alpha=.2);axs[i,j].set_title(gas+(' — one-row compression' if j==0 else ' — conditional long-cascade capture'),fontsize=11)
  axs[i,j].set_xlim(18,62);axs[i,j].set_xticks([20,30,45,60])
 axs[i,0].set_ylabel('K₀ = Tforward / Treverse')
 axs[i,1].set_ylabel('open fraction × (Tforward − Treverse)')
 axs[i,1].axhline(.1034 if gas=='He' else .0622,color='#555555',linestyle=':',linewidth=1.2,label='existing capture target')
 axs[i,0].axhline(1,color='#aaaaaa',linewidth=.7)
for ax in axs[1,:]:ax.set_xlabel('Blade tilt from pump axis (degrees)')
axs[0,0].legend(fontsize=8,loc='upper left')
for i,target in enumerate([.1034,.0622]): axs[i,1].text(60,target+.004,'capture reference',ha='right',fontsize=8,color='#555555')
fig.suptitle('Moving angled blades at 700 m/s and 1000 K',fontsize=15,y=.98)
fig.text(.5,.925,'Gap/pitch = 0.65. Solid: axial length/pitch 1. Dashed: 3. Error bars: 95% Monte Carlo sampling.',ha='center',fontsize=10)
fig.text(.5,.025,'100,000 molecules per direction per point. Stationary Maxwellian interstage reservoirs are assumed for cascade capture.\nAn actual stator/thermalization and bypass design is still required; these are unit-cell results.',ha='center',fontsize=9)
fig.subplots_adjust(top=.86,bottom=.14,hspace=.28,wspace=.27)
for ext in ['png','pdf','svg']:fig.savefig(p/('blade-screen.'+ext),dpi=180)
