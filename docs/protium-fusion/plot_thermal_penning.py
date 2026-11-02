from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from penning_equilibrium import equilibrium,at_wall,MP,MEV,E_CHARGE
ROOT=Path(__file__).parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'protium-thermal-penning-v1'})
a=json.loads((ROOT/'thermal-neutron-screen.json').read_text())['cases']
b=json.loads((ROOT/'penning-scan.json').read_text())['best'][0]
c=json.loads((ROOT/'thermal-recovery-sensitivity.json').read_text())['scenarios']
eq=equilibrium(b['h0'],b['c0'],b['xi'],b['xwall'],nout=1501)
fig,axs=plt.subplots(2,2,figsize=(12,8.5));fig.subplots_adjust(left=.075,right=.97,bottom=.12,top=.89,hspace=.43,wspace=.28)
fig.suptitle('Thermal carbon-13 neutron catalysis: conditional screens',fontsize=17,x=.075,ha='left',y=.97)
fig.text(.075,.926,'Common ion temperature · electron-free plasma · external neutron capture and deuterium fusion',fontsize=11,color='#555555')
ax=axs[0,0];t=np.array([v['ion_temperature_MeV'] for v in a]);use=(t>=.3)&(t<=1.2)
for idx,label,col in [(1,'80% heating / 90% ion recovery; free regeneration','#a8a8a8'),(2,'80% / 90%; 50 MeV per side product','#cf7441'),(3,'95% / 95%; 50 MeV per side product','#296a91')]:
 y=np.array([v['energy_cases'][idx]['net_MeV_per_pn'] for v in a]);ax.plot(t[use],y[use],'-o',ms=3,label=label,color=col)
ax.axhline(0,color='#777777',lw=.8);ax.set(ylim=(-15,6),xlabel='Ion temperature (MeV)',ylabel='Net energy per new neutron (MeV)',title='A. Nuclear process before trap costs')
ax.legend(fontsize=8,loc='lower right',frameon=False);ax.text(.98,.95,'40% heat conversion',ha='right',va='top',transform=ax.transAxes,color='#555555',fontsize=9)
ax=axs[0,1];l=np.sqrt(.6*MEV*MP)/(E_CHARGE*100)
ax.semilogy(eq['x']*l*1000,eq['h']/b['h0'],label='Proton density / central proton density',color='#296a91')
ax.semilogy(eq['x']*l*1000,eq['c']/b['c0'],label='Carbon density / central carbon density',color='#ad5939')
ax.axvline(b['radius_m']*1000,color='#777777',ls='--',lw=.8)
ax.set(ylim=(1e-20,5),xlabel='Radius (mm)',ylabel='Density relative to its axis value',title='B. One rotating equilibrium at 100 T')
ax.legend(fontsize=8,loc='lower left',frameon=False);ax.text(.95,.95,'0.6 MeV; wall at 18.2 mm',ha='right',va='top',transform=ax.transAxes,fontsize=9,color='#555555')
ax=axs[1,0]
colors={0.:'#356d97',10.:'#539871',50.:'#bb7933',230.:'#ad5266'}
for side in [0.,10.,50.,230.]:
 vals=[r for r in c if r['side_cost_MeV']==side]
 ax.plot([r['heat_efficiency']*100 for r in vals],[r['best']['field_mass_only_net_W_kg']*1000 for r in vals],'-o',label=f'{side:g} MeV / side product',color=colors[side],ms=4)
 ax.text(81,vals[-1]['best']['field_mass_only_net_W_kg']*1000,f'{side:g} MeV',va='center',color=colors[side],fontsize=8)
ax.axhline(1,color='#111111',ls='--',lw=1);ax.text(40.5,1.25,'Whole-plant target: 1 mW/kg',fontsize=9)
ax.set_yscale('symlog',linthresh=.003);ax.set(xlim=(38,94),ylim=(-.015,3),xlabel='Heat conversion efficiency (%)',ylabel='Field-structure output (mW/kg)',title='C. Regeneration and heat recovery matter')
ax.set_xticks([40,60,70,80]);ax.set_yticks([-.01,0,.01,.1,1],['−0.01','0','0.01','0.1','1'])
ax.axhline(0,color='#777777',lw=.7);ax.text(.02,.03,'Labels give net cost per side product',transform=ax.transAxes,fontsize=8,color='#555555')
ax=axs[1,1];sens=json.loads((ROOT/'thermal-penning-sensitivity.json').read_text())['cases']
for pot,col in [(1.,'#296a91'),(.999,'#539871'),(.99,'#ad5939')]:
 vals=[]
 for B in [50,100,200,400]:
  r=max([v for v in sens if v['B_T']==B and v['potential_in_efficiency']==pot],key=lambda v:v['field_mass_only_net_W_kg'])
  vals.append(r['field_mass_only_net_W_kg']*1000)
 ax.loglog([50,100,200,400],vals,'-o',ms=4,color=col,label=f'{100*pot:g}% potential-energy transfer')
ax.axhline(1,color='#111111',ls='--',lw=1);ax.set(ylim=(.0003,2),xlabel='Magnetic field (T)',ylabel='Field-structure output (mW/kg)',title='D. Field extrapolation needs a real coil')
ax.set_xticks([50,100,200,400],['50','100','200','400']);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());ax.legend(fontsize=8,loc='lower right',frameon=False)
ax.text(.02,.82,'40% heat; 50 MeV side-product cost',transform=ax.transAxes,va='top',fontsize=9,color='#555555')
fig.text(.075,.043,'C–D: 95% heating and ion recovery; specific structural strength 1 MJ/kg; mass factor 5; 70% column packing.\nField mass excludes the neutron blanket, deuterium plant, end plugs, converters and cooling. No complete reactor is established.',fontsize=9,color='#444444')
for ext in ['png','svg']:
 p=ROOT/f'thermal-penning.{ext}';fig.savefig(p,dpi=180,metadata={'Date':None} if ext=='svg' else None)
 if ext=='svg':p.write_text('\n'.join(v.rstrip() for v in p.read_text().splitlines())+'\n')
