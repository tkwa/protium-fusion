from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
x=json.loads((HERE/'delivery-screen.json').read_text())
a=x['cases'][0]
fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
ax=axes[0]
eps=np.logspace(-6,0,200)
rho=1000.;h=1e6;time=2e10/a['energy_rms_speed_m_s']
for temp,color in [(300,'#326b9c'),(500,'#c98b24'),(1000,'#ba454b')]:
    amin=3*eps*5.670374419e-8*temp**4*time/(rho*h)
    ax.loglog(eps,amin,label=f'{temp} K cavity',color=color,lw=2)
for j,ls in [(0,':'),(1,'--')]:
    p=a['cold_pellet_acceleration_survival_trade'][j]
    ax.axhline(p['maximum_accelerable_pellet_radius_m'],color='#333333',ls=ls,lw=1.7,
               label=f"{p['accelerator_gradient_V_m']/1e6:g} MV/m × {p['accelerator_length_m']/1e3:g} km")
ax.set(xlabel='Planck-mean absorption efficiency',ylabel='Pellet radius (m)',
       title='Cold-flight size versus acceleration limit',ylim=(1e-6,10))
ax.text(.02,.97,'N14 + p; 1 MJ/kg heat allowance\n1 GV/m grain surface field; ρ = 1000 kg/m³',
        transform=ax.transAxes,va='top',fontsize=9)
ax.legend(loc='lower right',fontsize=8,frameon=False)
ax.grid(alpha=.2,which='both')
ax=axes[1]
for index,label,color in [(0,'N14 + p','#326b9c'),(2,'N15 + p','#c98b24'),(5,'He factory','#ba454b')]:
    c=x['cases'][index]
    vals=c['launch_schedule']['mass_weighted_flight_time_quantiles_s']
    probs=np.array([float(k) for k in vals])
    times=np.array(list(vals.values()))/3600
    ax.plot(times,100*probs,'o-',color=color,label=label,lw=1.7)
ax.set(xlabel='Flight time before prescribed incoming state (hours)',
       ylabel='Cumulative fraction of pulse mass (%)',title='Kinematic launch schedule',ylim=(0,100))
ax.grid(alpha=.2)
ax.legend(frameon=False,loc='lower right')
fig.suptitle('Opaque-pulse driver: necessary delivery requirements',fontsize=14)
fig.savefig(HERE/'delivery-requirements.png',dpi=180)
fig.savefig(HERE/'delivery-requirements.svg')
