"""Render the bounded axial-drive discriminator; no full-reactor claim."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
HERE=Path(__file__).parent
main=json.loads((HERE/'axial-response-main.json').read_text());validation=json.loads((HERE/'axial-validation.json').read_text());capture=json.loads((HERE/'axial_capture.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(12,5.3))
fig.subplots_adjust(left=.075,right=.99,bottom=.23,top=.84,wspace=.20)
colors={40:'#2263a3',80:'#bd541b',160:'#238565'}
for k,color in colors.items():
 r=[r for r in main['rows'] if r['k_per_m']==k and r['beta_phase']<=.3]
 axs[0].plot([r['beta_phase'] for r in r],[r['median_reaction_transfer_abs'] for r in r],marker='.',color=color,label=f'k = {k} /m')
 r=[r for r in validation['derived_capture_costs'] if r['k_per_m']==k]
 axs[1].semilogy([r['beta_phase'] for r in r],[r['continuous_heat_MeV_per_primary'] for r in r],marker='o',color=color,label=f'k = {k} /m')
r=[r for r in capture['rows'] if r['k_per_m']==40 and .1<=r['beta_phase']<=.3]
axs[1].semilogy([r['beta_phase'] for r in r],[r['uniform_ion_heat_at_threshold_MeV_per_primary'] for r in r],color='#777777',linestyle='--',label='Uniform total potential, k = 40 /m')
axs[1].axhline(1.117140444,color='#282828',linestyle=':',label='Extra electrical margin (before heating credit)')
axs[0].set(xlabel='Wave phase velocity / c',ylabel='Interior potential / wall potential',title='Self-consistent linear field at median reaction radius',ylim=(0,2.2))
axs[1].set(xlabel='Wave phase velocity / c',ylabel='Initial ion heating [MeV per primary event]',title='Continuous field at cold-lepton trapping threshold',ylim=(.01,1e8))
for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best');ax.title.set_fontsize(11)
fig.suptitle('Axial waves: fast parallel motion does not remove radial ion resonances',fontsize=13)
fig.text(.5,.035,'Collisionless electrostatic model. Threshold admits zero phase measure; half-phase availability requires twice the voltage.\nHeating is not a net electrical cost; nonlinear saturation, pulse duty, usable heating and Maxwell fields remain separate.',ha='center',fontsize=9)
for ext in ['png','svg']:fig.savefig(HERE/('axial-summary.'+ext),dpi=180,bbox_inches='tight')
