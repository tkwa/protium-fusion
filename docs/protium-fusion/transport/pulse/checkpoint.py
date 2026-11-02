"""Produce compact resolution comparison and a scientific diagnostic figure."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('--out',default=str(Path(__file__).parent));a=ap.parse_args();out=Path(a.out)
main=json.loads((out/'pulse-transport-audit.json').read_text())
resolution=json.loads((out/'resolution-audit.json').read_text())
rows=[]
for coarse in resolution['rows']:
 fine=next(x for x in main['rows'] if x['stage']==coarse['stage'])
 keys=['conductive_sum_absolute_interface_transfer_over_nuclear','conductive_redistributed_energy_over_nuclear','conductive_max_interface_transfer_over_nuclear']
 rows.append({'stage':fine['stage'],'zones':[coarse['zones'],fine['zones']],
              **{k:{'coarse':coarse[k],'fine':fine[k],'fine_over_coarse':fine[k]/coarse[k]} for k in keys}})
(out/'resolution-comparison.json').write_text(json.dumps({'scope':'Frozen-trajectory conduction diagnostics. The interface sum and absolute divergence are resolution dependent; no converged transport correction is claimed.','rows':rows},indent=2)+'\n')
fig,axs=plt.subplots(1,3,figsize=(12,4.1),layout='constrained')
y=np.arange(len(main['rows']))
fields=[('Poutanen_kappa_m2_kg_within_domain',1,'Compton transport opacity\n(dashed: prescribed grey value)','m²/kg'),
        ('effective_optical_depth_per_cell_x5',.01,'Thermalization depth after 100×\nreduction of both coefficients','effective depth per cell, hν = 5kT'),
        ('conductive_diffusion_length_over_temperature_scale',10,'Conductive smoothing after\n100× increase in conductivity','diffusion length / temperature scale')]
for ax,(key,mult,title,label) in zip(axs,fields):
 data=np.array([r['burn_weighted_q01_q50_q99'][key] for r in main['rows']])*mult
 ax.errorbar(data[:,1],y,xerr=np.vstack([data[:,1]-data[:,0],data[:,2]-data[:,1]]),fmt='o',color='#175d87',capsize=4)
 ax.set_xscale('log');ax.set_yticks(y,[r['stage'] for r in main['rows']]);ax.invert_yaxis();ax.grid(axis='x',alpha=.25)
 ax.set_title(title,fontsize=10);ax.set_xlabel(label,fontsize=9)
axs[0].axvline(.001,color='#b64e22',ls='--',label='prescribed grey value')
axs[1].axvline(1,color='#b64e22',ls='--');axs[2].axvline(1,color='#b64e22',ls='--')
fig.suptitle('Saved spherical pulses: local microphysics during energy release',fontsize=13)
fig.savefig(out/'transport-margins.png',dpi=180)
