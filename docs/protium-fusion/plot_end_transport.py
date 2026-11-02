"""Plot the conditional end-transport screen; no extraction design is implied."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);args=ap.parse_args()
    s=json.loads((args.root/'end-transport-screen.json').read_text())
    budget=json.loads((args.root/'end-transport-budget.json').read_text())
    cases=s['transport_cases'];u=np.array([c['energy_per_charge_V'] for c in cases])/1e6
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    ax=axes[0,0]
    ax.loglog(u,[c['long_merged_beam_axis_potential_V'] for c in cases],'o-',label='All streams merged, no electrons',color='#ad493d')
    ax.loglog(u,[c['channel_details'][0]['channel_potential_V'] for c in cases],'o-',label='One screened titanium channel',color='#266c8b')
    ax.set(xlabel='Ion transport energy per charge (MV)',ylabel='Long-cylinder axis potential (V)',title='20.2 MA total; 16.3 mA per titanium channel')
    ax.legend(fontsize=8,loc='center right');ax.grid(alpha=.2)
    ax=axes[0,1]
    ax.loglog(u,[c['channel_details'][0]['sufficient_planar_gap_m_for_10percent_spacecharge'] for c in cases],'o-',color='#266c8b')
    ax.set(xlabel='Minimum ion energy per charge (MV)',ylabel='Sufficient planar gap length (m)',title='Titanium stream: space-charge depression ≤ 10%')
    ax.text(.04,.08,'Prescribed electrodes; transverse optics remain open.',transform=ax.transAxes,fontsize=9)
    ax.grid(alpha=.2)
    ax=axes[1,0];hi=[c for c in cases if c['energy_per_charge_V']>=40e6]
    ax.bar([f'{c["energy_per_charge_V"]/1e6:g} MV' for c in hi],[c['co_moving_electron_kinetic_MeV_per_primary'] for c in hi],color='#447f69',width=.5)
    margin=s['power_allowances']['extra_electrical_MeV_per_primary_at_1mW_kg']
    ax.axhline(margin,color='#ad493d',linestyle='--',label='Remaining extra electrical allowance')
    ax.set(ylabel='MeV per primary titanium reaction',title='External co-moving electrons: kinetic energy only',ylim=(0,1.35))
    ax.legend(fontsize=8,loc='upper center')
    ax.text(.03,.43,'Neutralizer hardware, matching and instabilities\nare additional requirements.',transform=ax.transAxes,fontsize=9)
    ax=axes[1,1];r=budget['cases'][:6]
    ax.plot([c['additional_W_per_channel'] for c in r],[c['net_W_kg']*1e3 for c in r],'o-',color='#266c8b')
    ax.axhline(1,color='#ad493d',linestyle='--',label='Requested target')
    ax.axvline(s['power_allowances']['extra_electrical_W_per_channel_at_1mW_kg'],color='#666',linestyle=':',label='475 W/channel allowance')
    ax.set(xlabel='Additional load averaged over all channels (W)',ylabel='Net useful output (mW/kg)',title='Complete mass, waste heat and replacement recalculated',ylim=(0,1.65))
    ax.legend(fontsize=8,loc='upper right');ax.grid(alpha=.2)
    fig.suptitle('Conditional titanium reactor: finite end-transport constraints',fontsize=15)
    for ext in ['png','svg']:fig.savefig(args.root/f'end-transport-summary.{ext}',dpi=170)
