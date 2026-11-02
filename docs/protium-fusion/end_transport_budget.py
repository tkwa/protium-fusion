"""Translate specified additional end-system loads into complete-plant output.

Electron-kinetic cases charge the calculated co-moving electron kinetic energy
as an additional unrecovered electrical load. Injector, neutralizer and sorting
hardware remain unspecified; these cases are sensitivities, not plant designs.
"""
from pathlib import Path
import argparse,json,gzip
from network_plant import plant
ROOT=Path(__file__).parent
MEV=1.602176634e-13

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--screen-dir',type=Path,default=ROOT);args=ap.parse_args()
    energy=json.loads((args.root/'network-ion-budget.json').read_text())
    family=json.loads(gzip.decompress((args.root/'network-channel-family.json.gz').read_bytes()))
    screen=json.loads((args.screen_dir/'end-transport-screen.json').read_text())
    scale=screen['reference_primary_rate_s']*MEV;nch=screen['power_allowances']['total_channels']
    cases=[]
    variants=[(f'Additional {w:g} W per channel',w*nch/scale) for w in [0.,50.,100.,250.,screen['power_allowances']['extra_electrical_W_per_channel_at_1mW_kg'],1000.]]
    variants += [(f'Unrecovered co-moving electron kinetic energy at {c["energy_per_charge_V"]/1e6:g} MV per charge',c['co_moving_electron_kinetic_MeV_per_primary'])
                 for c in screen['transport_cases'] if c['energy_per_charge_V']>=40e6]
    for name,load in variants:
        r=plant(energy,family,extra_processing_MeV_per_primary=load)
        cases.append({'variant':name,'additional_MeV_per_primary':load,'additional_W_per_channel':load*scale/nch,
          'net_W_kg':r['net_W_kg'],'net_W':r['net_W'],'mass_kg':r['total_mass_kg'],
          'extra_processing_W':r['extra_processing_W'],'mass_ledger_error_kg':r['mass_ledger_error_kg'],
          'energy_ledger_error_W':r['energy_ledger_error_W']})
        print(name, 'mW/kg',r['net_W_kg']*1e3)
    (args.screen_dir/'end-transport-budget.json').write_text(json.dumps({'scope':__doc__,'cases':cases},indent=2,allow_nan=False)+'\n')
