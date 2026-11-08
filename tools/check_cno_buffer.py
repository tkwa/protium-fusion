#!/usr/bin/env python3
"""Replay the inherited linear buffer model on the selected large-N15 calendar.

All primary pump, conduit and cold-panel speeds scale 16-fold; buffer volume
stays 0.5 m3/m2 at 1000 K. This is a capacity screen, with apparatus mass and
surface qualification remaining acceptance requirements of the plant account.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'docs/protium-fusion/second-generation'
F=ROOT/'docs/protium-fusion/component-design/catcher-followup'
PAGE=ROOT/'public/blog/ordinary-hydrogen-fusion'


def main():
    src=F/'warm-mechanical-cleanup/buffer.py'
    ns={'__file__':str(src)}
    definitions=src.read_text().split('outrows=[];')[0]
    definitions=definitions.replace("Path('/home/tkwa/.codex/worktrees/protium-fusion/website')",'Path('+repr(str(ROOT))+')')
    exec(compile(definitions,str(src),'exec'),ns)
    x=json.loads((P/'publication-inputs.json').read_text());d=ns['d'];rec=x['record'];st=x['probe']['stage']
    z=np.load(PAGE/'plasma-simulation/data/n15.npz');names=z['species'].tolist();nn=len(names)
    state=z['state'].reshape(int(z['zones']),int(z['variables_per_zone']),-1)[:,:,-1]
    mass=z['cell_mass_kg'];y=np.maximum(state[:,:nn],0.)
    a=np.array([1 if n in ('p','n') else int(''.join(filter(str.isdigit,n))) for n in names])
    y/=(y@a)[:,None]
    new=d['stages']['n15']
    new.update(radius_bounds_m=np.r_[0,state[:,nn+1]*rec['R0_m']].tolist(),
               velocity_bounds_m_s=np.r_[0,state[:,nn+2]*float(z['velocity_unit_cm_s'])/100].tolist(),
               endpoint_time_s=rec['elapsed_s'],pulse_mass_kg=rec['mass_kg'])
    for n,key in [('p','H'),('he4','He')]:
        expected=st['cold_output'][n]*st['attempts_per_pulse']*d['atomic_mass_kg'][n]
        raw=mass*y[:,names.index(n)];raw*=expected/sum(raw)
        new['HHe_cell_mass_kg'][n]=raw.tolist();new['cold_'+key+'_atomic_kg']=expected
    # The calendar is exported with its exact focus times, including delays.
    with (P/'publication-events.csv').open() as f:
        d['events']=[dict(cat=e['stage'],focus_s=float(e['focus_s'])) for e in csv.DictReader(f)]
    selected=json.loads((PAGE/'selected-plant.json').read_text())
    periods={s['catalyst']:s['average_interval_s'] for s in selected['stages']}
    ns['qmean']=np.array([sum(s['cold_'+key+'_atomic_kg']/periods[n] for n,s in d['stages'].items())
                         for key in ['H','He']])*ns['KB']*ns['T']/(ns['A']*ns['m'])
    ns['Se']*=16
    joint=json.loads((F/'gas/joint-ports.json').read_text())
    joint=next(r for r in joint['scenarios'] if r['downstream_speed_target_m3_s_per_wall_m2_each_species']==.001
               and r['shared_mouth_area_fraction']==.0001 and r['net_molecular_transmission']==.1)
    areas=np.array(joint['required_active_panel_area_fraction_Hgetter_Hecryo'])*16
    panel=areas*np.array([10.,13.])*np.sqrt(300/np.array([300.,273.15]))
    outlet=ns['Sc']*16;ns['Sc']=outlet/(1+outlet*np.sqrt(300/1000)/panel);ns['floor']=np.zeros(2)
    rows=[]
    for dt,points in [(50.,1024),(25.,2048)]:
        _,q,meta=ns['sources'](dt,points);r,*_=ns['response'](q,dt,.5,1e6)
        r['source_discretization']=meta;rows.append(r)
    assert abs(rows[0]['total_buffer_pressure_max_Pa']/rows[1]['total_buffer_pressure_max_Pa']-1)<.001
    out={'scope':__doc__,'cases':rows,'checks':ns['checks'],
         'input_files_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in [src,P/'publication-inputs.json',P/'publication-events.csv',PAGE/'plasma-simulation/data/n15.npz']}}
    (P/'publication-buffer.json').write_text(json.dumps(out,indent=2)+'\n')
    print([(r['total_buffer_pressure_min_Pa'],r['total_buffer_pressure_max_Pa']) for r in rows])


if __name__=='__main__':main()
