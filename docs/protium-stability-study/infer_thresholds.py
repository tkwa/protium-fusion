"""Summarize deterministic brackets; this script does not assign credible intervals."""
import json
from pathlib import Path
from packet_geometry import meridional

root=Path(__file__).resolve().parent
base=Path('/var/tmp/tkwa-chatgpt/protium-stability/runs')

def read(name):
    p=base/name/'diagnostics.json'
    if not p.exists():p=root/'results/cases'/(name+'.json')
    if not p.exists():return None
    d=json.loads(p.read_text());e=d.get('execution') or {}
    return d if e.get('reached_time_limit') and not e.get('fatal_error') else None

def energy(d):return d.get('nuclear_output_J',d['gain']*d['initial_energy_J'])

def curve(control,points,kind):
    c=read(control)
    if not c:return None
    rows=[]
    for x,name in points:
        d=read(name)
        if not d:continue
        rows.append(dict(x=x,name=name,retained_output=energy(d)/energy(c),
                         late_output_fraction=d.get('output_fraction_produced_in_last_2us'),
                         initial_mass_ratio=d['initial_mass_kg']/c['initial_mass_kg'],
                         initial_energy_ratio=d['initial_energy_J']/c['initial_energy_J'],
                         max_energy_residual=d['maximum_energy_residual_fraction_input']))
    rows.sort(key=lambda r:r['x'])
    out=dict(control=control,kind=kind,points=rows,threshold=None)
    for a,b in zip(rows[:-1],rows[1:]):
        ya,yb=a['retained_output'],b['retained_output']
        crossing=(ya<.9<=yb) if kind=='packet_count' else (ya>=.9>yb)
        if not crossing:continue
        transform=(lambda x:1/x) if kind=='packet_count' else (lambda x:x)
        inverse=(lambda x:1/x) if kind=='packet_count' else (lambda x:x)
        frac=(.9-ya)/(yb-ya)
        estimate=inverse(transform(a['x'])+frac*(transform(b['x'])-transform(a['x'])))
        out['threshold']=dict(bracket=[a['x'],b['x']],interpolated=estimate,
                              interpolation='linear in 1/N' if kind=='packet_count' else 'linear between bracket endpoints')
        break
    return out

curves={}
for ell in [1,2,8,16]:
    nt=128 if ell==16 else 64
    names=[]
    for sigma in [0,1,2,3,4,5,5.5,6]:
        name=f'A-384x{nt}-ref-cos{ell}-s{sigma}us'+('' if ell==16 else '-hlle')
        if sigma==0:name='A-384x64-ref-cos2-s0us-hlle'
        names.append((sigma,name))
    curves[f'timing_2d_order{ell}']=curve('A-384x64-ref-cos2-s0us-hlle',names,'timing_us')
curves['packets_2d']=curve('A-384x128-smooth-hllc',[(meridional(k)['circular_cap_equivalent_N'],f'A-384x128-packets-k{k}') for k in [4,5,6,8,16,32]],'packet_count')
for grid in [32,48,64]:
    curves[f'packets_3d_{grid}']=curve(f'A-cart{grid}-smooth',[(n,f'A-cart{grid}-packets-n{n}') for n in [32,64,96,128,256]],'packet_count')
    curves[f'timing_3d_{grid}']=curve(f'A-cart{grid}-wide-ref0-gaussian-s0us',[(s,f'A-cart{grid}-wide-ref0-gaussian-s{s}us') for s in [0,1,2,3,4,4.5,5,6]],'timing_us')
curves['timing_2d_earlier_order8']=curve('A-384x64-ref0-cos8-s0us',[(s,f'A-384x64-ref0-cos8-s{s}us') for s in [0,3,4,4.5,5,5.5,6]],'timing_us')
for mode in ['preheat','preheatfloor']:
    points=[(0,'A-768-1d-preheat-f0')]+[(f,f'A-768-1d-{mode}-f{f}') for f in [.01,.02,.025,.03,.05,.075,.1,.3,.6]]
    curves[mode]=curve('A-768-1d-preheat-f0',points,'kinetic_energy_fraction')
curves={k:v for k,v in curves.items() if v is not None}
(root/'results/deterministic-thresholds.json').write_text(json.dumps(curves,indent=2)+'\n')
for name,c in curves.items():print(name,c['threshold'])
