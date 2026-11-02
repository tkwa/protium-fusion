#!/usr/bin/env python3
import json, math
from pathlib import Path
import calculate as c
checks=[]
def check(name,test):
 checks.append({'name':name,'passed':bool(test)})
 if not test:raise AssertionError(name)
for q in [.001,.01,.1,.95,1.]:
 for n in [.01,.1,.5,.9]:
  for fun in [c.G_muto,c.G_hasebe]:
   check(f'inverse {fun.__name__} q{q} n{n}',abs(c.inverse_survival(fun(n,q),q,fun)-n)<1e-12)
for f in [.1,.5,.9,.99]:
 check(f'erfc inverse {f}',abs(math.erfc(c.erfc_inverse(f))-f)<1e-14)
 for L in [1e-7,1e-5,1e-4]:
  D=c.required_diffusivity(L,100,f)
  check(f'diffusion threshold {f} {L}',abs(c.diffusion_release(D*100,L)-f)<1e-14)
for T in [723.15,1000,1700]:
 out=c.integrate_clock([0,100,7698.198541942018],[T]*3,8.2,1)
 check(f'isothermal clock {T}',abs(out[-1]/(c.arrhenius(8.2,1,T)*7698.198541942018)-1)<2e-15)
# Independent256/1024-panel Simpson for largest prescribed temperature range.
for E in [.38,.9,1.,3.48,6.93]:
 clock=c.integrate_clock([0,100],[1000,1700],1,E)[-1]
 n=1024;vals=[c.arrhenius(1,E,1000+700*i/n) for i in range(n+1)]
 finer=100/n/3*(vals[0]+vals[-1]+4*sum(vals[1:-1:2])+2*sum(vals[2:-1:2]))
 check(f'variableT convergence E{E}',abs(clock/finer-1)<1e-5)
config={'times_s':[0,1,100,1000,7698.198541942018], 'temperature_K':[1000,1000,1700,1000,1000],
 'initial_H_atoms_m2':1e19,'populations':[{'name':'prompt','model':'prompt','fraction':.1},
 {'name':'held','model':'held','fraction':.1},
 {'name':'mobile','model':'diffusion','fraction':.3,'D0_m2_s':1.2e-4,'E_eV':.9,'distance_m':100e-6},
 {'name':'trap','model':'first_order','fraction':.2,'prefactor_s':1e13,'E_eV':3.},
 {'name':'Muto','model':'muto','fraction':.3,'q':.01}],
 'cavity_fraction':.25,'atomic_H_fraction':.1,'methane_H_fraction':.2}
d=c.response(config)
for i,t in enumerate(config['times_s']):
 released=d['released_H_atoms_m2'][i]; held=d['retained_H_atoms_m2'][i]
 check(f'H closure t{t}',abs((released+held)/1e19-1)<1e-14)
 check(f'route closure t{t}',abs((d['released_cavity_H_atoms_m2'][i]+d['released_collector_H_atoms_m2'][i]-released)/1e19)<1e-14)
 gasH=d['atomic_H_particles_m2'][i]+2*d['H2_molecules_m2'][i]+4*d['CH4_molecules_and_C_atoms_m2'][i]
 check(f'chemistry H closure t{t}',abs((gasH-released)/1e19)<1e-14)
check('monotonic release',all(a<=b for a,b in zip(d['released_H_atoms_m2'],d['released_H_atoms_m2'][1:])))
check('held population survives',d['retained_H_atoms_m2'][-1]>=1e18)
snap=json.loads((c.ROOT/'input-snapshot.json').read_text());res=c.tables(snap)
check('H snapshot baryon mass',all(abs(x['H_kg_m2']/(x['H_atoms_m2']*c.MU)-1)<1e-12 for x in snap['stages'].values()))
check('nucleus-based work',abs(res['power_W_m2_per_eV_per_H']/(res['mean_H_atoms_m2_s']*c.EV)-1)<1e-14)
check('repair linear depth',abs(res['broad_layer_refurbishment_sensitivity'][2]['W_m2_per_MJ_kg_input']/res['broad_layer_refurbishment_sensitivity'][0]['W_m2_per_MJ_kg_input']-100)<1e-12)
(c.ROOT/'validation.json').write_text(json.dumps({'scope':'numerical identities and bookkeeping; not experimental validation of reactor coefficients','all_passed':True,'check_count':len(checks),'checks':checks},indent=2)+'\n')
print(len(checks),'checks passed')
