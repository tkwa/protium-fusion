"""Diagnostic charged-particle screening factors at selected pulse states."""
import json
import math
from pathlib import Path
import pynucastro as p
from pynucastro.screening import make_plasma_state,make_screen_factors,chugunov_2009

rows=[]
for cat,np_,rho,temps in [('n14',1,1e8,[3e8,1.45e9]),
                          ('n15',1,1e6,[2e8,4.5e8]),
                          ('c12',2,1e8,[3e8,1.45e9]),
                          ('o16',2,1e10,[5e8,1.6e9])]:
    n=p.Nucleus(cat);comp=p.Composition(['p',cat])
    comp['p']=np_/(n.A+np_);comp[cat]=n.A/(n.A+np_)
    fac=make_screen_factors('p',cat)
    for temp in temps:
        state=make_plasma_state(temp,rho/1000,comp.get_molar())
        # Version 3.1 returns ln(enhancement), not enhancement itself.
        log_factor=float(chugunov_2009(state,fac))
        rows.append({'composition':f'{cat}+{np_}p','rho_kg_m3':rho,'T_K':temp,
                     'log_factor':log_factor,'factor':math.exp(log_factor)})
out={'source':'pynucastro 3.1.0 implementation of Chugunov & DeWitt 2009',
     'scope':'Forward-rate diagnostics at fixed initial density/composition; not a screened pulse integration or Coulomb-corrected EOS.',
     'rows':rows}
Path(__file__).with_name('screening-screen.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
