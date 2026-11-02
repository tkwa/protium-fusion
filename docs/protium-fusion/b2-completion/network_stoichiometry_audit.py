"""Audit reaction metadata independently of the runner's stoichiometry arrays."""
import sys,json
from pathlib import Path
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE))
import pulse_model as p
rows=[]
for j,r in enumerate(p.RATES):
    db=sum(n.A for n in r.products)-sum(n.A for n in r.reactants)
    dz=sum(n.Z for n in r.products)-sum(n.Z for n in r.reactants)
    q=sum(n.dm for n in r.reactants)-sum(n.dm for n in r.products)
    rows.append({'rate':str(r),'weak_type':r.weak_type,'baryon_delta':db,'charge_delta':dz,'mass_excess_Q_MeV':q,'runner_Q_MeV':float(p.Q[j]),'reaclib_Q_MeV':r.Q,'reaclib_minus_mass_excess_Q_MeV':r.Q-q})
out={'rates':len(rows),'weak_rates':sum(bool(r['weak_type']) for r in rows),
 'max_baryon_delta':max(abs(r['baryon_delta']) for r in rows),
 'max_strong_charge_delta':max(abs(r['charge_delta']) for r in rows if not r['weak_type']),
 'max_Q_mass_excess_minus_runner_MeV':max(abs(r['mass_excess_Q_MeV']-r['runner_Q_MeV']) for r in rows),
 'max_reaclib_minus_mass_excess_Q_MeV':max(abs(r['reaclib_minus_mass_excess_Q_MeV']) for r in rows),'rates_metadata':rows,
 'limits':['Shared isotope library and rate metadata; validates accounting rather than measured cross sections.',
 'Small Reaclib Q versus modern mass-table differences do not by themselves quantify reverse-rate detailed-balance error.']}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rates_metadata'},indent=2))
