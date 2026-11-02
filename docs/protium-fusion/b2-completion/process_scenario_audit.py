"""Independent checks of mesh families and net structural-carbon reserves."""
import sys,json,hashlib
from pathlib import Path
from independent_cycle_audit import BASE,STAGES,load_records,open_external_carbon
from open_process_comparison import FIELDS

def run(path):
 data=json.loads(Path(path).read_text());rows=[]
 for family,f in data['families'].items():
  records=load_records({stage:str(BASE/f['stage_map'][cat]) for cat,(stage,p) in STAGES.items()})
  for label,p in f['cases'].items():
   ours=open_external_carbon(records,1-p['carbon_loss_per_cold_stage_output'],1-p['hydrogen_loss_per_cold_stage_output'],p['C13_number_fraction_in_fresh_carbon'])
   rel={k:(ours[k]-p[v])/max(abs(ours[k]),abs(p[v]),1e-300) for k,v in FIELDS.items() if ours[k] is not None}
   total=p['total_ongoing_C_H_mass_budget'];target=p['fuel_C_H_target']
   reserve=(total-target)*ours['fresh_hydrogen_kg_per_helium'];perJ=reserve/(ours['selected_60driver_40heat_99capture_net_MeV_per_helium']*1.602176634e-13)
   checks={'independent_process_fields':max(abs(x) for x in rel.values())<1e-8,
    'fuel_C_H_target':abs(ours['fresh_C_to_H_mass_ratio']-target)<1e-10,
    'reserve_same_H_denominator':abs(reserve-p['net_structural_carbon_allowance_kg_per_He'])<1e-8*max(reserve,1e-300),
    'net_reserve_per_J':abs(perJ-p['net_structural_carbon_allowance_kg_per_preaux_J'])<1e-8*max(perJ,1e-300),
    'total_C_budget':abs((ours['fresh_carbon_kg_per_helium']+reserve)/ours['fresh_hydrogen_kg_per_helium']-total)<1e-10,
    'nuclear_external_Q':abs(ours['mass_excess_Q_sum_minus_net_MeV'])<1e-8,
    'baryons':abs(ours['baryon_residual'])<1e-9}
   rows.append({'family':family,'fuel_target':label,'relative_process_differences':rel,'checks':{k:bool(v) for k,v in checks.items()},'independent':ours})
 return {'scope':'Six-variable independent atom and energy balance, including net material-carbon reserve with the same fresh-H denominator. No material recovery qualification.',
  'source':str(path),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'cases':rows,
  'checks_passed':sum(sum(x['checks'].values()) for x in rows),'checks_failed':sum(sum(not c for c in x['checks'].values()) for x in rows)}
if __name__=='__main__':
 r=run(sys.argv[1]);Path(sys.argv[2]).write_text(json.dumps(r,indent=2)+'\n');print({k:r[k] for k in ['checks_passed','checks_failed']})
