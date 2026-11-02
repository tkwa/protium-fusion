"""Conditional small-stage material/energy ledger using existing accounting code."""
import argparse,hashlib,json,sys
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--factory',default='compression-He-64-r30-burning.json');ap.add_argument('--name',default='cycle-check.json');a=ap.parse_args();sys.path.insert(0,str(a.root));import cycle_ledger as cl
base={'n14':('scaled-hot-A-64-R30.json',1),'n15':('scaled-hot-B1-64-R30.json',1),'c13':('scaled-hot-C13-64-R30.json',1)}
variants=[]
for b2,scope in [('scaled-hot-B2-64-R30-tight.json','Conditional freeze at failed radiating solver endpoint; not a completed radiating pulse.'),('scaled-hot-B2-64-R30-norad.json','Radiation-off endpoint control; read B2_success and B2_termination. Not a completed radiating solution.')]:
 p=a.output/b2
 if not p.exists():continue
 d=json.loads(p.read_text())
 files={k:(str(a.output/v),n) for k,(v,n) in base.items()};files['c12']=(str(p),2)
 cases=[cl.ledger(collection_loss=loss,files=files,factory_file=a.factory) for loss in [0.,1e-5,1e-4,1e-3]]
 for case in cases:
  for stage in case['stage_data']:stage['source']=Path(stage['source']).name
 variants.append({'B2_input':b2,'B2_success':d['success'],'B2_termination':d.get('effective_termination_reason',d['termination_reason']),'scope':scope,'cases':cases})
 print(b2,[(x['collection_loss_per_stage'],x['available_heat_MeV_per_helium'],x['fuel_preparation_MeV_per_helium'],x['net_energy_cases'][2]['net_MeV_per_helium_before_other_auxiliaries']) for x in cases],flush=True)
out={'scope':'Existing steady-state isotope accounting applied to small prescribed hot inflows with ideal separation and the named helium factory; B2 variants are explicitly conditional. No engineered driver or collection demonstrated.','factory_source':Path(a.factory).name,'factory_sha256':hashlib.sha256((a.root/a.factory).read_bytes()).hexdigest(),'ledger_code_sha256':hashlib.sha256((a.root/'cycle_ledger.py').read_bytes()).hexdigest(),'input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in a.output.glob('scaled-hot-*-64-R30*.json') if 'initial' not in p.name and 'inputs' not in p.name and 'failure' not in p.name},'variants':variants}
(a.output/a.name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
