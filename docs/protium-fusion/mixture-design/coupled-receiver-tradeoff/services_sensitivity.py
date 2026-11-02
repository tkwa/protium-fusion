"""Explicit service-allocation scenarios, with the complete plant re-solved.

The 400/700/1000/1500 kg/m² choices are not a probability distribution or a
completed component bill of materials. Source and other named masses remain.
"""
from pathlib import Path
import argparse,hashlib,json,sys,types
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import receiver_tradeoff as rt
import plant_compare as pc
import station_generalized as generalized

def main():
    a=argparse.ArgumentParser();a.add_argument('--family',type=Path,required=True);a.add_argument('--tradeoff',type=Path,required=True);a.add_argument('--sigma',type=float,default=250);a.add_argument('--out',type=Path,required=True);z=a.parse_args()
    family=json.loads(z.family.read_text());trade=json.loads(z.tradeoff.read_text())
    base=next(r for r in trade['rows'] if r['receiver_sigma_kg_m2']==z.sigma and r['prompt_policy']=='EM20' and r['missed_hot_routing']=='cold_equipment' and r['additional_missed_hot_fraction_HplusUprecursor']==0)
    selected=family['plant_selected'];cycle=family['precursor_cycle'];rows=[]
    for service_sigma in [400.,700.,1000.,1500.]:
        def evaluate(baseline,cy,**kwargs):
            b=dict(baseline);b['mass_breakdown_kg']=dict(baseline['mass_breakdown_kg'])
            b['mass_breakdown_kg']['processing_structure_thermal_service']=service_sigma*4*generalized.math.pi*b['radius_m']**2
            result=generalized.evaluate(b,cy,**kwargs)
            result['parameters']['service_areal_allocation_kg_m2']=service_sigma
            return result
        pc.cold=types.SimpleNamespace(MEV=generalized.MEV,evaluate=evaluate)
        p=rt.optimized(selected,cycle,z.sigma,1.,base['total_hot_fraction_HplusUprecursor'])
        assert p['feasible'],p.get('reasons')
        assert max(abs(x) for x in p['closure_relative_residuals'].values())<1e-10
        if service_sigma==1000:
            assert abs(p['net_W_kg']/base['plant']['net_W_kg']-1)<1e-12
            assert abs(p['mass_kg']/base['plant']['mass_kg']-1)<1e-12
        rows.append({'service_kg_m2':service_sigma,'mW_kg':1000*p['net_W_kg'],'net_W':p['net_W'],'mass_kg':p['mass_kg'],'plant':p})
    out={'scope':__doc__,'receiver_sigma_kg_m2':z.sigma,'family':family.get('family'),'rows':rows,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [z.family,z.tradeoff,Path(__file__),Path(pc.__file__),Path(rt.__file__),Path(generalized.__file__)]},'interpretation':'The 1000 kg/m² row reproduces the main account. Smaller allowances require a completed recovery, valve, actuator, joint and service equipment design fitting within the reduced budget. Larger allowances illustrate mass sensitivity. No loss of retained carbon is implied.'}
    z.out.parent.mkdir(parents=True,exist_ok=True);z.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k!='plant'} for r in rows],indent=2))
if __name__=='__main__':main()
