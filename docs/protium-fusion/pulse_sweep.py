"""Find approximate 90%-conversion pulse radii in the reduced model."""
import json
from pathlib import Path
import numpy as np
import pulse_model as pm
from fast_eos import FastEOS

pm.EOS=FastEOS()
OUT=Path(__file__).with_suffix('.json')
records=[];solutions=[]
def save():
    OUT.write_text(json.dumps({'criterion':'90% desired product per initial catalyst',
      'solutions':solutions,'all_runs':records},indent=2,allow_nan=False)+'\n')
def run(stage,temp,rho,radius):
    r=pm.simulate(stage,temp,rho,radius,rtol=2e-6)
    r['initial_electron_eta']=None
    records.append(r)
    print(stage,temp,rho,f'R={radius:.4g}',f'burn={r["desired_product_per_initial_catalyst"]:.5f}',
      f'Q/u={r["deposited_nuclear_over_initial_EOS_energy"]:.3f}',flush=True)
    save()
    return r

if __name__=='__main__':
    cases=[(s,3e8,rho) for rho in [1e7,1e8,1e9] for s in ['A','B2']]
    cases += [(s,t,1e8) for t in [2e8,5e8] for s in ['A','B2']]
    for stage,temp,rho in cases:
        guess=(100 if stage=='A' else 300)*(1e8/rho)
        low=guess/30;high=guess*30
        rlo=run(stage,temp,rho,low);rhi=run(stage,temp,rho,high)
        for _ in range(4):
            if rlo['desired_product_per_initial_catalyst']<.9: break
            high=low;rhi=rlo;low/=10;rlo=run(stage,temp,rho,low)
        for _ in range(4):
            if rhi['desired_product_per_initial_catalyst']>=.9: break
            low=high;rlo=rhi;high*=10;rhi=run(stage,temp,rho,high)
        if rlo['desired_product_per_initial_catalyst']<.9<=rhi['desired_product_per_initial_catalyst']:
            for _ in range(7):
                mid=np.sqrt(low*high);rmid=run(stage,temp,rho,mid)
                if rmid['desired_product_per_initial_catalyst']>=.9: high=mid;rhi=rmid
                else: low=mid;rlo=rmid
            solutions.append({'stage':stage,'T0_K':temp,'rho0_kg_m3':rho,
              'radius_bracket_m':[low,high],'upper_bracket_result':rhi})
        else:
            solutions.append({'stage':stage,'T0_K':temp,'rho0_kg_m3':rho,'status':'not bracketed'})
        save()
