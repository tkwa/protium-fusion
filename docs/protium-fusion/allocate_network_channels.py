"""Finite family of equilibria, with global area/energy optimization.

For fixed processor flows, choose one configuration per isotope to maximize
(Ebare - sum flow*loss)/(sum flow*area). Dinkelbach iteration selects the minimum
loss + ratio*area independently for each isotope, then updates the ratio.
This optimizes the saved finite family, not all physical plasma states.
"""
from pathlib import Path
import json,gzip,argparse,time
from math import pi
import numpy as np
from scipy.optimize import minimize_scalar
from catalyst_network import load_data,HOT,label
from catalyst_catalog import nucmass
from catalyst_penning import equilibrium
from bundle_solenoid import source_equilibrium
from network_ion_budget import nuclear_energy,geometry
ROOT=Path(__file__).parent

def family(energy,B=140):
    allrows={};stats={};t0=time.monotonic();T=energy['T_MeV'];ref=source_equilibrium()
    for s in energy['stages']:
        za=tuple(s['ZA']);z,a=za;flow=s['pn_per_primary']
        if flow<=0:continue
        scale=z*nucmass(1,1)/(2*nucmass(z,a));rows=[];tried=0;valid=0
        params=[(xi,scale*frac*share,scale*frac*(1-share)/z)
          for xi in [.3,.4,.5] for frac in [.65,.75,.85,.9,.95,1.,1.05,1.1]
          for share in [.65,.75,.8,.85,.9,.95]]
        if za==(22,49):params.append((ref['xi'],ref['h0'],ref['c0']))
        for xi,h,c in params:
            tried+=1;eq=equilibrium(z,a,h,c,xi,nout=6001,rtol=1e-10)
            if not eq['success']:continue
            valid+=1
            for x in [18.,20.,22.,24.,26.,28.,30.,34.]:
                rows.append(geometry(eq,x,s,T=T,B=B))
        assert rows,za
        allrows[label(za)]=rows;stats[label(za)]={'attempted':tried,'closed_equilibria':valid,'wall_configurations':len(rows)}
        print(label(za),stats[label(za)],'seconds',round(time.monotonic()-t0,2),flush=True)
    return {'scope':__doc__,'T_MeV':T,'B_T':B,'statistics':stats,'configurations':allrows}


def choose(energy,candidates):
    ratio=0.;history=[];groups=candidates['configurations'];flows=energy['network']['processor_pn_rates']
    for _ in range(100):
        selected={za:min(rows,key=lambda r:r['total_loss_MeV_per_pn']+ratio*r['area_per_pn_rate_m3_s']) for za,rows in groups.items()}
        loss=sum(flows[za]*r['total_loss_MeV_per_pn'] for za,r in selected.items())
        area=sum(flows[za]*r['area_per_pn_rate_m3_s'] for za,r in selected.items())
        net=energy['bare_electricity_MeV_per_primary']-loss;new=net/area
        history.append({'ratio_MeV_m3_s':new,'net_MeV_per_primary':net,'area_m3_s':area,'fixed_point_residual_MeV':net-ratio*area})
        if abs(new-ratio)<=1e-12*max(abs(new),1):break
        ratio=new
    else:raise RuntimeError('fractional optimization did not settle')
    return {'scope':__doc__,'selected':selected,'history':history,'net_source_MeV_per_primary':net,
            'channel_volume_per_primary_rate_m3_s':area,'source_power_W_m3_at_full_channel_fill':net*1.602176634e-13/area,
            'area_fractions':{za:flows[za]*r['area_per_pn_rate_m3_s']/area for za,r in selected.items()}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reuse-family',action='store_true');args=p.parse_args()
    energy=json.loads((ROOT/'network-ion-budget.json').read_text())
    path=ROOT/'network-channel-family.json.gz'
    if args.reuse_family:fam=json.loads(gzip.decompress(path.read_bytes()))
    else:
        fam=family(energy);path.write_bytes(gzip.compress(json.dumps(fam,separators=(',',':'),allow_nan=False).encode(),mtime=0))
    result=choose(energy,fam);(ROOT/'network-channel-allocation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print({k:v for k,v in result.items() if k!='selected'})
