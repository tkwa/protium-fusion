"""Small-radius diagnostic of the old prescribed hot inflow, not a driver."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1];sys.path.insert(0,str(ROOT))
from radial_pulse import simulate

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--zones',type=int,default=64)
    ap.add_argument('--radius',type=float,default=30.)
    args=ap.parse_args()
    initial=json.loads((ROOT/f'compression-A-{args.zones}-initial.json').read_text())
    initial['outer_radius_m']=(np.array(initial['outer_radius_m'])*args.radius/300).tolist()
    initial['description']='Spatially scaled prescribed hot A inflow; density, velocity and temperature unchanged. No source, pellet or driver model.'
    stem=f'scaled-hot-A-{args.zones}-R{args.radius:g}'
    (HERE/(stem+'-initial.json')).write_text(json.dumps(initial,indent=2)+'\n')
    result=simulate('A',3e8,1e8,args.radius,zones=args.zones,nuclear=True,
                    reduced_a=True,radiation=True,initial_profile=initial,final_time_units=4,
                    save_trace=HERE/(stem+'.npz'),rtol=1e-6,atol=1e-11)
    result['diagnostic_scope']='Spatially scaled prescribed hot inflow; no delivery or multidimensional instability validation.'
    (HERE/(stem+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['mass_kg','initial_energy_J','nuclear_energy_J',
        'deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K',
        'energy_residual_fraction_initial','termination_reason','wall_seconds']},indent=2))

if __name__=='__main__':main()
