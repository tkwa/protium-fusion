"""Recompute numerical comparison and terminal kinematics from saved results."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument('--raw',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args()
names=[
 'B2-64-compression-cap', 'B2-64-compression-cap-tight',
 'B2-64-compression-cap-t8', 'B2-64-divergence-sensitivity',
 'B2-128-compression-cap', 'B2-128-compression-cap-energy-matched',
 'B2-256-compression-cap-energy-matched',
]
rows={}
base=json.loads((a.raw/(names[0]+'.json')).read_text())
for name in names:
 p=a.raw/(name+'.json')
 if not p.exists():continue
 r=json.loads(p.read_text());v=np.asarray(r['profile']['velocity_m_s'])
 m=np.asarray(r['profile']['mass_fraction_of_total'])
 node_mass=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
 row={k:r[k] for k in ['zones','success','termination_reason','elapsed_sound_times',
  'deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst',
  'peak_temperature_K','energy_residual_fraction_initial','initial_energy_J',
  'final_internal_J','final_kinetic_J','radiation_escaped_J','wall_seconds']}
 row['result_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
 row['relative_gain_change_from_64']=r['deposited_nuclear_over_initial_energy']/base['deposited_nuclear_over_initial_energy']-1
 row['relative_yield_change_from_64']=r['desired_product_per_initial_catalyst']/base['desired_product_per_initial_catalyst']-1
 row['final_internal_fraction_of_internal_plus_kinetic']=r['final_internal_J']/(r['final_internal_J']+r['final_kinetic_J'])
 row['terminal_outward_and_strictly_monotone']=bool(np.all(v>0) and np.all(np.diff(v)>0))
 row['minimum_adjacent_velocity_increase_m_s']=float(np.min(np.diff(v)))
 row['maximum_velocity_m_s']=float(np.max(v))
 row['outer_radius_m']=r['profile']['outer_radius_m'][-1]
 row['minimum_final_temperature_K']=min(r['profile']['T_K'])
 row['regridding']=r['regridding']
 if name.endswith('-t8'):
  v0=np.asarray(base['profile']['velocity_m_s']);dv=v/v0-1
  row['velocity_growth_from_four_times_maximum']=float(np.max(dv))
  row['velocity_growth_from_four_times_node_mass_weighted']=float(node_mass@dv/node_mass.sum())
  row['target_absolute_change_from_four_times']=r['desired_product_per_initial_catalyst']-base['desired_product_per_initial_catalyst']
 rows[name]=row
out={'scope':'Completed specified-hot-inflow spherical model; numerical errors and model approximations remain distinct.','cases':rows}
if names[-2] in rows:
 fine=json.loads((a.raw/(names[-2]+'.json')).read_text())
 x=np.asarray(base['profile']['desired_product_per_initial_catalyst'])
 y=np.asarray(fine['profile']['desired_product_per_initial_catalyst'])
 m=np.asarray(base['profile']['mass_fraction_of_total'])
 n=np.asarray(fine['profile']['mass_fraction_of_total'])
 rebinned=(y*n).reshape(64,2).sum(axis=1)/m
 contribution=m*(rebinned-x)
 out['matched_128_yield_change_by_original_cell']={
  'weighted_target_contribution_by_parent_cell':contribution.tolist(),
  'all_parent_bins_increase':bool(np.all(contribution>=0)),
  'outer_eight_parent_cells_fraction_of_gain':float(contribution[-8:].sum()/contribution.sum()),
  'outer_eight_parent_cells_mass_fraction':float(m[-8:].sum()),
 }
if names[-1] in rows and names[-2] in rows:
 lo=rows[names[-2]];hi=rows[names[-1]]
 out['matched_128_to_256']={
  'relative_gain_change':hi['deposited_nuclear_over_initial_energy']/lo['deposited_nuclear_over_initial_energy']-1,
  'relative_yield_change':hi['desired_product_per_initial_catalyst']/lo['desired_product_per_initial_catalyst']-1,
 }
a.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
