#!/usr/bin/env python3
"""Verify the current publication against its frozen research input snapshot."""
import hashlib
import json
from math import isclose
from pathlib import Path

root=Path(__file__).resolve().parents[1]
base=root/'docs/protium-fusion/second-generation'
snapshot_path=base/'publication-inputs.json'
snapshot=json.loads(snapshot_path.read_text())
page=root/'public/blog/ordinary-hydrogen-fusion'
interface=json.loads((page/'catcher-interface.json').read_text())
selected=json.loads((page/'selected-plant.json').read_text())
digest=hashlib.sha256(snapshot_path.read_bytes()).hexdigest()
assert interface['publication_input_sha256']==selected['publication_input_sha256']==digest
assert selected['release_version']==(root/'VERSION').read_text().strip()
case=interface['case_rows'][0];source=snapshot['plant'];energy=case['energy_balance']
for key in ('net_W','mass_kg'):
    assert isclose(case[key],source[key],rel_tol=1e-12),key
    assert isclose(selected['plant'][key],source[key],rel_tol=1e-12),key
assert isclose(case['net_W']/case['mass_kg']*1000,case['net_mW_kg'],rel_tol=1e-12)
assert isclose(energy['nuclear_total_W'],energy['neutrino_W']+energy['cold_rejection_W']+case['net_W'],rel_tol=1e-12)
assert isclose(energy['captured_pulse_W'],energy['nuclear_non_neutrino_W']+energy['prepared_pulse_W'],rel_tol=1e-12)
assert isclose(energy['gross_electric_W']-energy['driver_W']-energy['ordinary_auxiliary_W']-energy['replacement_W']-energy['support_W'],case['net_W'],rel_tol=1e-12)
assert isclose(sum(source['mass_breakdown_kg_m2'].values()),source['installed_kg_m2'],rel_tol=1e-12)
stages={s['catalyst']:s for s in selected['stages']}
assert stages['n15']['cold_outputs_per_initial_catalyst']==snapshot['probe']['stage']['cold_output']
assert isclose(stages['n15']['mass_kg']/stages['n14']['mass_kg'],1,rel_tol=1e-12)
assert selected['cycle']['stage_attempts_per_He']==snapshot['probe']['cycle']['stage_attempts_per_He']
assert source['direct_fraction_of_all_K']==0
assert selected['plant']['radiator_area_fraction']<.9
print(f"v{selected['release_version']} source identity, energy, mass and isotope accounts verified: {case['net_mW_kg']:.4f} mW/kg")
