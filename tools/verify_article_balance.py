#!/usr/bin/env python3
"""Check the published warm-buffer figures against their archived plant account."""
import hashlib
import json
from math import isclose
from pathlib import Path

root = Path(__file__).resolve().parents[1]
account_path = root / 'docs/protium-fusion/component-design/catcher-followup/plant-account/results.json'
interface_path = root / 'public/blog/ordinary-hydrogen-fusion/catcher-interface.json'
account = json.loads(account_path.read_text())
interface = json.loads(interface_path.read_text())
assert hashlib.sha256(account_path.read_bytes()).hexdigest() == interface['source_sha256']['plant-account/results.json']
case = next(row for row in interface['case_rows'] if row['label'] == 'warm_compression_buffer_reference')
source = next(row for row in account['rows'] if row['label'] == case['label'])
assert all(account['checks'].values())
for key in ('net_W', 'mass_kg'):
    assert isclose(case[key], source['plant'][key], rel_tol=1e-12), key
for key, value in case['energy_balance'].items():
    if isinstance(value, (float, int)):
        assert isclose(value, source['plant'][key], rel_tol=1e-12), key
energy = case['energy_balance']
assert isclose(energy['nuclear_total_W'], energy['neutrino_W'] + energy['cold_rejection_W'] + case['net_W'], rel_tol=1e-12)
assert isclose(case['net_W'] / case['mass_kg'] * 1000, case['net_mW_kg'], rel_tol=1e-12)
print(f"Warm-buffer source hash and energy balance verified: {case['net_mW_kg']:.4f} mW/kg")
