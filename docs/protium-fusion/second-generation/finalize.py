#!/usr/bin/env python3
"""Write a compact derived summary and a noncircular artifact/provenance manifest."""
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()

def main():
    data=json.loads((P/'plant-results.json').read_text());c=next(r for r in data['cases'] if r['name']=='spherical_large_N15_16')
    oxygen=json.loads((P/'oxygen-results.json').read_text());pen=oxygen['large_N15_16x_source_comparison']['allocated_incremental_net_W_bounds']
    summary=dict(status='Conditional engineering candidate; 16x specific output remains unqualified',
        selected_reference_W_kg=data['reference'],candidate_128_W_kg=c['net_W_kg'],
        candidate_128_power_W=c['net_W'],candidate_128_mass_kg=c['mass_kg'],
        allocated_O_return_net_W_bounds=[c['net_W']+v for v in pen],
        allocated_O_return_W_kg_bounds=[(c['net_W']+v)/c['mass_kg'] for v in pen],
        bounds_meaning='Deterministic weak-neutrino-energy bracket, with fixed apparatus. Not a design credible interval.',
        oxygen_combination_scope='First-order allocated oxygen return at fixed primary quartet rate and hardware. Full coupled eight-state material/dispatch equilibrium and oxygen-only or mixed-source apparatus remain unqualified. The exact separate oxygen reaction ledger gives its additional protons, alpha products and nuclear energy.',
        blocks=['1600K receiver spectral transients, carbon emission and material service life',
            'Cold-packet formation, reflective capsule/neutralizer survival and launcher preparation rate',
            'Completed cold-intercept, refrigerator and short-territory equipment bill',
            'Current O16/O17/F19 rate replacement and coupled mixed-pulse network/hydro'],
        canon='Article and selected operating point preserved')
    (P/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Preserve exact numerical source versions used by completed runs.
    for f in P.glob('*.py'):
        if f.name in ['network.py','pulse_probe.py']:(P/'source-versions'/(sha(f)+'.py')).write_bytes(f.read_bytes())
    files={}
    for f in sorted(P.rglob('*')):
        if not f.is_file() or '__pycache__' in f.parts or f.name in ['manifest.json','validation.json','equipment-run.txt']:continue
        files[str(f.relative_to(P))]=dict(bytes=f.stat().st_size,sha256=sha(f))
    traces={}
    for f in P.glob('pulse-n15-*.json'):
        r=json.loads(f.read_text());tr=Path(r['trace_path'])
        traces[str(tr)]=dict(bytes=tr.stat().st_size,sha256=sha(tr))
    versions={}
    for package in ['numpy','scipy','pynucastro','numba','llvmlite','matplotlib']:
        try:versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:versions[package]='not on this process import path'
    rev=lambda ref:subprocess.check_output(['git','rev-parse',ref],cwd=ROOT,text=True).strip()
    out=dict(run_time_UTC=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        runtime_packages=versions,repository_head=rev('HEAD'),refreshed_origin_main=rev('origin/main'),
        files=files,persistent_raw_traces=traces,
        exclusions='validation.json regenerates checks and is excluded to avoid circular hashes; Python caches and the duplicate equipment stdout log are excluded.',
        scope='New research packet only. No article or original selected point changed. The primary oxygen/F19 data review is newer than the frozen numerical rate set.',
        reproduced_controls='Archived N15 conversion, exact baseline dispatch, abundance/Q replays, mass/charge/nuclear/electric/radiator conservation.',
        qualifications=summary['blocks'])
    (P/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(artifacts=len(files),raw_traces=len(traces),run_time_UTC=out['run_time_UTC'],summary=summary),indent=2))

if __name__=='__main__':main()
