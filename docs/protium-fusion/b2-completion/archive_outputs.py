"""Archive exact run metadata and compact, unmodified samples of full traces.

All calculations/audits use the full raw trace. The compact trace is a storage
artifact: no interpolation, rounding, or state renormalization is performed.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

ap = argparse.ArgumentParser()
ap.add_argument('--raw', type=Path, required=True)
ap.add_argument('--destination', type=Path, required=True)
ap.add_argument('--names', nargs='+', required=True)
args = ap.parse_args()
args.destination.mkdir(parents=True, exist_ok=True)
entries = []
for name in args.names:
    result_path = args.raw / f'{name}.json'
    result = json.loads(result_path.read_text())
    assert result['model'] == 'radial Lagrangian pulse expansion'
    for suffix in ('.json', '.initial.json', '.inputs.json'):
        source = args.raw / (name + suffix)
        if source.exists():
            destination = args.destination / source.name
            shutil.copyfile(source, destination)
            assert sha(source) == sha(destination)
    source = args.raw / (name + '.npz')
    with np.load(source) as data:
        temp = data['state'][len(data['species'])::int(data['variables_per_zone'])]
        samples = data['state'].shape[1]
        selected = np.unique(np.r_[
            np.rint(np.linspace(0, samples - 1, min(48, samples))).astype(int),
            np.unravel_index(np.argmax(temp), temp.shape)[1],
            np.unravel_index(np.argmin(temp), temp.shape)[1],
        ])
        content = {key: data[key] for key in data.files}
        content['state'] = content['state'][:, selected]
        content['time_s'] = content['time_s'][selected]
        content['full_trace_sample_indices'] = selected
        dest = args.destination / (name + '.trace.npz')
        np.savez_compressed(dest, **content)
        with np.load(dest) as check:
            assert np.array_equal(check['state'], data['state'][:, selected])
            assert np.array_equal(check['time_s'], data['time_s'][selected])
            assert selected[0] == 0 and selected[-1] == samples - 1
    entry = {
        'result': result_path.name,
        'result_sha256': sha(result_path),
        'full_trace_path': str(source),
        'full_trace_sha256': sha(source),
        'compact_trace': dest.name,
        'compact_trace_sha256': sha(dest),
        'full_sample_count': samples,
        'compact_sample_count': len(selected),
        'selection': '48 evenly spaced saved-state indices plus saved global temperature minimum/maximum; exact original first and final states',
        'full_trace_retained': True,
    }
    (args.destination / (name + '-trace-reference.json')).write_text(json.dumps(entry, indent=2) + '\n')
    entries.append(entry)
print(json.dumps(entries, indent=2))
