"""Archive the exact binary tables and profiles used by the Athena calculations."""
import hashlib,json,tarfile
from pathlib import Path
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
out=Path(__file__).resolve().parent/'results';out.mkdir(exist_ok=True)
paths=[]
for directory in ['inputs-yi2','inputs-backward']:
    paths.extend(p for p in (base/directory).iterdir() if p.is_file())
paths.extend(base/name for name in ['phase-reference.bin','phase-signed.bin','phase-signed-wide.bin'])
paths.extend(base.glob('packet-times-*'))
paths=sorted(set(paths))
archive=out/'reproduction-inputs.tar.gz'
with tarfile.open(archive,'w:gz',compresslevel=6) as tar:
    for p in paths:tar.add(p,arcname=str(p.relative_to(base)),recursive=False)
manifest=dict(archive=archive.name,archive_bytes=archive.stat().st_size,
              archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
              files=[dict(path=str(p.relative_to(base)),bytes=p.stat().st_size,
                          sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths])
(out/'reproduction-inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:manifest[k] for k in ['archive','archive_bytes','archive_sha256']}))
