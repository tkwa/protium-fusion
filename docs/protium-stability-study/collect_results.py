"""Retain small, auditable research outputs alongside the analysis code."""
import gzip,hashlib,json,shutil
from pathlib import Path
from analyze import analyze
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
out=Path(__file__).resolve().parent/'results';out.mkdir(exist_ok=True)
cases=out/'cases';cases.mkdir(exist_ok=True)
logs=out/'logs';logs.mkdir(exist_ok=True)
inventory=[]
for d in sorted((base/'runs').iterdir()):
    execution=d/'execution.json'
    if not execution.exists():continue
    e=json.loads(execution.read_text());cfg=json.loads((d/'config.json').read_text())
    row=dict(name=d.name,config=cfg,execution=e)
    log=(d/'stdout.log').read_text()
    with gzip.open(logs/(d.name+'.log.gz'),'wb') as f:f.write(log.encode())
    shutil.copy2(d/'athinput',logs/(d.name+'.athinput'))
    completed=e.get('reached_time_limit','Terminating on time limit' in log)
    failed=e.get('fatal_error','FATAL ERROR' in log or e.get('returncode',1)!=0)
    if completed and not failed and cfg['test']==0:
        diag=d/'diagnostics.json'
        if not diag.exists() or diag.stat().st_mtime < max(execution.stat().st_mtime,Path(__file__).with_name('analyze.py').stat().st_mtime):analyze(d)
        shutil.copy2(diag,cases/(d.name+'.json'))
        with gzip.open(cases/(d.name+'.hst.gz'),'wb') as f:
            f.write((d/'pulse.hst').read_bytes())
    elif (d/'diagnostics.json').exists():
        partial=out/'incomplete';partial.mkdir(exist_ok=True)
        shutil.copy2(d/'diagnostics.json',partial/(d.name+'.json'))
    row['input_sha256']=hashlib.sha256((d/'athinput').read_bytes()).hexdigest()
    row['log_tail']=(d/'stdout.log').read_text()[-1500:]
    inventory.append(row)
(out/'run-inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
for name in ['packet-geometry.json','validation/adapter-validation.json','validation/noh-validation.json','validation/limiter-validation.json','validation/phase-roundtrip.json']:
    if (base/name).exists():shutil.copy2(base/name,out/Path(name).name)
for source in (base/'validation').glob('timing-quadrature*.json'):
    shutil.copy2(source,out/source.name)
for name in ['lagrangian-yi2','lagrangian-physical','lagrangian-physical-trapped','lagrangian-fixed-diffusion']:
    target=out/name;target.mkdir(exist_ok=True)
    for f in (base/name).glob('A-*.json'):shutil.copy2(f,target/f.name)
target=out/'inputs';target.mkdir(exist_ok=True)
for name in ['profile.txt','rates.txt','manifest.json']:shutil.copy2(base/'inputs-yi2'/name,target/name)
for source in base.glob('packet-times-*'):
    if source.suffix in ['.txt','.json']:shutil.copy2(source,target/source.name)
for version in (base/'versions').iterdir():
    target=out/'source-versions'/version.name;target.mkdir(parents=True,exist_ok=True)
    for source in version.glob('pulse*'):shutil.copy2(source,target/source.name)
    if (version/'framework.patch').exists():shutil.copy2(version/'framework.patch',target/'framework.patch')
print('Retained',len(inventory),'run records in',out)
