"""Summarize recorded Athena resource limits, completion, and time-window compliance."""
import datetime as dt,json,subprocess
from pathlib import Path
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
out=Path(__file__).resolve().parent/'results'
deadline=dt.datetime(2026,10,31,8,22,56,tzinfo=dt.timezone.utc)
mapping={}
for row in subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True).splitlines():
    if row.startswith('#'):continue
    cpu,core,socket=map(int,row.split(','));mapping[cpu]=(socket,core)
records=[];used=set()
for p in sorted((base/'runs').iterdir()):
    if not (p/'config.json').exists():continue
    c=json.loads((p/'config.json').read_text())
    cores={int(x) for x in c['cores'].split(',')};used|=cores
    e=json.loads((p/'execution.json').read_text()) if (p/'execution.json').exists() else None
    start=dt.datetime.fromisoformat(c['started_utc'])
    end=start+dt.timedelta(seconds=e['elapsed_wall_s']) if e and 'elapsed_wall_s' in e else None
    records.append(dict(name=p.name,logical_cpus=sorted(cores),physical_cores=len({mapping[x] for x in cores}),
                        started_utc=start.isoformat(),approximate_end_utc=end.isoformat() if end else None,
                        reached_time_limit=bool(e and e.get('reached_time_limit') and not e.get('fatal_error')),
                        incomplete=e is None,ended_after_deadline=bool(end and end>deadline)))
report=dict(deadline_utc=deadline.isoformat(),logical_cpus_used=sorted(used),
            union_physical_cores=len({mapping[x] for x in used}),
            note='Athena runs with recorded launch metadata. All launches also use an independent timeout ending before the deadline; approximate end times add recorded wall duration to configuration time.',
            total_attempts=len(records),completed=sum(r['reached_time_limit'] for r in records),
            still_running=[r['name'] for r in records if r['incomplete']],
            deadline_exceptions=[r['name'] for r in records if r['ended_after_deadline']],runs=records)
assert report['union_physical_cores']<8
(out/'campaign-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='runs'},indent=2))
