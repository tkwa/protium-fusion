"""Sequential bounded cases on one disjoint CPU set; never outlive authorization."""
import argparse,datetime,json,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('queue',type=Path)
p.add_argument('--cores',required=True)
p.add_argument('--wait-case')
a=p.parse_args()
base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
deadline=datetime.datetime(2026,10,31,8,22,56,tzinfo=datetime.timezone.utc).timestamp()
while a.wait_case and not (base/'runs'/a.wait_case/'execution.json').exists():
    if a.queue.with_suffix('.stop').exists():sys.exit('Queue stopped while waiting')
    if time.time()>=deadline:sys.exit('Authorization ended while waiting')
    time.sleep(3)
for args in json.loads(a.queue.read_text()):
    if time.time()>=deadline or a.queue.with_suffix('.stop').exists():break
    print('START',args,flush=True)
    binary_args=[] if '--binary' in args else ['--binary',str(base/'versions/v6/athena')]
    r=subprocess.run([sys.executable,str(Path(__file__).with_name('run_case.py')),
        *args,'--cores',a.cores,*binary_args])
    print('END',r.returncode,args[0],flush=True)
