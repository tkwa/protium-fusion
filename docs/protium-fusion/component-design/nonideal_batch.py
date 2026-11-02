"""Run the eight authorized independent Coulomb stress cases, at most two at once."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import argparse,hashlib,json,os,subprocess,time
HERE=Path(__file__).resolve().parent
PYTHON='/tmp/protium-python-interpreters/cpython-3.10.20-linux-x86_64-gnu/bin/python3.10'
CASES=[('A',1),('A',3),('B1',1),('B1',2),('B2',2),('B2',2.5),('C13',1),('C13',4)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONPATH='/tmp/protium-python',MPLCONFIGDIR='/tmp/protium-mpl',OPENBLAS_NUM_THREADS='1',NUMBA_CACHE_DIR=str(a.out/'.numba-cache'))
 start=time.monotonic();results=[]
 def one(stage,h):
  name=f'{stage}-H{h}';cmd=[PYTHON,str(HERE/'nonideal_run.py'),'--stage',stage,'--h-per-catalyst',str(h),'--out',str(a.out)]
  with (a.out/(name+'.log')).open('w') as log:
   p=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
  return {'stage':stage,'h_per_catalyst':h,'returncode':p.returncode,'elapsed_from_batch_start_s':time.monotonic()-start,'log':str(a.out/(name+'.log'))}
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures={pool.submit(one,*case):case for case in CASES}
  for f in as_completed(futures):
   try:r=f.result()
   except Exception as e:r={'case':futures[f],'error':repr(e)}
   results.append(r);(a.out/'batch-progress.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(r),flush=True)
 (a.out/'batch-result.json').write_text(json.dumps({'cases':results,'max_parallel':2,'wall_s':time.monotonic()-start,'batch_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
if __name__=='__main__':main()
