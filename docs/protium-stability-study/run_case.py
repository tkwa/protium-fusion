"""Run one bounded, CPU-affined Athena++ case and retain its complete log.

All processes end by the authorization deadline, 2026-10-31 08:22:56 UTC.
Use disjoint subsets of cores 0..6 when launching concurrent runs.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import signal
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser()
p.add_argument('name')
p.add_argument('--nr',type=int,default=256)
p.add_argument('--nt',type=int,default=1)
p.add_argument('--ell',type=int,default=2)
p.add_argument('--perturbation',type=int,choices=[0,1,2,3,4,5],default=0,help='0 radial speed; 1 density; 2 transverse speed; 3 ballistic time translation; 4 common-direction meridional sectors; 5 evolved reference phase')
p.add_argument('--timing-model',choices=['ballistic','reference'],default='ballistic')
p.add_argument('--phase-table',type=Path,default=Path('/var/tmp/tkwa-chatgpt/protium-stability/phase-reference.bin'))
p.add_argument('--reference-time-us',type=float,default=11.)
p.add_argument('--phase-taper-m',type=float,default=0.)
p.add_argument('--divisions',type=int,default=8,help='Meridional intervals for packet-direction surrogate')
p.add_argument('--angular-family',type=int,choices=[0,1],default=0,help='0 Legendre; 1 bounded cosine, both zero mean and unit angular RMS')
p.add_argument('--timing-us',type=float,default=None,help='RMS angular arrival-time error in microseconds; selects perturbation 3')
p.add_argument('--amplitude',type=float,default=0.)
p.add_argument('--test',type=int,default=0)
p.add_argument('--burn',type=int,default=1)
p.add_argument('--source-accuracy',type=float,default=.02)
p.add_argument('--preheat-fraction',type=float,default=0.)
p.add_argument('--preheat-floor-fraction',type=float,default=0.)
p.add_argument('--ambient',type=float,default=1e-9)
p.add_argument('--density-floor',type=float,default=1e-14)
p.add_argument('--cfl',type=float,default=.3)
p.add_argument('--tlim',type=float,default=2.28)
p.add_argument('--rmax',type=float,default=24.)
p.add_argument('--cores',default='0')
p.add_argument('--inputs',type=Path,default=Path('/var/tmp/tkwa-chatgpt/protium-stability/inputs-yi2'))
p.add_argument('--base',type=Path,default=Path('/var/tmp/tkwa-chatgpt/protium-stability'))
p.add_argument('--binary',type=Path,default=None)
p.add_argument('--output-interval',type=float,default=.1)
p.add_argument('--xorder',type=int,default=2)
p.add_argument('--max-wall-seconds',type=int,default=None)
p.add_argument('--geometry',choices=['spherical','cartesian'],default='spherical')
p.add_argument('--grid',type=int,default=32,help='Cartesian base grid cells per axis; divisible by four')
p.add_argument('--levels',type=int,default=4)
p.add_argument('--packet-count',type=int,default=128)
p.add_argument('--phase-packet-directions',type=int,choices=[0,1],default=0,help='Use one bulk direction per packet in the 3D phase-error calculation')
p.add_argument('--packet-times',type=Path,default=Path('/var/tmp/tkwa-chatgpt/protium-stability/packet-times-128-12345.txt'))
p.add_argument('--memory-limit-mib',type=int,default=1800)
a=p.parse_args()
if a.timing_us is not None:
    a.perturbation=5 if a.timing_model=='reference' else 3
    a.amplitude=a.timing_us/30.
cores=[int(x) for x in a.cores.split(',')]
assert cores and len(set(cores))==len(cores) and set(cores)<=set(range(7))
assert a.geometry=='cartesian' or a.nt>1 or a.amplitude==0, 'Angular perturbations require a multidimensional mesh'
out=a.base/'runs'/a.name
out.mkdir(parents=True,exist_ok=False)
boundary='polar_wedge' if a.nt>1 else 'reflecting'
inp=f'''<job>
problem_id = pulse
<output1>
file_type = hst
dt = 0.005
data_format = %24.16e
<output2>
file_type = vtk
variable = cons
dt = {a.output_interval}
<time>
cfl_number = {a.cfl}
nlim = -1
tlim = {a.tlim}
integrator = rk2
xorder = {a.xorder}
ncycle_out = 100
<mesh>
nx1 = {a.nr}
x1min = 0
x1max = {a.rmax}
ix1_bc = reflecting
ox1_bc = outflow
nx2 = {a.nt}
x2min = 0
x2max = 3.141592653589793
ix2_bc = {boundary}
ox2_bc = {boundary}
nx3 = 1
x3min = 0
x3max = 6.283185307179586
ix3_bc = periodic
ox3_bc = periodic
num_threads = {len(cores)}
<meshblock>
nx1 = {a.nr}
nx2 = {a.nt if len(cores)==1 else max(4,a.nt//len(cores))}
nx3 = 1
<hydro>
gamma = 1.6666666666666667
dfloor = {a.density_floor}
pfloor = 1e-18
efloor = 1e-18
sfloor = 0
<problem>
eos_path = {a.inputs/'eos.bin'}
profile = {a.inputs/'profile.txt'}
rates = {a.inputs/'rates.txt'}
burn = {a.burn}
test = {a.test}
ideal = {int(a.test>0)}
amplitude = {a.amplitude}
ell = {a.ell}
perturbation = {a.perturbation}
divisions = {a.divisions}
angular_family = {a.angular_family}
phase_table = {a.phase_table}
reference_time = {a.reference_time_us/30.}
phase_taper_radius = {a.phase_taper_m/30.}
ambient = {a.ambient}
source_accuracy = {a.source_accuracy}
preheat_fraction = {a.preheat_fraction}
preheat_floor_fraction = {a.preheat_floor_fraction}
'''
(out/'athinput').write_text(inp)
if a.geometry=='cartesian':
    assert a.grid%4==0 and a.grid>=16
    mesh=f'<mesh>\nrefinement = static\nnum_threads = {len(cores)}\n'
    for axis in [1,2,3]:
        mesh+=f'nx{axis} = {a.grid}\nx{axis}min = {-a.rmax}\nx{axis}max = {a.rmax}\nix{axis}_bc = outflow\nox{axis}_bc = outflow\n'
    mesh+='<meshblock>\n'+''.join(f'nx{axis} = {a.grid//4}\n' for axis in [1,2,3])
    for level in range(1,a.levels+1):
        extent=a.rmax/2**level
        mesh+=f'<refinement{level}>\nlevel = {level}\n'
        mesh+=''.join(f'x{axis}min = {-extent}\nx{axis}max = {extent}\n' for axis in [1,2,3])
    inp=inp[:inp.index('<mesh>')]+mesh+inp[inp.index('<hydro>'):]
    inp+=f'packet_count = {a.packet_count}\npacket_times = {a.packet_times}\nphase_packet_directions = {a.phase_packet_directions}\n'
    (out/'athinput').write_text(inp)
binary=a.binary or a.base/'athena/bin/athena'
source_dir=binary.parent if (binary.parent/'pulse.cpp').exists() else Path(__file__).parent
args={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
args.update(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
            started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            source_hashes={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in source_dir.glob('pulse*')})
args['input_hashes']={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in
                     [a.inputs/'eos.bin',a.inputs/'profile.txt',a.inputs/'rates.txt']}
if a.perturbation==5:
    args['phase_table_sha256']=hashlib.sha256(a.phase_table.read_bytes()).hexdigest()
    if a.geometry=='cartesian':args['packet_times_sha256']=hashlib.sha256(a.packet_times.read_bytes()).hexdigest()
(out/'config.json').write_text(json.dumps(args,indent=2)+'\n')
deadline=dt.datetime(2026,10,31,8,22,56,tzinfo=dt.timezone.utc).timestamp()
seconds=deadline-time.time()-10  # leave time for the forced-kill grace interval
if a.max_wall_seconds is not None:seconds=min(seconds,a.max_wall_seconds)
assert seconds>0,'Simulation authorization deadline has passed'
env=os.environ.copy()
env.update(OMP_NUM_THREADS=str(len(cores)),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_PROC_BIND='true',OMP_WAIT_POLICY='PASSIVE')
start=time.monotonic()
result=dict(returncode=None,reached_time_limit=False,fatal_error=True)
try:
    with (out/'stdout.log').open('w') as f:
        def memory_limit():
            import resource
            limit=a.memory_limit_mib*1024**2
            resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
        process=subprocess.Popen(['taskset','-c',a.cores,'timeout','--signal=TERM','--kill-after=10',str(int(seconds)),str(binary),'-i','athinput'],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=memory_limit if a.geometry=='cartesian' else None)
        stopped=False
        while process.poll() is None:
            if (out/'STOP').exists():
                stopped=True
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL)
                break
            time.sleep(1)
        returncode=process.wait()
    log=(out/'stdout.log').read_text()
    result=dict(returncode=returncode,elapsed_wall_s=time.monotonic()-start,operator_stopped=stopped,
                reached_time_limit=('Terminating on time limit' in log),
                fatal_error=('FATAL ERROR' in log or returncode!=0))
finally:
    (out/'execution.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(name=a.name,**result)),flush=True)
