"""Install the retained problem/EOS adapters into the pinned Athena++ checkout."""
import argparse,os,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('checkout',type=Path);p.add_argument('--cores',default='6');p.add_argument('--flux',choices=['hllc','hlle'],default='hllc');p.add_argument('--source-dir',type=Path);p.add_argument('--geometry',choices=['spherical_polar','cartesian'],default='spherical_polar');a=p.parse_args()
cores=[int(c) for c in a.cores.split(',')]
assert cores and len(set(cores))==len(cores) and set(cores)<=set(range(7))
root=Path(__file__).resolve().parent
source_dir=a.source_dir or root
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.checkout,text=True).strip()
assert commit=='053bc738d5c85dac717326b4338b34c7db4e715d',commit
patch=source_dir/'framework.patch' if (source_dir/'framework.patch').exists() else root/'athena-framework.patch'
already=subprocess.run(['git','apply','--reverse','--check',str(patch)],cwd=a.checkout,capture_output=True).returncode==0
if not already:
    subprocess.run(['git','apply','--check',str(patch)],cwd=a.checkout,check=True)
    subprocess.run(['git','apply',str(patch)],cwd=a.checkout,check=True)
for name,dest in [('pulse.cpp','src/pgen'),('pulse_eos.cpp','src/eos/general'),('pulse_eos.hpp','src/eos/general')]:
    shutil.copy2(source_dir/name,a.checkout/dest/name)
for name in ['src/pgen/pulse.cpp','src/eos/general/pulse_eos.cpp']:
    os.utime(a.checkout/name,None)  # Athena's makefile does not track our new header.
subprocess.run(['python3','configure.py','--prob=pulse',f'--coord={a.geometry}',
    '--eos=general/pulse_eos',f'--flux={a.flux}','--nscalars=2','-omp'],cwd=a.checkout,check=True)
subprocess.run(['taskset','-c',a.cores,'make',f'-j{len(cores)}'],cwd=a.checkout,check=True)
