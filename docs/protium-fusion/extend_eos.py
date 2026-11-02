"""Extend the same source EOS to colder, dilute expanded surface zones."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import time
import numpy as np
from fast_eos import TABLE,build_row

if __name__=='__main__':
    old=np.load(TABLE)
    oldld=old['log_density'];oldlt=old['log_temperature']
    lt_low=np.linspace(np.log(1e5),oldlt[0],82)[:-1]
    ld_low=np.linspace(np.log(1e-6),oldld[0],33)[:-1]
    ld=np.concatenate([ld_low,oldld]);lt=np.concatenate([lt_low,oldlt])
    ep=np.empty((len(ld),len(lt)));ee=ep.copy()
    ep[len(ld_low):,len(lt_low):]=old['log_pressure']
    ee[len(ld_low):,len(lt_low):]=old['log_energy_density']
    tasks=[(d,np.exp(lt)) for d in np.exp(ld_low)]
    tasks += [(d,np.exp(lt_low)) for d in np.exp(oldld)]
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,(pp,uu) in enumerate(pool.map(build_row,tasks)):
            if i<len(ld_low): ep[i]=pp;ee[i]=uu
            else: ep[i,:len(lt_low)]=pp;ee[i,:len(lt_low)]=uu
            if i%10==0: print(f'{i+1}/{len(tasks)} extension rows; {time.monotonic()-start:.1f}s',flush=True)
    backup=TABLE.with_name('electron-eos-original.npz')
    if not backup.exists(): backup.write_bytes(TABLE.read_bytes())
    np.savez_compressed(TABLE,log_density=ld,log_temperature=lt,
      log_pressure=ep,log_energy_density=ee,
      source='pynucastro 3.1.0 ElectronEOS; extended low-density/temperature range; no Coulomb corrections')
