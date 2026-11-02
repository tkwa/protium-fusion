"""Interpolated electron/positron EOS for repeated pulse calculations.

The source of the table is pynucastro's ElectronEOS, including positron
rest-mass creation energy. Density coordinate is rho*Ye (g/cm^3).
All internal quantities are cgs. Ions and radiation are analytic.
Interpolation is an approximation; validate against the source EOS.
"""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.interpolate import RectBivariateSpline

TABLE=Path(__file__).with_name('electron-eos.npz')
KB_CGS=1.380649e-16
NA=6.02214076e23
ARAD_CGS=7.565733250033928e-15


class FastEOS:
    def __init__(self,path=TABLE):
        data=np.load(path)
        self.ld=data['log_density'];self.lt=data['log_temperature']
        self.ep=RectBivariateSpline(self.ld,self.lt,data['log_pressure'],s=0)
        self.ee=RectBivariateSpline(self.ld,self.lt,data['log_energy_density'],s=0)
    def pe_state(self,rho,T,comp):
        ye=comp.ye; yi=1/comp.abar
        return self.state(rho,T,ye,yi)
    def state(self,rho,T,ye,yi):
        """Scalar or array state from electron and ion molar abundances."""
        ld=np.log(rho*ye);lt=np.log(T)
        if np.any((ld<self.ld[0])|(ld>self.ld[-1])|(lt<self.lt[0])|(lt>self.lt[-1])):
            raise ValueError(f'EOS table exceeded: rho={rho}, Ye={ye}, T={T}')
        p_e=np.exp(self.ep.ev(ld,lt));v_e=np.exp(self.ee.ev(ld,lt))
        pe_rho=p_e*self.ep.ev(ld,lt,dx=1)/rho
        pe_T=p_e*self.ep.ev(ld,lt,dy=1)/T
        ue=v_e/rho
        ue_rho=ue*(self.ee.ev(ld,lt,dx=1)-1)/rho
        ue_T=ue*self.ee.ev(ld,lt,dy=1)/T
        pi=rho*KB_CGS*NA*T*yi;ui=1.5*pi/rho
        pr=ARAD_CGS*T**4/3;ur=3*pr/rho
        pressure=p_e+pi+pr;energy=ue+ui+ur
        p_rho=pe_rho+pi/rho
        p_T=pe_T+pi/T+4*pr/T
        u_rho=ue_rho-ur/rho
        u_T=ue_T+ui/T+4*ur/T
        gamma=rho/pressure*(p_rho+p_T*(pressure/rho**2-u_rho)/u_T)
        return SimpleNamespace(p=pressure,e=energy,dp_drho=p_rho,dp_dT=p_T,
          de_drho=u_rho,de_dT=u_T,c_v=u_T,gamma1=gamma,eta=float('nan'))


def build_row(args):
    d,temps=args
    import pynucastro as p
    from pynucastro.eos.electron_eos import ElectronEOS
    eos=ElectronEOS();comp=p.Composition(['p']);comp.set_all(1.)
    pressures=[];energies=[]
    for t in temps:
        e,pos=eos.pe_state(d,t,comp,compute_derivs=False)
        pressures.append(e.p+pos.p)
        energies.append(d*(e.e+pos.e))
    return np.log(pressures),np.log(energies)


if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    import time
    ld=np.linspace(np.log(.01),np.log(1e7),73)
    lt=np.linspace(np.log(1e7),np.log(3e9),101)
    ep=[];ee=[];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,(pp,uu) in enumerate(pool.map(build_row,[(d,np.exp(lt)) for d in np.exp(ld)])):
            ep.append(pp);ee.append(uu)
            if i%8==0: print(f'{i+1}/{len(ld)} EOS rows; {time.monotonic()-start:.1f}s',flush=True)
    np.savez_compressed(TABLE,log_density=ld,log_temperature=lt,
      log_pressure=ep,log_energy_density=ee,
      source='pynucastro 3.1.0 ElectronEOS; ions/radiation analytic; no Coulomb corrections')
    print(TABLE,flush=True)
