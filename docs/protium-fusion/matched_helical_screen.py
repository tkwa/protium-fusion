"""Local fixed-total-log matching screen for the weak helical collision model.

Uses unmagnetized total minus filtered unmagnetized plus filtered helical.
Both physical long/short contributions are positive for a positive remaining
short-range logarithm. Filter convergence does not qualify static screening,
straight-ion orbits, source homogeneity or the chosen total logarithm.
"""
from pathlib import Path
from math import pi,sqrt
import argparse,json,time
import numpy as np
from finite_gyro_spectrum import helical_coefficients,unmagnetized_gyro_average
from electron_wall_orbits import Cylinder,ME,E_CHARGE as EV,EPS0
from lepton_collision_screen import local_state


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent)
    ap.add_argument('--outdir',type=Path,default=Path(__file__).parent)
    ap.add_argument('--cut',type=float,default=1.);ap.add_argument('--cold',action='store_true')
    ap.add_argument('--nk',type=int,default=96);ap.add_argument('--nmu',type=int,default=64);ap.add_argument('--nphi',type=int,default=64)
    ap.add_argument('--log-total',type=float,default=15.);args=ap.parse_args();start=time.monotonic()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049'])
    cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])
    r=float(np.interp(.5,cdf,c.eq['x'])*c.scale);local=local_state(c,r)
    Omega=EV*c.B/ME;U=local['relative_azimuthal_m_s'];kappa=1/local['ion_static_screening_m']
    reference_E=3467.60505;vp=0. if args.cold else sqrt(4*reference_E*EV/(3*ME));vz=0. if args.cold else sqrt(2*reference_E*EV/(3*ME))
    mTi,_,_,TTi=local['species'][1];vscale=sqrt(2*TTi*EV/mTi)
    kcut=args.cut*Omega/vscale;rows=[]
    totalA=np.zeros(3);totalD=np.zeros((3,3));baseA=totalA.copy();baseD=totalD.copy()
    def measures(A,D):
        S=np.diag([1/vscale,1/vscale,1/vscale**2]);ev=np.linalg.eigvalsh(S@D@S)
        return {'drift':A.tolist(),'diffusion':D.tolist(),
            'perpendicular_heating_eV_s':float(ME/EV*A[2]),
            'parallel_heating_eV_s':float(ME/EV*(vz*A[1]+.5*D[1,1])),
            'total_heating_eV_s':float(ME/EV*(A[2]+vz*A[1]+.5*D[1,1])),
            'electron_radial_drift_m_s':float(A[0]/(-Omega)),
            'scaled_covariance_eigenvalues_s-1':ev.tolist()}
    for mass,density,z,T in local['species']:
        vth=sqrt(2*T*EV/mass);G0=4*pi*(z*EV**2/(4*pi*EPS0))**2*density/ME**2
        helix=helical_coefficients(vp,vz,U,vth,Omega,ME/mass,kappa,kcut,G0,args.nk,args.nmu,args.nphi)
        logarithm=helix['spectral_log']
        if logarithm>args.log_total:raise ValueError('Filtered log exceeds the prescribed total.')
        total=unmagnetized_gyro_average(vp,vz,U,vth,ME/mass,args.log_total,G0)
        short=unmagnetized_gyro_average(vp,vz,U,vth,ME/mass,args.log_total-logarithm,G0)
        A=short['drift']+helix['drift'];D=short['diffusion']+helix['diffusion']
        totalA+=A;totalD+=D;baseA+=total['drift'];baseD+=total['diffusion']
        rows.append({'charge':z,'mass_kg':mass,'density_m-3':density,'vth_m_s':vth,
            'filtered_logarithm':logarithm,'short_logarithm':args.log_total-logarithm,'harmonic_limit':helix['harmonic_limit'],
            'matched':measures(A,D),'unmagnetized':measures(total['drift'],total['diffusion']),
            'filtered_helical':measures(helix['drift'],helix['diffusion'])})
    out={'scope':__doc__,'local_state':local,'vperp_m_s':vp,'vz_m_s':vz,'gyrofrequency_rad_s':Omega,
        'cut_factor':args.cut,'kappa_m-1':kappa,'kcut_m-1':kcut,'log_total':args.log_total,
        'quadrature':{'nk':args.nk,'nmu':args.nmu,'nphi':args.nphi},'species':rows,
        'matched':measures(totalA,totalD),'unmagnetized':measures(baseA,baseD),'seconds':time.monotonic()-start}
    state='cold' if args.cold else 'warm'
    name=f'matched-helical-{state}-cut{args.cut:g}-{args.nk}-{args.nmu}-{args.nphi}.json'
    (args.outdir/name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'file':name,'seconds':out['seconds'],'matched':out['matched'],'unmagnetized':out['unmagnetized']}),flush=True)


if __name__=='__main__':main()
