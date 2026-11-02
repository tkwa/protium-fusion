"""Sensitivity to consistently removing a distant collision logarithmic interval.

Increasing kappa by f and reducing the total logarithm by ln(f) leaves the
short-range contribution almost unchanged. This diagnoses the influence of
long wavelengths; it does not calculate dynamic screening or bound its effect.
"""
from pathlib import Path
from math import sqrt,pi,log
import argparse,json,time
import numpy as np
from streaming_helical_spectrum import helical_coefficients
from finite_gyro_spectrum import unmagnetized_gyro_average
from electron_wall_orbits import ME,E_CHARGE as EV,EPS0


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args();start=time.monotonic()
    base=json.loads((args.root/'matched-streaming-warm-cut2-512-192-32.json').read_text());l=base['local_state'];rows=[]
    for cold in [True,False]:
        vp=0. if cold else base['vperp_m_s'];vz=0. if cold else base['vz_m_s']
        for factor,nk in [(1.,512),(2.,512),(10.,512),(10.,768)]:
            tick=time.monotonic();A=np.zeros(3);D=np.zeros((3,3));species=[]
            for mass,density,z,T in l['species']:
                vth=sqrt(2*T*EV/mass);G0=4*pi*(z*EV**2/(4*pi*EPS0))**2*density/ME**2
                h=helical_coefficients(vp,vz,l['relative_azimuthal_m_s'],vth,base['gyrofrequency_rad_s'],ME/mass,base['kappa_m-1']*factor,base['kcut_m-1'],G0,nk,192,32)
                remaining=base['log_total']-log(factor)-h['spectral_log']
                if remaining<0:raise ValueError('Negative complementary logarithm')
                s=unmagnetized_gyro_average(vp,vz,l['relative_azimuthal_m_s'],vth,ME/mass,remaining,G0)
                A+=h['drift']+s['drift'];D+=h['diffusion']+s['diffusion'];species.append({'charge':z,'filtered_log':h['spectral_log'],'complementary_log':remaining})
            item={'cold':cold,'kappa_multiplier':factor,'total_log':base['log_total']-log(factor),'nk':nk,'nmu':192,'nphi':32,'drift':A.tolist(),'diffusion':D.tolist(),'perpendicular_heating_eV_s':ME/EV*A[2],'parallel_heating_eV_s':ME/EV*(vz*A[1]+.5*D[1,1]),'total_heating_eV_s':ME/EV*(A[2]+vz*A[1]+.5*D[1,1]),'electron_radial_drift_m_s':A[0]/(-base['gyrofrequency_rad_s']),'species':species,'seconds':time.monotonic()-tick};rows.append(item)
            print('cold',cold,'kappa multiplier',factor,'nk',nk,'seconds',item['seconds'],'Pperp/Ppar',item['perpendicular_heating_eV_s']/1e6,item['parallel_heating_eV_s']/1e6,flush=True)
    out={'scope':__doc__,'reference_source':'matched-streaming-warm-cut2-512-192-32.json','cases':rows,'seconds':time.monotonic()-start}
    (args.outdir/'screening-sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
