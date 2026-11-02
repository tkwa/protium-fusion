"""Encounter-ordering probes in the saved rigidly rotating ion equilibrium.

The co-rotating thermal equation has effective cyclotron frequency Omega_i+2w
and force q E + m w(Omega_i+w)r. These are local ordering diagnostics, not a
curved-orbit collision operator or a dynamic dielectric calculation.
"""
from pathlib import Path
from math import sqrt
import argparse,json,time
import numpy as np
from electron_wall_orbits import Cylinder,ME,E_CHARGE as EV
from lepton_collision_screen import local_state


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent)
    ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args();start=time.monotonic()
    plant=next(p for p in json.loads((args.root/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);cdf=c.eq['ihc']/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc']);rows=[]
    warm_vz=sqrt(2*3467.60505*EV/(3*ME))
    for fraction in [.1,.5,.9,.99,.999]:
        r=float(np.interp(fraction,cdf,c.eq['x'])*c.scale);l=local_state(c,r);w=l['ion_azimuthal_m_s']/r
        wE=l['ExB_azimuthal_m_s']/r;relative_w=wE-w
        E=l['radial_E_V_m'];Eprime=-c.T*1e6/c.scale**2*float(c.spline(r/c.scale,2))
        row={'reaction_fraction_inside':fraction,'local_state':l,'ion_rotation_rad_s':w,'electron_GC_rotation_rad_s':wE,'relative_GC_rotation_rad_s':relative_w,'electron_cyclotron_rotation_fraction':2*w/(-EV*c.B/ME),'species':[]}
        for mass,density,z,T in l['species']:
            Omega=z*EV*c.B/mass;effective=Omega+2*w;vth=sqrt(2*T*EV/mass)
            a0=z*EV/mass*E+w*(Omega+w)*r
            fr=z*EV/mass*Eprime+w*(Omega+w);ft=a0/r
            lnprime=mass*a0/(T*EV);lnsecond=mass*fr/(T*EV)
            item={'charge':z,'mass_kg':mass,'laboratory_cyclotron_rad_s':Omega,'rotating_thermal_cyclotron_rad_s':effective,'thermal_speed_m_s':vth,'effective_radial_acceleration_m_s2':a0,'force_gradient_eigenvalues_s-2':[fr,ft],'log_density_first_derivative_m-1':lnprime,'log_density_second_derivative_m-2':lnsecond,'density_scale_m':1/abs(lnprime),'encounters':[]}
            for label,vz in [('cold',0.),('warm',warm_vz)]:
                speed=sqrt(l['relative_azimuthal_m_s']**2+vz*vz+vth*vth)
                for mult in [.1,1.]:
                    b=mult*l['ion_static_screening_m'];tau=b/speed
                    item['encounters'].append({'electron_probe':label,'screening_length_fraction':mult,'distance_m':b,'characteristic_relative_speed_m_s':speed,'time_s':tau,'laboratory_ion_turn_rad':Omega*tau,'rotating_thermal_turn_rad':effective*tau,'force_gradient_time_squared':[fr*tau*tau,ft*tau*tau],'constant_ion_acceleration_displacement_over_distance':abs(a0)*tau*tau/(2*b),'electron_GC_curvature_displacement_over_distance':relative_w**2*r*tau*tau/(2*b),'log_density_change_linear_over_distance':abs(lnprime)*b})
            row['species'].append(item)
        rows.append(row)
    out={'scope':__doc__,'coordinate_convention':'Signed rigid rotation w<0; positive-ion laboratory Omega>0; thermal u=v-w cross r. Electron GC follows the leading local E cross B orbit.','rows':rows,'seconds':time.monotonic()-start}
    (args.outdir/'rotating-collision-scales.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('seconds',out['seconds'])
    for row in rows:
        print('fraction',row['reaction_fraction_inside'],'rotation',row['ion_rotation_rad_s'],'electron correction',row['electron_cyclotron_rotation_fraction'])
        for sp in row['species']:
            print(' Z',sp['charge'],'Omega effective',sp['rotating_thermal_cyclotron_rad_s'],'gradients',sp['force_gradient_eigenvalues_s-2'])
            for e in sp['encounters']:
                if e['screening_length_fraction']==1:print(' ',e['electron_probe'],'rot turn',e['rotating_thermal_turn_rad'],'ion/electron curvature ratios',e['constant_ion_acceleration_displacement_over_distance'],e['electron_GC_curvature_displacement_over_distance'],'density log variation',e['log_density_change_linear_over_distance'])


if __name__=='__main__':main()
