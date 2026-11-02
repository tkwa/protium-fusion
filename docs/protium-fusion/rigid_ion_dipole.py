"""Two-species rigid-centroid approximation for transverse dipole forcing.

Mutual and conducting-wall image forces follow exact linear force moments for
rigid translations of the two equilibrium density profiles. A rigid closure
omits profile distortion and kinetic damping. Its poles are model frequencies,
not the full warm Vlasov spectrum.
"""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.optimize import brentq
from electron_wall_orbits import Cylinder
from nonneutral_screen import E_CHARGE as EV,EPS0
from catalyst_penning import MP
ROOT=Path(__file__).parent

def modes(mass_line,charge_line,cyclotron,Kcross,R):
    image=np.outer(charge_line,charge_line)/(2*np.pi*EPS0*R*R)
    cross=Kcross*np.array([[1.,-1.],[-1.,1.]])
    A=cross+image;M=mass_line
    # With z=X+iY proportional to exp(i omega t):
    # [diag(M omega(omega+Omega)) + A] X=0.
    scale=max(cyclotron)
    C=np.block([[np.zeros((2,2)),np.eye(2)],[-A/M[:,None]/scale**2,-np.diag(cyclotron)/scale]])
    vals,vec=np.linalg.eig(C)
    result=[]
    for i in np.argsort(vals.real):
        v=vec[:2,i];v=v/np.max(np.abs(v));w=vals[i]*scale
        result.append({'omega_real_rad_s':float(w.real),'omega_imag_rad_s':float(w.imag),
                       'relative_centroid_amplitude_real':v.real.tolist(),'relative_centroid_amplitude_imag':v.imag.tolist()})
    return result,A,image,cross

def run(root):
    plant=next(x for x in json.loads((root/'network-plant.json').read_text())['cases'] if x['assumptions']['B']==140 and x['assumptions']['heat']==.7 and x['assumptions']['end_in']==.999)
    c=Cylinder(plant['selected_source_profiles']['22049']);N=np.array([c.src['proton_inventory_per_m'],c.src['catalyst_inventory_per_m']]);q=EV*np.array([1.,22.]);mass=MP*np.array([1.,c.eq['M']]);lam=N*q;ml=N*mass;Om=q*c.B/mass
    overlap=2*np.pi*c.scale**2*c.ns**2*np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc']);K=q.prod()*overlap/(2*EPS0)
    roots,A,image,cross=modes(ml,lam,Om,K,c.R)
    def omegaE(r):return -(-c.T*1e6/c.scale*float(c.spline(r/c.scale,1)))/(c.B*r)
    rgrid=np.linspace(c.R*1e-6,c.R,2001);wgrid=np.array([omegaE(r) for r in rgrid])
    for row in roots:
        w=row['omega_real_rad_s'];brackets=np.flatnonzero((wgrid[:-1]-w)*(wgrid[1:]-w)<0);rs=[brentq(lambda r:omegaE(r)-w,rgrid[i],rgrid[i+1],xtol=1e-15) for i in brackets]
        row['lepton_ExB_matching_radii_m']=rs
        row['reaction_fraction_inside_matching_radius']=[float(np.interp(r/c.scale,c.eq['x'],c.eq['ihc'])/np.interp(c.src['xwall'],c.eq['x'],c.eq['ihc'])) for r in rs]
    response=[]
    for w in [-1.56e9,-1.51e9,-1e9,-.15e9]:
        x=np.linalg.solve(-np.diag(ml*w*(w+Om))-A,lam)
        response.append({'omega_rad_s':w,'centroid_m_per_V_m':x.tolist()})
    return {'scope':__doc__,'source':c.src,'line_number_m-1':N.tolist(),'line_charge_C_m':lam.tolist(),'line_mass_kg_m':ml.tolist(),
            'cyclotron_rad_s':Om.tolist(),'density_overlap_m-5':float(overlap),'cross_force_derivative_N_m2':float(K),
            'image_force_matrix_N_m2':image.tolist(),'total_force_matrix_N_m2':A.tolist(),'modes':roots,'responses':response,
            'lepton_ExB_rotation_axis_rad_s':float(wgrid[0]),'lepton_ExB_rotation_wall_rad_s':float(wgrid[-1]),
            'limits':['Rigid translations only; nonuniform warm plasma profiles can distort and phase mix.',
                      'The mutual force coefficient follows overlap of both number densities, without a fitted spring.',
                      'The image force is the exact leading dipole field of the grounded cylindrical boundary.',
                      'A frequency crossing is a reason to calculate the coupled response; it is not a proof of prohibitive damping.']}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=ROOT/'rigid-ion-dipole.json');a=ap.parse_args();d=run(a.root);a.out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(json.dumps(d['modes'],indent=2))
