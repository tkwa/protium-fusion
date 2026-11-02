"""Independent bounded check of the fresh-H2 proton-ionization screen.

Run with the research Python environment, for example:
  PYTHONPATH=/tmp/protium-python OPENBLAS_NUM_THREADS=1 python3 independent_gas_check.py

Uses the saved cylindrical equilibrium through Cylinder, but separately evaluates
its shifted Maxwell average by a scalar speed integral, and replaces random
molecular flights with Legendre angular and generalized-Laguerre speed quadrature.
The ionization fit is independently transcribed from Padovani et al. 2009,
A&A 501,619, equations 5--6 (Rudd et al. 1985). These are conditional, prescribed
fields and first-ionization calculations, not a gas/electron equilibrium.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import scipy
from scipy.integrate import quad,simpson
from scipy.special import roots_legendre,roots_genlaguerre
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from electron_wall_orbits import Cylinder
from catalyst_penning import MP
from nonneutral_screen import MEV,E_CHARGE
KB=1.380649e-23
M_E=9.1093837139e-31
A0=5.29177210544e-11


def cross_section(E_MeV):
    x=np.maximum(M_E/MP*E_MeV/(13.598e-6),1e-30)
    low=.51*x**1.24
    high=(.71*np.log1p(x)+1.63)/x
    return 4*np.pi*A0*A0*low*high/(low+high)


def independent_rates(c,r):
    thermal_speed=np.sqrt(c.T*MEV/MP)
    omega=c.eq['kappa']*c.eq['xi']*E_CHARGE*c.B/MP
    result=[]
    # The speed of a 3D unit-normal vector displaced by u has radial density
    # sqrt(2/pi)*w/(2u)*[exp(-(w-u)^2/2)-exp(-(w+u)^2/2)].
    # At u=0 this tends to the ordinary Maxwell density.
    for rr in r:
        u=omega*rr/thermal_speed
        def integrand(w):
            if u<1e-9:
                density=np.sqrt(2/np.pi)*w*w*np.exp(-w*w/2)
            else:
                density=np.sqrt(2/np.pi)*w/(2*u)*(np.exp(-(w-u)**2/2)-np.exp(-(w+u)**2/2))
            return w*density*cross_section(.5*MP*(thermal_speed*w)**2/MEV)/1e-20
        val,error=quad(integrand,0,u+12,epsabs=1e-9,epsrel=1e-10)
        result.append(thermal_speed*1e-20*val)
    return np.array(result)


def angle_quadrature(R,nt,nphi):
    q,w=roots_legendre(nt)
    mu=np.sqrt((q+1)/2)[:,None]  # mu^2 uniform for cosine-law wall emission
    wu=w/2
    q,w=roots_legendre(nphi)
    ph=np.pi*(q+1)[None,:]
    wp=w/2
    transverse=mu**2+(1-mu**2)*np.cos(ph)**2
    impact2=R*R*(1-mu**2)*np.cos(ph)**2/transverse
    half_chord=R*mu/np.sqrt(transverse)
    return mu,transverse,impact2,half_chord,wu[:,None]*wp[None,:]


def deterministic_survival(r,nu,R,nt,nphi,nspeed,nline):
    mu,a2,b2,zmax,angle_weight=angle_quadrature(R,nt,nphi)
    ql,wl=roots_legendre(nline)
    integrated_rate=np.zeros_like(a2)
    for xx,ww in zip(ql,wl):
        rr=np.sqrt(b2+(zmax*xx)**2)
        integrated_rate+=ww*np.interp(rr,r,nu)
    integrated_rate*=zmax/np.sqrt(a2)
    # t=m*v^2/(2kT) has density t*exp(-t) in a Maxwell wall flux.
    x,w=roots_genlaguerre(nspeed,1)
    speeds=np.sqrt(KB*300/MP*x)
    survive=sum(ww*np.exp(-integrated_rate/v) for ww,v in zip(w,speeds))
    ell=2*R*mu/a2
    return {'n_mu_squared':nt,'n_azimuth':nphi,'n_speed':nspeed,'n_chord':nline,
        'first_ionization_probability':float(np.sum(angle_weight*(1-survive))),
        'mean_chord_over_2R':float(np.sum(angle_weight*ell)/(2*R))}


def birth_potential(c,r,nu,nt=64,nphi=128,nspeed=48,nline=128):
    # Uniform path cells preserve chronological hazard integration independently
    # of the Monte Carlo birth histogram. Each cell is represented at its midpoint.
    mu,a2,b2,zmax,aw=angle_quadrature(c.R,nt,nphi)
    a2=a2.ravel();b2=b2.ravel();zmax=zmax.ravel();aw=aw.ravel()
    ell=2*zmax/np.sqrt(a2)
    x,sw=roots_genlaguerre(nspeed,1);speed=np.sqrt(KB*300/MP*x)
    survival=np.ones((a2.size,nspeed));total=0.;weighted_potential=0.
    for j in range(nline):
        z=zmax*(-1+2*(j+.5)/nline);rr=np.sqrt(b2+z*z)
        hazard=np.interp(rr,r,nu)[:,None]*ell[:,None]/(nline*speed[None,:])
        birth=survival*(-np.expm1(-hazard))
        birth_per_angle=birth@sw
        total+=float(aw@birth_per_angle)
        weighted_potential+=float(aw@(birth_per_angle*c.potential_MV(rr)))
        survival*=np.exp(-hazard)
    return {'n_mu_squared':nt,'n_azimuth':nphi,'n_speed':nspeed,'n_chord':nline,
            'first_ionization_probability':total,
            'conditional_birth_weighted_potential_MV':weighted_potential/total}


def main():
    source_file=HERE/'gas-supply.json';saved=json.loads(source_file.read_text())
    c=Cylinder(saved['source']);profile=saved['profile']
    r=np.array(profile['radius_m']);n=np.array(profile['proton_density_m3'])
    ks=np.array(profile['ionization_rate_coefficient_m3_s']);ki=independent_rates(c,r);nu=n*ki
    flights=[deterministic_survival(r,nu,c.R,*args) for args in
             [(40,80,32,80),(80,160,48,128),(160,320,64,192)]]
    births=[birth_potential(c,r,nu,nline=count) for count in [128,256]]
    vbar=np.sqrt(8*KB*300/(np.pi*2*MP))
    rr=np.linspace(1e-10,c.R,10001)
    E=-c.T*MEV/(E_CHARGE*c.scale)*c.spline.derivative()(rr/c.scale)
    alpha=8.9568e-41 # measured H2 molecular polarizability, NIST 2009, 273 K
    output={
      'scope':__doc__,
      'versions':{'numpy':np.__version__,'scipy':scipy.__version__},
      'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      [source_file,Path(__file__),ROOT/'electron_wall_orbits.py',ROOT/'catalyst_penning.py']},
      'rates':{'axis_independent_m3_s':float(ki[0]),'axis_saved_m3_s':float(ks[0]),
               'maximum_relative_saved_difference':float(np.max(abs(ks/ki-1))),
               'independent_volume_mean_frequency_s':float(simpson(2*r*nu,x=r)/(c.R*c.R)),
               'radial_rate_coefficients_m3_s':ki.tolist()},
      'deterministic_flight_refinements':flights,
      'birth_potential_refinements':births,
      'monte_carlo_comparison':{'saved_probability':saved['first_flight']['first_flight_ionization_probability'],
               'saved_standard_error':saved['first_flight']['standard_error'],
               'difference_in_saved_standard_errors':(saved['first_flight']['first_flight_ionization_probability']-flights[-1]['first_ionization_probability'])/saved['first_flight']['standard_error']},
      'kinetic_theory_moments':{'mean_chord_m':2*c.R,'mean_time_s':2*c.R/vbar,
               'saved_mean_time_s':saved['first_flight']['mean_flight_s']},
      'neutral_polarization_scale':{
          'source':'Schmidt, Moldover & May 2009, Reference measurements of hydrogen\u2019s dielectric permittivity; https://www.nist.gov/publications/reference-measurements-hydrogens-dielectric-permittivity',
          'polarizability_F_m2':alpha,'maximum_E_V_m':float(E.max()),
          'maximum_E_radius_m':float(rr[np.argmax(E)]),
          'max_induced_dipole_energy_eV':float(.5*alpha*E.max()**2/E_CHARGE),
          'wall_to_max_field_energy_over_kT_300':float(.5*alpha*(E.max()**2-E[-1]**2)/(KB*300))},
      'interpretation':[
          'The birth potential is a gross energy-per-unit-charge scale. Co-born positive and negative potential energies cancel initially; this is not an automatic net electrical loss.',
          'Positive H2+/H+ products require consistent inventory, rotation, extraction and recovery accounting. An assumed complete cost per electron is not validated by this screen.',
          'The Knudsen no-plasma profile requires gas-gas mean free path much larger than channel diameter, and cannot be extrapolated indiscriminately to larger outgassing.',
          'Ordinary induced-dipole deflection is small at 300 K in this saved smooth field. Dissociation, high-Z ionization, charge exchange, secondary particles and wall recycling remain separate.',
          'Proton-only first ionization is not an unconditional lower bound on free-electron births when competing molecular chemistry is admitted. Two free electrons per fresh H2 is a fresh-inventory ceiling, not a recycling avalanche ceiling.'
      ]}
    (HERE/'independent-gas-check.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['scope','rates']},indent=2))
    print(json.dumps({k:v for k,v in output['rates'].items() if k!='radial_rate_coefficients_m3_s'},indent=2))

if __name__=='__main__':main()
