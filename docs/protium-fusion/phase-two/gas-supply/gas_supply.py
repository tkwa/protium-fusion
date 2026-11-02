"""Fresh H2 transport screen in the saved Ti source; no integrated gas equilibrium.

Proton-only first ionization uses the Rudd fit as reproduced in Padovani et al.
2009, eqs.5--6. Fully stripped Ti, secondary ionization, charge exchange,
photoionization, wall recycling and neutral polarizability are omitted.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.special import roots_hermitenorm
from scipy.integrate import simpson,quad
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from electron_wall_orbits import Cylinder
from catalyst_penning import MP
from nonneutral_screen import MEV,E_CHARGE
KB=1.380649e-23;ME=9.1093837139e-31;A0=5.29177210544e-11

def sigma(E):
    x=np.maximum(ME/MP*E/(13.598e-6),1e-30)
    lo=.51*x**1.24;hi=(.71*np.log1p(x)+1.63)/x
    return 4*np.pi*A0*A0/(1/lo+1/hi)

def rate_profile(c,nr=257,ng=24):
    r=np.linspace(0,c.R,nr);eq=c.eq
    density=c.ns*np.interp(r/c.scale,eq['x'],eq['h'])
    x,w=roots_hermitenorm(ng);w=w/np.sqrt(2*np.pi)
    xx,yy,zz=np.meshgrid(x,x,x,indexing='ij');ww=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
    vv=np.stack([xx.ravel(),yy.ravel(),zz.ravel()])*np.sqrt(c.T*MEV/MP)
    omega=eq['kappa']*eq['xi']*E_CHARGE*c.B/MP
    rates=[]
    for rr in r:
        v=np.sqrt(vv[0]**2+(vv[1]+omega*rr)**2+vv[2]**2)
        rates.append(float(np.sum(ww*v*sigma(.5*MP*v*v/MEV))))
    return r,density,np.array(rates)

def exact_rate_profile(c,nr=257):
    """One-dimensional noncentral-Maxwell integral, independently suggested review."""
    r=np.linspace(0,c.R,nr);density=c.ns*np.interp(r/c.scale,c.eq['x'],c.eq['h'])
    thermal=np.sqrt(c.T*MEV/MP)
    omega=c.eq['kappa']*c.eq['xi']*E_CHARGE*c.B/MP
    def one(rr):
        u=omega*rr/thermal
        def integrand(w):
            if u<1e-9: f=np.sqrt(2/np.pi)*w*w*np.exp(-w*w/2)
            else: f=np.sqrt(2/np.pi)*w/(2*u)*(np.exp(-(w-u)**2/2)-np.exp(-(w+u)**2/2))
            return w*f*sigma(.5*MP*(thermal*w)**2/MEV)/1e-20
        return thermal*1e-20*quad(integrand,0,u+12,epsabs=1e-9,epsrel=1e-10)[0]
    return r,density,np.array([one(rr) for rr in r])

def flight_screen(r,nu,R,T=300.,n=131072,nquad=64,seed=813):
    rng=np.random.default_rng(seed)
    # Diffuse thermal emission: cosine-law angles and flux-weighted Maxwell speed.
    mu=np.sqrt(rng.random(n));az=2*np.pi*rng.random(n)
    vt=np.sqrt(1-mu*mu)*np.cos(az);vz=np.sqrt(1-mu*mu)*np.sin(az)
    speed=np.sqrt(-2*KB*T/(2*MP)*np.log(rng.random(n)*rng.random(n)))
    ell=2*R*mu/(mu*mu+vt*vt);dt=ell/speed
    # Midpoint quadrature allows a separate first-birth radius distribution.
    survive=np.ones(n);bins=np.linspace(0,R,81);hist=np.zeros(80)
    for j in range(nquad):
        l=ell*(j+.5)/nquad
        rr=np.sqrt(np.maximum(0,R*R-2*R*mu*l+(mu*mu+vt*vt)*l*l))
        h=np.interp(rr,r,nu)*dt/nquad
        birth=survive*(-np.expm1(-h));hist+=np.histogram(rr,bins,weights=birth)[0]
        survive*=np.exp(-h)
    a=float(np.mean(survive));p=1-a
    # Independent identical diffuse flights in an infinite cylinder. Getter acts
    # only upon reaching the next wall. It cannot remove an emitted molecule
    # before its first flight. End paths are negligible only for long channels.
    pumps=[]
    for sticking in [0.,.002,.01,.1,1.]:
        prob=p/(1-a*(1-sticking))
        pumps.append(dict(sticking=sticking,first_ionization_probability=prob,
                          first_birth_electrons_per_fresh_H2=prob,
                          conditional_full_stripping_electrons_per_fresh_H2=2*prob))
    hist=hist/hist.sum();cdf=np.cumsum(hist);cent=(bins[:-1]+bins[1:])/2
    return dict(samples=n,quadrature=nquad,seed=seed,mean_flight_s=float(dt.mean()),
        mean_flight_length_m=float(ell.mean()),mean_abs_axial_displacement_m=float(np.mean(abs(ell*vz))),
        first_flight_ionization_probability=p,standard_error=float(np.std(1-survive,ddof=1)/np.sqrt(n)),
        radius_quantiles_m={str(q):float(np.interp(q,cdf,cent)) for q in [.1,.5,.9]},
        radius_bins_m=bins.tolist(),first_birth_probability_distribution=hist.tolist(),pump_cases=pumps)

def diffusion(R,L,q,T,nu,alpha):
    v=np.sqrt(8*KB*T/(np.pi*2*MP));D=2*R*v/3;g=2*q/R;wall=alpha*v/(2*R);lam=nu+wall
    if lam==0:mean=g*L*L/(12*D);end=1.
    else:
        z=L/2*np.sqrt(lam/D);end=np.tanh(z)/z
        mean=g/lam*(1-end)
    return dict(mean_pressure_Pa=float(mean),end_fraction=float(end),ion_fraction=float(nu*mean/g),getter_fraction=float(wall*mean/g),
                Knudsen_diffusivity_m2_s=float(D),end_only_diffusive_residence_s=float(L*L/(12*D)))

def main():
    plant=next(p for p in json.loads((ROOT/'network-plant.json').read_text())['cases'] if p['assumptions']['B']==140 and p['assumptions']['heat']==.7 and p['assumptions']['end_in']==.999)
    src=plant['selected_source_profiles']['22049'];c=Cylinder(src);r,dens,rate=rate_profile(c)
    r2,d2,k2=exact_rate_profile(c)
    nu=dens*k2
    flights=flight_screen(r,nu,c.R);fine=flight_screen(r,dens*k2,c.R,nquad=128)
    F=src['reaction_rate_per_m_s']/(2*np.pi*c.R);area=plant['primary_pn_rate_s']/F
    gas=[]
    for q in [1e-12,1e-11,1e-10,1e-9,1e-8,3.815e-8]:
        molecules=q/(KB*300)/F
        gas.append(dict(q_Pa_m3_s_m2=q,H2_per_primary=molecules,
                    full_stripping_electrons_per_primary=2*molecules,
                    assumed_point1MeV_cost_per_primary=.1*2*molecules,
                    H2_mass_flow_kg_s=q*area/(KB*300)*2*MP))
    L=json.loads((ROOT/'lepton-escape-main.json').read_text())['active_length_m']
    rows=[]
    nuavg=simpson(2*r*nu,x=r)/(c.R*c.R)
    for length in [1.,10.,100.,L]:
        for alpha in [0.,.002,.01,.1]:
            rows.append(dict(length_m=length,sticking=alpha,**diffusion(c.R,length,1e-10,300,nuavg,alpha)))
    checks={
        'inventory_profile_relative_error':abs(simpson(2*np.pi*r*dens,x=r)/src['proton_inventory_per_m']-1),
        'GH24_rate_relative_error_against_1D_integral':float(np.max(abs(rate/k2-1))),
        'flight_midpoint_refinement_absolute_change':abs(flights['first_flight_ionization_probability']-fine['first_flight_ionization_probability']),
        'probability_closure_max_error':max(abs(row['end_fraction']+row['ion_fraction']+row['getter_fraction']-1) for row in rows),
        'diffusion_zero_sink_matches_analytic':abs(diffusion(c.R,L,1e-10,300,0,0)['mean_pressure_Pa']/(1e-10*L*L/(4*c.R*c.R*np.sqrt(8*KB*300/(np.pi*2*MP))))-1),
        'mean_chord_over_2R':fine['mean_flight_length_m']/(2*c.R),
    }
    out=dict(scope=__doc__,source=src,active_length_m=L,Ti_wall_area_m2=area,
        fit_source='Padovani, Galli & Glassgold 2009, A&A501,619,eqs5-6; underlying Rudd et al.1985 RMP57,965',
        cross_section_at_1MeV_m2=float(sigma(1.)),radial_volume_average_ionization_frequency_s=nuavg,
        profile=dict(radius_m=r.tolist(),proton_density_m3=dens.tolist(),ionization_rate_coefficient_m3_s=k2.tolist()),
        first_flight=fine,first_flight_coarse={k:v for k,v in flights.items() if k not in ['radius_bins_m','first_birth_probability_distribution']},
        source_budget=gas,uniform_diffusion_sensitivity=rows,
        no_plasma_no_getter=diffusion(c.R,L,1e-10,300,0,0),checks=checks,
        qualifications=['Proton velocity distribution is nonrelativistic Maxwellian with rigid rotation; finite-speed correction is not included.',
          'Only first H2 ionization by protons is modeled; it supplies one electron. Two times this first-ionization probability is conditional on neutral H2 persisting until that event. The unconditional fresh-H2 ceiling is two electrons per molecule; recycled wall gas and nuclear pair births are separate.',
          'The positive molecular ion and its dissociation products also require potential, rotation and collection accounting. The assigned .1MeV/e electrical cost is not established by this screen.',
          'No self-consistent gas chemistry or evolving lepton cloud. Wall outgassing rates are illustrative, not measured carbon-reactor properties.',
          'Repeated-flight model assumes diffuse fully accommodated reflection and an infinite uniform cylinder; it does not model saturation or desorption.',
          'Uniform diffusion is only a conductance illustration; radial ballistic calculation is used for first flight.',
          'A neutral leaving the high-field wall can polarize or dissociate; those effects are not included.'])
    here=Path(__file__).parent;(here/'gas-supply.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['cross_section_at_1MeV_m2','radial_volume_average_ionization_frequency_s','no_plasma_no_getter','checks']},indent=2))
    print(json.dumps({k:v for k,v in fine.items() if k not in ['radius_bins_m','first_birth_probability_distribution']},indent=2))
if __name__=='__main__':main()
