"""Selected-mixture material collection screen; no new hydrodynamics.
Actual endpoint species and uniform-dm/dv cells are propagated with the saved
network's one-body weak decays. Measured electronic stopping and ZBL elastic
transfer are reused unchanged from the earlier component screen.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math, sys
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import expm

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

def quantile(x, weight):
    order = np.argsort(x)
    c = np.cumsum(weight[order]); x = x[order]
    return {str(p): float(x[min(np.searchsorted(c, p*c[-1]), len(x)-1)])
            for p in [.5, .9, .99, .999, .9999]}

def setup(root):
    base = root/'docs/protium-fusion'
    sys.path.insert(0, str(base))
    import pulse_model as pm
    path = base/'membrane-receiver/stopping.py'
    spec = importlib.util.spec_from_file_location('selected_stopping', path)
    st = importlib.util.module_from_spec(spec); spec.loader.exec_module(st)
    data = base/'component-design/collection-stopping-data.json'
    st.DATA += json.loads(data.read_text())['selected_data']
    st.EL['Al'] = (13, 26.9815385, 166.)
    st.MAT.update({'aluminium': {'rho': 2300., 'parts': {'Al': 1}, 'db': 'Al', 'I': 166.},
                   'alumina': {'rho': 3900., 'parts': {'Al': 2, 'O': 3}, 'db': 'Al2O3', 'I': 145.2}})
    W = np.zeros((pm.NN, pm.NN)); channels = []
    for j, r in enumerate(pm.RATES):
        if r.weak_type:
            assert len(r.reactants) == 1 and r.dens_exp == 0
            parent = pm.NUCS.index(r.reactants[0]); rate = r.eval(1e8)
            W[:, parent] += pm.STOICH[:, j]*rate
            channels.append({'name': r.fname, 'rate_s_inverse': rate})
    eig, V = np.linalg.eig(W); invV = np.linalg.inv(V)
    assert np.max(np.abs(eig.imag)) < 1e-14
    eig, V, invV = eig.real, V.real, invV.real
    def decay(y, time):
        return np.maximum((np.exp(np.asarray(time)[..., None]*eig)*(invV@y))@V.T, 0.)
    test = np.arange(1, pm.NN+1, dtype=float); test /= pm.AA@test
    for t in [1., 100., 1000.]:
        assert np.max(abs(decay(test, t)-expm(W*t)@test)) < 1e-13
    hashes = {str(p.relative_to(root)): sha(p) for p in [path, data,
        base/'membrane-receiver/stopping-data.json', base/'pulse_model.py']}
    return pm, st, decay, channels, hashes

def spectrum(record, pm, decay, radius, order):
    d = json.loads(record.read_text()); trace = Path(d['trace_path']); z = np.load(trace)
    nz, nv = int(z['zones']), int(z['variables_per_zone'])
    last = z['state'].reshape(nz, nv, -1)[:, :, -1]
    names = d['state_variables']; assert names[:pm.NN] == pm.SPECIES
    Y = np.maximum(last[:, :pm.NN], 0.); m = z['cell_mass_kg']
    v = np.r_[0., last[:, names.index('v/cs0')]*float(z['velocity_unit_cm_s'])/100]
    r = np.r_[0., last[:, names.index('r/R0')]*d['R0_m']]
    assert np.all(np.diff(v)>0) and np.all(np.diff(r)>0)
    x, w = np.polynomial.legendre.leggauss(order)
    speed = .5*(v[1:, None]+v[:-1, None])+np.diff(v)[:, None]*x/2
    rr = .5*(r[1:, None]+r[:-1, None])+np.diff(r)[:, None]*x/2
    mass = m[:, None]*np.broadcast_to(w/2, speed.shape)
    proper = (radius-rr)/speed*np.sqrt(1-(speed/pm.C)**2)
    wallY = np.stack([decay(Y[i], proper[i]) for i in range(nz)])
    n = (mass[:, :, None]/pm.MU*wallY).reshape(-1, pm.NN)
    eu = (.5*pm.MU*speed**2/pm.MEV).ravel()
    ncat = d['mass_kg']/pm.MU*sum(y for s,y in d['initial_Y'].items() if s!='p')
    k = float(np.sum(n*pm.AA*eu[:, None])*pm.MEV)
    smooth = float(np.sum(m/np.diff(v)*np.diff(v**3)/6))
    checks = {'baryon_mass_relative': float(np.sum(n*pm.AA)*pm.MU/d['mass_kg']-1),
              'smooth_kinetic_relative': k/smooth-1}
    assert max(abs(x) for x in checks.values()) < 1e-9
    return d, n, eu, ncat, checks, {'record': sha(record), 'trace': sha(trace)}, v[-1]

def reduce(record, root, radius, pm, st, decay, table_cache):
    d,n,eu,ncat,checks,hashes,vmax = spectrum(record,pm,decay,radius,64)
    _,n128,eu128,_,_,_,_ = spectrum(record,pm,decay,radius,128)
    checks['64_vs_128_wall_species_L1_relative_nuclei'] = float(np.sum(abs(n.sum(axis=0)-n128.sum(axis=0)))/n.sum())
    area=4*math.pi*radius**2; Z=pm.ZZ; A=pm.AA; cno=(Z>=6)&(Z<=8); rows=[]
    for material in ['carbon','aluminium','alumina']:
        atom_mass=st.target(material)[0]*st.MU/st.target(material)[2]
        rho=st.MAT[material]['rho']
        for factor in [1., .5, 2.]:
            ranges=np.zeros_like(n); carbon_fusion=np.zeros_like(n)
            for j in range(pm.NN):
                if Z[j]==0 or n[:,j].sum()<n.sum()*1e-12: continue
                key=(material,int(Z[j]),int(A[j]),factor)
                if key not in table_cache: table_cache[key]=st.table(*key)
                eg,rg,_,se,sn=table_cache[key]
                assert np.all(se+sn>0) and eu.max()*A[j]<eg[-1]
                ranges[:,j]=np.interp(A[j]*eu,eg,rg)/100
                if material=='carbon' and Z[j]==6 and A[j] in [12,13]:
                    ec=eg*12/(12+A[j]); ew=np.maximum(ec,4.9)
                    sig=2.1*7.4**2/(2*ew)*np.logaddexp(0,2*np.pi*(ew-5.7)/2.1)/100
                    g=.9895*36*np.sqrt(12*A[j]/(12+A[j]))
                    sig*=np.where(ec<4.9,4.9/ec*np.exp(np.maximum(-700.,-g*(1/np.sqrt(ec)-1/np.sqrt(4.9)))),1.)
                    tau=cumulative_trapezoid(sig/(se+sn)/100/atom_mass*1e-28,eg,initial=0)
                    carbon_fusion[:,j]=-np.expm1(-np.interp(A[j]*eu,eg,tau))
            q=quantile(ranges[:,cno].ravel(),n[:,cno].ravel())
            if factor == 1.:
                r128=np.zeros_like(n128)
                for j in np.where(cno)[0]:
                    key=(material,int(Z[j]),int(A[j]),factor)
                    if key not in table_cache: table_cache[key]=st.table(*key)
                    eg,rg,_,_,_=table_cache[key]
                    r128[:,j]=np.interp(A[j]*eu128,eg,rg)/100
                q128=quantile(r128[:,cno].ravel(),n128[:,cno].ravel())
                checks['CNO_range_quantile_64_vs_128_max_relative_'+material]=max(abs(q[p]/q128[p]-1) for p in q)

            charged=Z>0; event=[]
            for sigma in [.1,1.,3.]:
                p=-np.expm1(-sigma*1e-28*ranges/atom_mass)
                event.append({'sigma_barn':sigma,'CNO_event_fraction':float(np.sum(n[:,cno]*p[:,cno])/n[:,cno].sum()),
                              'all_projectile_events_per_initial_catalyst':float(np.sum(n*p)/ncat)})
            coating=[]
            for thick in [10.,25.,50.,100.]:
                inside=ranges<=rho*thick*1e-6
                coating.append({'thickness_um':thick,'CNO_fraction_beyond_nominal_CSDA_range':float(np.sum(n[:,cno]*(~inside[:,cno]))/n[:,cno].sum()),
                                'all_charged_mass_fraction_beyond_nominal_CSDA_range':float(np.sum(n[:,charged]*A[charged]*(~inside[:,charged]))/np.sum(n[:,charged]*A[charged]))})
            rows.append({'material':material,'electronic_stopping_multiplier':factor,'density_kg_m3':rho,
                         'CNO_range_quantiles_kg_m2':q,'CNO_range_quantiles_um':{p:x/rho*1e6 for p,x in q.items()},
                         'maximum_packet_stopping_column_kg_m2':float(np.max(ranges)),
                         'constant_cross_section_event_sensitivities':event,'coatings':coating,
                         'carbon_fusion_fit_events_per_initial_catalyst':float(np.sum(n*carbon_fusion)/ncat) if material=='carbon' else None})
    return {'case':record.stem,'record_path':str(record),'radius_m':radius,'input_sha256':hashes,
            'hydrogen_per_initial_catalyst':d['mixture_inputs']['hydrogen_per_catalyst'],'initial_catalysts':ncat,
            'source_maximum_MeV_u':float(.5*pm.MU*vmax*vmax/pm.MEV),
            'kinetic_fluence_J_m2_smooth':float(np.sum(n*A*eu[:,None])*pm.MEV/area),
            'CNO_nuclei_at_wall_per_initial_catalyst':float(n[:,cno].sum()/ncat),
            'wall_nuclei_per_initial_catalyst':{s:float(n[:,j].sum()/ncat) for j,s in enumerate(pm.SPECIES) if n[:,j].sum()/ncat>1e-10},
            'rows':rows,'checks':checks}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--record',type=Path,action='append',required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--radius',type=float,default=1e9);a=ap.parse_args();root=a.root.resolve()
    pm,st,decay,channels,hashes=setup(root);cache={}
    cases=[reduce(p.resolve(),root,a.radius,pm,st,decay,cache) for p in a.record]
    out={'scope':'Species-resolved selected-mixture CSDA range and nuclear-event sensitivity. Coating penetration means redistribution into backing or another retained receiver reservoir, not escape from the closed station. No implantation-release yield or achieved catalyst return is inferred.',
         'primary_stopping_source':'https://www-nds.iaea.org/stopping/',
         'nuclear_C_C_source':'https://doi.org/10.1103/PhysRevC.85.014607',
         'nuclear_scope':'Only C12/C13-on-carbon fusion uses the illustrative measured-fit extension documented in component-design/collection.md. Constant cross sections are sensitivity functions, not cross-section predictions or terminal catalyst destruction. Recoverable isotope conversion remains in the material loop.',
         'range_scope':'CSDA path, not projected depth/straggling guarantee. Electronic multipliers 0.5 and 2 are sensitivity choices, not rigorous bounds or 80% credible intervals. Actual charge equilibration, damage and material phase need qualification.',
         'velocity_quadrature_order_per_cell':64,'finite_flight_scope':'Actual one-body network decay, parent velocity inherited by daughter; recoil and charged decay products separate. All baryon mass retained; no normalization to nodal kinetic energy.',
         'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'weak_channels':channels,'cases':cases}
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{'case':x['case'],'max_MeV_u':x['source_maximum_MeV_u'],'depth_um':[{r['material']:r['CNO_range_quantiles_um']} for r in x['rows'] if r['electronic_stopping_multiplier']==1.],'checks':x['checks']} for x in cases],indent=2))
if __name__=='__main__': main()
