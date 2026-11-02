"""Independent mathematical and archived-input checks for the diagnostic."""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import kve
from microphysics import *
ap=argparse.ArgumentParser()
ap.add_argument('--repo',default=str(Path(__file__).resolve().parents[2]))
ap.add_argument('--trace-dir',default='/tmp/protium-research')
ap.add_argument('--out',default=str(Path(__file__).parent))
a=ap.parse_args();ROOT=Path(a.repo);OUT=Path(a.out)
sys.path.insert(0,str(ROOT))
import pynucastro as pn
from pynucastro.eos.electron_eos import ElectronEOS
from pynucastro.screening.screen import f0
checks=[]
def check(name,error,tolerance,detail=None):
    checks.append({'name':name,'error':float(error),'tolerance':float(tolerance),
                   'pass':bool(error<=tolerance),'detail':detail})

# The ideal classical relativistic gas has an independent Bessel-function
# number integral. Keep pair production negligible for this control.
pref=(ME*C/HB)**3/np.pi**2
for theta in [.0001,.001,.01]:
    T=theta*ME*C*C/KB;eta0=-15.
    ne=pref*theta*kve(2,1/theta)*np.exp(eta0)
    eta,nm,np_=fermions(T,ne,order=128)
    check('Maxwell-Juttner number '+str(theta),abs(float(eta)-eta0),2e-7)

# Density and temperature span sampled pulse states. The independent source
# EOS integrates relativistic Fermi functions with a different variable/method.
eos=ElectronEOS();comp=pn.Composition(['p']);comp.set_all(1.)
points=[(1e7,1e4),(1e8,1e6),(3e8,5e7),(1e9,5e7),
        (1.5e9,1e8),(1.2e9,5e9),(1e8,1e9),(1e7,1e6)]
for T,rhoye in points:
    eta,nm,np_=fermions(T,rhoye/MU,order=128)
    finer=fermions(T,rhoye/MU,order=256)
    ele,pos=eos.pe_state(rhoye/1000,T,comp,compute_derivs=False)
    check('eta vs independent EOS '+str((T,rhoye)),abs(float(eta)-ele.eta),2e-5)
    check('Fermi quadrature refinement '+str((T,rhoye)),abs(float(eta-finer[0])),2e-5)
    check('charge neutrality '+str((T,rhoye)),abs(float((nm-np_)*MU/rhoye-1)),1e-9)

# Kirchhoff consistency integrates alpha_nu times the Planck function, rather
# than comparing two copies of a frequency-integrated coefficient.
for T in [1e7,1e8,1e9]:
    ne=1e32;z2ni=3e32
    def integrand(x):
        alpha=freefree(T,ne,z2ni,x=x)[0]
        nu=x*KB*T/H
        bnu=2*H*nu**3/C**2/np.expm1(x)
        return 4*np.pi*alpha*bnu*KB*T/H
    actual=quad(integrand,1e-8,60,epsabs=0,epsrel=2e-10)[0]
    expected=freefree(T,ne,z2ni)[1]
    check('Kirchhoff emissivity '+str(T),abs(actual/expected-1),5e-4)

# Published nondegenerate fit differs slightly from the all-eta fit; table
# accuracy is percent-level, not machine precision.
for kev in [2.,10.,50.,100.,300.]:
    T=kev*1000*E/KB;rho=1.;lep=1/MU
    actual,valid=poutanen(T,np.asarray(-10.),lep,rho)
    expected=SIGT*lep/rho/(1+(kev/41.5)**.9)
    check('Poutanen nondegenerate fit '+str(kev),abs(float(actual/expected-1)),.085,
          'Published all-eta fit error <=6.5%, separate nondegenerate fit <=2%; 8.5% comparison budget.')

# OCP energy is a derivative of the free energy; finite differences provide
# an independent check on the analytic derivative and thermodynamic scaling.
for g in np.geomspace(1e-3,100.,12):
    dh=1e-4
    numeric=(f0(g*np.exp(dh))-f0(g*np.exp(-dh)))/(2*dh)
    actual=ocp_internal(g)
    check('OCP free-energy derivative '+str(g),abs(actual/numeric-1),1e-6)
check('OCP weak-coupling limit',abs(float(ocp_internal(1e-8))/(-np.sqrt(3)/2*1e-12)-1),.001)

# Spitzer-Harm Lorentz coefficient reproduces the standard hydrogen result
# after inserting its published epsilon*delta_T instead of the Lorentz 0.4.
T=1e7;ne=1e25
_,kc,kd=conductivity_envelope(T,ne,ne,lnlambda=10.)
hydrogen=kc*(.419*.225)/.4
standard=1.84e-10*T**2.5/10
check('Spitzer-Harm hydrogen normalization',abs(hydrogen/standard-1),.005)
for rho,T in [(1e8,1e9),(1e10,1e9)]:
    ne=rho*.5/MU;z2=6*ne
    _,_,kd=conductivity_envelope(T,ne,z2)
    pf=HB*(3*np.pi**2*ne)**(1/3);mf=np.sqrt(ME**2+(pf/C)**2)
    # Algebraically eliminate nu and vf from Cassisi et al eqs2 and8.
    independent=np.pi*KB**2*T*ne*pf**3/(12*EC**2*z2*mf**2)
    check('Degenerate conductivity algebra '+str(rho),abs(float(kd/independent)-1),1e-12)

# NRL Formulary Maxwellian equal-temperature electron-ion exchange normalization.
TeV=100.;T=TeV*E/KB;ni=1e25;mi=MU
rate=8*np.sqrt(2*np.pi)*ni*EC**2/(3*ME*mi)*(KB*T/ME+KB*T/mi)**(-1.5)
nrl=3.2e-9*(ni/1e6)*TeV**(-1.5)
check('NRL electron-ion exchange normalization',abs(rate/nrl-1),.02)

archive=json.loads((ROOT/'compression-traces.json').read_text())['traces']
audit=json.loads((OUT/'pulse-transport-audit.json').read_text())
for row in audit['rows']:
    stem=row['stem'];raw=Path(a.trace_dir)/(stem+'.npz')
    actual=hashlib.sha256(raw.read_bytes()).hexdigest()
    check('trace hash '+stem,0 if actual==archive[stem]['source_npz_sha256']==row['trace_sha256'] else 1,0)
    q=np.load(OUT/(stem+'-diagnostics.npz'))
    check('shell mass positivity '+stem,0 if np.all(q['mass_kg']>0) else 1,0)
    check('time ordering '+stem,0 if np.all(np.diff(q['time_s'])>0) else 1,0)
    check('radius ordering '+stem,0 if np.all(np.diff(q['radius_m'],axis=0)>0) else 1,0)
    check('nonnegative nuclear statistical weights '+stem,0 if np.all(q['nuclear_weight_J']>=0) else 1,0)
    check('nuclear weighting sum '+stem,abs(q['nuclear_weight_J'].sum()/row['positive_sampled_nuclear_increments_J']-1),1e-12)
    # Test a deterministic grid of real sampled states against 256-point rule.
    pick=np.linspace(0,q['T_K'].size-1,32,dtype=int)
    tt=q['T_K'].ravel()[pick];rr=q['rho_kg_m3'].ravel()[pick]
    yy=q['normalized_molar']
    import pulse_model as pm
    ye=(yy*pm.ZZ).sum(axis=-1).ravel()[pick]
    eta=fermions(tt,rr*ye/MU,order=256)[0]
    check('real-state Fermi quadrature '+stem,np.max(np.abs(eta-q['eta'].ravel()[pick])),2e-5)

out={'checks':checks,'passed':sum(x['pass'] for x in checks),'total':len(checks),
     'all_pass':all(x['pass'] for x in checks),
     'scope':'Mathematical verification and archived input identity, not experimental validation.'}
(OUT/'pulse-transport-validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
for c in checks:
    if not c['pass']:print(c)
if not out['all_pass']:raise SystemExit(1)
