"""Local counterstream/guide linear dynamics and conditional feedback screen.
This is not a discrete-carrier, field solver, or certified control design.
Fourier modes retain equal-and-opposite guide/stream actuator reaction.
"""
from pathlib import Path
import json, math, hashlib
import numpy as np
from scipy.linalg import solve_continuous_are
P=Path(__file__).parent
paths={k:P/v for k,v in {'guide':'guide-component.json','topology':'support-topology.json','sector':'sector-account.json'}.items()}
raw={k:json.loads(v.read_text()) for k,v in paths.items()}
omega=raw['guide']['confinement']['transverse_omega_rad_s']
v=raw['guide']['parameters']['equatorial_speed_m_s']
sigma=raw['sector']['support']['total_areal_mass_kg_m2']
sigmastream=max(raw['topology']['geometry']['axis_surface_stream_kg_m2'].values())
epsilon=sigmastream/sigma
T=1e4 # N/m declared weak mesh pretension, not its in-plane elastic modulus
margin=raw['guide']['confinement']['available_transverse_margin_m']
yref=.001
# Dimensionless time omega*t, displacement yref, velocity omega*yref.
def matrices(wavelength,kappa_scale=1.,stream_mass_scale=1.):
    k=2*math.pi/wavelength; p=k*v/omega; c=kappa_scale**2
    e=epsilon*stream_mass_scale; wg2=T/sigma*k*k/omega**2
    A=np.zeros((6,6),complex); A[0,1]=A[2,3]=A[4,5]=1
    A[1,0]=A[3,2]=p*p-c
    A[1,1]=-2j*p; A[3,3]=2j*p
    A[1,4]=A[3,4]=c
    A[5,0]=A[5,2]=e*c/2; A[5,4]=-wg2-e*c
    B=np.zeros((6,2),complex); B[1,0]=B[3,1]=1; B[5,:]=-e/2
    return A,B,{'p':p,'c':c,'epsilon':e,'wg2':wg2}

def growth(A):
    z=np.linalg.eigvals(A); return float(max(z.real)*omega)

def closed(A,B,K,tau,gain=1.):
    tp=tau*omega
    return np.block([[A,B],[-gain*K/tp,-np.eye(2)/tp]])

critical=2*math.pi*v/omega
wavelengths=np.unique(np.r_[np.geomspace(1e3,1e6,61),critical*np.linspace(.9,1.1,41)])
open_rows=[]; Ks={}; riccati_max=0.; poly_max=0.; reaction_max=0.
for w in wavelengths:
    A,B,p=matrices(w)
    ev=np.linalg.eigvals(A)
    # Independent polynomial in x=s²: D(x)(x+wg²+eps*c)-eps*c²(x+c-b).
    b=p['p']**2; c=p['c']; e=p['epsilon']; aa=c-b
    coeff=np.polymul([1,2*(c+b),aa*aa],[1,p['wg2']+e*c])
    coeff[-2:]-=e*c*c*np.array([1,aa])
    for s in ev:
        x=s*s; scale=sum(abs(co)*max(1.,abs(x))**(3-i) for i,co in enumerate(coeff))
        poly_max=max(poly_max,float(abs(np.polyval(coeff,x))/scale))
    reaction_max=max(reaction_max,float(np.max(abs(np.array([0,e/2,0,e/2,0,1])@B))))
    S=solve_continuous_are(A,B,np.eye(6),np.eye(2)); K=B.conj().T@S; Ks[float(w)]=K
    residual=A.conj().T@S+S@A-S@B@B.conj().T@S+np.eye(6)
    riccati_max=max(riccati_max,float(np.linalg.norm(residual)/(1+np.linalg.norm(S@B@B.conj().T@S))))
    open_rows.append({'wavelength_m':float(w),'growth_per_s':float(max(ev.real)*omega),'nominal_instantaneous_feedback_growth_per_s':growth(A-B@K)})

robust=[]
for tau in [.001,.01,.1,1.]:
    worst={'growth_per_s':-math.inf}; stable_count=0; total=0; band_total=0; band_stable=0; band_worst=-math.inf
    for ks in [.8,1.,1.2]:
      for es in [.5,1.,2.]:
       for gain in [.8,1.,1.2]:
        vals=[]
        for w in wavelengths:
            A,B,_=matrices(w,ks,es); gr=growth(closed(A,B,Ks[float(w)],tau,gain)); total+=1
            if gr<0: stable_count+=1
            if w>=2e4:
                band_total+=1; band_stable+=int(gr<0); band_worst=max(band_worst,gr)
            if gr>worst['growth_per_s']:
                worst={'growth_per_s':gr,'wavelength_m':float(w),'kappa_scale':ks,'stream_mass_scale':es,'actuator_gain_scale':gain}
        
    robust.append({'first_order_actuator_time_s':tau,'tested_combinations':total,'stable_combinations':stable_count,'all_tested_stable':stable_count==total,'worst':worst,'active_band_20_to_1000km':{'tested':band_total,'stable':band_stable,'worst_growth_per_s':band_worst}})

# Position-only axis cancellation is not equivalent to closed-loop full-state
# control. Record a narrow counterexample rather than crediting perfect tracking.
def axis_servo(w,tau=.01,position_gain=1.,velocity_lead=0.):
    A,_,_=matrices(w); z=np.zeros((7,7),complex);z[:6,:6]=A
    z[1,6]=z[3,6]=1.;z[5,6]=-epsilon
    tp=tau*omega;z[6,4]=-position_gain/tp;z[6,5]=-velocity_lead;z[6,6]=-1/tp
    return z
servo=[]
for gain,lead in [(1.,0.),(1.,1.),(.99,1.),(1.01,1.)]:
    best=max((growth(axis_servo(w,position_gain=gain,velocity_lead=lead)),float(w)) for w in wavelengths)
    servo.append({'position_gain':gain,'velocity_lead':lead,'actuator_time_s':.01,'max_growth_per_s':best[0],'wavelength_m':best[1]})

# Modal transient requirements. Input yref is an amplitude normalization, not
# a measured disturbance. Real physical fields combine k and -k conjugates.
transient=[]
for w in [1e3,1e4,4e4,critical,5e4,1e5,1e6]:
    A,B,_=matrices(w); S=solve_continuous_are(A,B,np.eye(6),np.eye(2));K=B.conj().T@S
    C=closed(A,B,K,.01)
    ev,V=np.linalg.eig(C); Vi=np.linalg.inv(V)
    times=np.r_[0.,np.geomspace(1e-5,1e5,700)]
    oe,ov=np.linalg.eig(A);j=np.argmax(oe.real);unstable=ov[:,j];unstable/=max(abs(unstable[[0,2,4]]))
    seeds={'guide_offset_1mm':np.array([0,0,0,0,1,0],complex),'stream_plus_offset_1mm':np.array([1,0,0,0,0,0],complex),'open_loop_fastest_mode_positions_at_most_1mm':unstable}
    for name,seed in seeds.items():
        x0=np.r_[seed,[0.,0.]]
        states=V@((Vi@x0)[:,None]*np.exp(ev[:,None]*omega*times[None,:]))
        initial_res=float(np.max(abs(states[:,0]-x0)))
        assert initial_res<1e-7
        acc=float(np.max(abs(states[6:]))*omega**2*yref)
        axis_shift=acc/omega**2
        rel=float(max(np.max(abs(states[0]-states[4])),np.max(abs(states[2]-states[4])))*yref)
        transient.append({'wavelength_m':float(w),'seed':name,'maximum_stream_relative_to_guide_m':rel,'maximum_guide_displacement_m':float(np.max(abs(states[4]))*yref),'maximum_command_acceleration_m_s2':acc,'equivalent_quadrupole_axis_trim_m':axis_shift,'equivalent_dipole_trim_T':axis_shift*raw['guide']['parameters']['gradient_T_m'],'within_position_and_trim_margin_separately':rel<margin and axis_shift<margin,'maximum_eigenvector_reconstruction_residual':initial_res})

fast=max(open_rows,key=lambda r:r['growth_per_s'])
result={'scope':'One local straight rail family with paired equal counterstreams, flexible guide coupled to full receiver tributary mass, weak mesh tension. Linear Fourier screen; full-state feedback is mode-dependent and presumes sensors/communication. No certified three-family controller, finite-carrier result, saturation proof or failure recovery.',
'input_sha256':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
'parameters':{'stream_to_receiver_mass_ratio':epsilon,'heaviest_family_kg_m2':sigmastream,'receiver_total_kg_m2':sigma,'stream_m_s':v,'quadrupole_transverse_omega_rad_s':omega,'mesh_pretension_N_m':T,'critical_wavelength_m':critical,'reference_displacement_m':yref,'available_transverse_margin_m':margin},
'open_loop_fastest':{**fast,'efold_s':1/fast['growth_per_s']},'open_loop_rows':open_rows,'position_only_axis_servo_counterexamples':servo,'full_state_feedback_robustness':robust,'modal_transients':transient,'transient_interpretation':'The 1 km full-state controller unnecessarily damps already stable short waves and exceeds the reference trim in two 1 mm seeds. The 40–1000 km samples all fit separate position and trim margins. Do not convert an arbitrary short-wave LQR choice into a necessary alignment or latency requirement.',
'checks':{'independent_polynomial_scaled_residual_max':poly_max,'actuator_reaction_residual_max':reaction_max,'CARE_relative_residual_max':riccati_max},
'sources':{'fluid_conveying_pipe_theory_experiment':'https://journals.sagepub.com/doi/10.1243/JMES_JOUR_1970_012_017_02','electromagnetic_control_computation':'https://doi.org/10.1016/j.jsv.2018.06.045','semi_active_experimental_control':'https://www.ippt.pan.pl/repository/open/o6415.pdf'},
'next_step':'Implement a real-space finite rail/carrier/actuator model with shared cross-braced deck, measurement noise, latency, field saturation, and discrete faults; compare hardware and AC losses with the existing service/auxiliary allowance. The full-state screen establishes a possible stabilization direction, not its installed cost.'}
assert poly_max<1e-10 and reaction_max<1e-14 and riccati_max<1e-5
Path(__file__).with_name('guide-local-control.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'open_loop_fastest':result['open_loop_fastest'],'robustness':robust,'checks':result['checks'],'max_transient_acceleration':max(x['maximum_command_acceleration_m_s2'] for x in transient),'max_relative_displacement':max(x['maximum_stream_relative_to_guide_m'] for x in transient),'transients_within_margin':sum(x['within_position_and_trim_margin_separately'] for x in transient)},indent=2))
