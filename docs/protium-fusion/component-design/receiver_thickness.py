"""Receiver thickness: finite conduction, finite heat storage and local plate loads.
Actual selected pulse mixture, with all source heat retained for the thermal screen.
Photon-containment fractions are a separate study; no capture credit is inferred.
"""
from pathlib import Path
import json, math, hashlib, sys
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erfc, zeta
P=Path(__file__).parent;ROOT=P.parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(P))
from fast_eos import FastEOS
import pulse_model as pm
import material_cycle as mc
src=Path('/home/tkwa/src/protium-fusion/output/mixture-runs/plant-scan-003.json')
packet=json.loads(src.read_text())['top_detailed'][0];plant=packet['plant'];cycle=packet['cycle'];R=plant['radius_m'];AREA=plant['area_m2']
RHO=2000.;DT=.25;CUT=20000.;GUARD=86400.;C=299792458.;YEAR=365.25*86400
stages={};hashes={str(src):mc.hashfile(src),str(Path(mc.__file__)):mc.hashfile(mc.__file__)}
for dep in ['fast_eos.py','electron-eos.npz','pulse_model.py']:
    hashes[str(ROOT/dep)]=mc.hashfile(ROOT/dep)
rates={'n13':.001159113256783153,'o14':.009817978825996571,'o15':.005681556796002065}
daughters={'n13':'c13','o14':'n14','o15':'n15'}
CAT={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
eos=FastEOS();edges=np.arange(round(CUT/DT)+1)*DT;times=edges[1:]

def pulse_bins(s,ee,delay_bound=False):
    q=np.zeros(len(ee)-1);qu=q.copy();qp=q.copy();dt=np.diff(ee)
    for j in range(len(s['m'])):
        vl=np.maximum(s['v'][j],s['b'][j]/np.maximum(ee[1:]-s['off'][j],1e-100))
        vh=np.minimum(s['v'][j+1],s['b'][j]/np.maximum(ee[:-1]-s['off'][j],1e-100))
        ok=(vh>vl)&(ee[1:]>s['off'][j]);dv=np.where(ok,vh-vl,0.)
        q+=np.where(ok,s['slope'][j]*(vh**3-vl**3)/(6*AREA),0.)/dt
        qu+=s['slope'][j]*s['u'][j]*dv/(AREA*dt)
        qp+=np.where(ok,s['slope'][j]*(vh**2-vl**2)/(2*AREA),0.)/dt
    qdec=np.zeros_like(q)
    lag=2*R/C if delay_bound else 0.
    for n,E in s['decay_fluence'].items():
        rate=rates[n]
        # A conservative exponential envelope over a bounded photon arrival lag;
        # used for temperature only. Actual pulse-energy accounting remains E.
        left=ee[:-1];right=ee[1:]
        integ=rate*np.maximum(0.,np.minimum(right,lag)-left)
        aa=np.maximum(left,lag);bb=np.maximum(right,lag)
        integ+=np.exp(-rate*(aa-lag))-np.exp(-rate*(bb-lag))
        qdec+=E*integ/dt
    return q,qu,qdec,qp

for name,case in plant['selected_cases'].items():
    cachepath=ROOT/'mixture-design/receiver-results'/f'{case}.receiver.json';cache=json.loads(cachepath.read_text())
    assert math.isclose(cache['radius_m'],1e9,rel_tol=1e-12) and cache['completed_requested_duration'] and cache['usable_ballistic']
    rp=Path(cache['source_paths']['record']);tp=Path(cache['source_paths']['trace']);d=json.loads(rp.read_text());tr=np.load(tp)
    for path in [cachepath,rp,tp]:hashes[str(path)]=mc.hashfile(path)
    assert cache['input_sha256'][str(rp)]==hashes[str(rp)]
    nz=int(tr['zones']);nv=int(tr['variables_per_zone']);st=tr['state'].reshape(nz,nv,-1);names=d['state_variables'];nn=names.index('T/T0')
    m=tr['cell_mass_kg'];ro=np.r_[0.,st[:,names.index('r/R0'),-1]*d['R0_m']]
    v=np.r_[0.,st[:,names.index('v/cs0'),-1]*float(tr['velocity_unit_cm_s'])/100]
    Y=np.maximum(st[:,:nn,-1],0.);bar=Y@pm.AA;rho=m/(4*np.pi/3*np.diff(ro**3))
    u=eos.state(rho/1000,st[:,names.index('T/T0'),-1]*d['T0_K'],(Y@pm.ZZ)/bar,Y.sum(axis=1)/bar).e*1e-4
    drdv=np.diff(ro)/np.diff(v);b=R-ro[:-1]+drdv*v[:-1];off=float(tr['time_s'][-1])-drdv-R/C
    assert np.all(off<0) and np.all(u>0)
    stage=mc.load_stage(rp);decay={n:(mc.nuc(n).mass-mc.nuc(daughters[n]).mass-mc.DEFAULT_COLD_NU[n])*stage['hot_output'].get(n,0)*stage['attempts_per_pulse']*mc.MEV/AREA for n in rates}
    decsum=sum(decay.values());allcold=(stage['Q_cold_MeV']-stage['nu_cold_MeV'])*stage['attempts_per_pulse']*mc.MEV/AREA
    assert abs(decsum-allcold)<1e-6*max(1.,decsum)
    K=float(np.sum(m/np.diff(v)*np.diff(v**3)/6));U=float(m@u);prompt=float(d['radiation_escaped_J'])/AREA
    s={'m':m,'v':v,'b':b,'off':off,'slope':m/np.diff(v),'u':u,'decay_fluence':decay,'fluence':(K+U)/AREA+prompt+decsum,
       'K_fluence':K/AREA,'K_nodal_fluence':float(d['final_kinetic_J'])/AREA,'U_fluence':U/AREA,'prompt_fluence':prompt,'period':plant['pulse_interval_s'][CAT[name]],
       'prompt_peak_reference_K':cache['prompt']['radial_ray_arrival_surface_screen']['peak_surface_rise_K']*(cache['radius_m']/R)**2,
       'prompt_duration_s':cache['prompt']['radial_ray_arrival_surface_screen']['duration_s'],
       'minimum_slot':cache['timing_at_radius']['mass999']['occupied_s']*R/cache['radius_m'],
       'impulse':cache['ballistic']['momentum_fluence_Pa_s']*(cache['radius_m']/R)**2,
       'peak_pressure':cache['ballistic']['peak_momentum_flux_Pa']*(cache['radius_m']/R)**3}
    q,qu,qd,qp=pulse_bins(s,edges,True);s['q']=q+qu+qd
    s['actual_early_energy_fraction']=float(sum(sum(z) for z in pulse_bins(s,edges,False)[:3])*DT/(s['fluence']-prompt))
    stages[name]=s

# Exact step response at the front: front flux, either fixed-temperature back
# (D), or insulated back with the uniform bulk mode removed (N_nonuniform).
def front_step(t,L,k,cp,boundary):
    t=np.asarray(t);alpha=k/(RHO*cp);Carea=RHO*cp*L;out=np.zeros_like(t,dtype=float)
    short=(t>0)&(alpha*t/L**2<.2);ts=t[short]
    val=2*np.sqrt(alpha*ts)/(k*np.sqrt(np.pi))
    for j in range(1,8):
        zz=j*L/np.sqrt(alpha*ts)
        val+=4*np.sqrt(alpha)/(k*np.sqrt(np.pi))*(1 if boundary=='N' else (-1)**j)*(np.sqrt(ts)*np.exp(-zz*zz)-np.sqrt(np.pi)*j*L/np.sqrt(alpha)*erfc(zz))
    if boundary=='N':val-=ts/Carea
    out[short]=val
    long=(t>0)&~short;tt=t[long]
    n=np.arange(1,65) if boundary=='N' else np.arange(64)+.5
    lam=alpha*(n*np.pi/L)**2
    lim=L/(3*k) if boundary=='N' else L/k
    out[long]=lim-2/Carea*np.sum(np.exp(-lam[:,None]*tt)/lam[:,None],axis=0)
    return out

def history(L,k,cp,period,cut,boundary):
    alpha=k/(RHO*cp);n=np.arange(1,1001) if boundary=='N' else np.arange(1000)+.5;lam=alpha*(n*np.pi/L)**2
    return float(2/(RHO*cp*L)*np.sum(np.exp(-lam*(period-cut))/(-np.expm1(-lam*period))))

def tail_bound(s,period,cut):
    ok=s['v'][:-1]<=s['b']/(cut-s['off'])
    # Actual advected U has a t^-2 envelope; it must not be assigned K's t^-4.
    c4=float(np.sum(.5*s['slope'][ok]*s['b'][ok]**3/AREA))
    c2=float(np.sum(s['slope'][ok]*s['u'][ok]*s['b'][ok]/AREA))
    material=c4/period**4*zeta(4,cut/period)+c2/period**2*zeta(2,cut/period)
    decay=sum(E*rates[n]*math.exp(-rates[n]*max(0.,cut-2*R/C))/(-math.expm1(-rates[n]*period)) for n,E in s['decay_fluence'].items())
    return float(material+decay)

def back_pdf(t,L,k,cp):
    t=np.asarray(t);alpha=k/(RHO*cp);tau=alpha*t/L**2;out=np.zeros_like(t)
    small=(tau>0)&(tau<.2);z=tau[small];val=np.zeros_like(z)
    for n in range(8):val+=(-1)**n*(2*n+1)*np.exp(-(2*n+1)**2/(4*z))/(np.sqrt(np.pi)*z**1.5)
    out[small]=val
    large=tau>=.2;z=tau[large];n=np.arange(64)+.5
    out[large]=2*np.pi*np.sum((n*(-1)**np.arange(64))[:,None]*np.exp(-(n[:,None]*np.pi)**2*z),axis=0)
    return out*alpha/L**2

big=[n for n in stages if n!='B1'];qmean=sum(s['fluence']/s['period'] for s in stages.values());big_rate=sum(1/stages[n]['period'] for n in big)
assert big_rate*GUARD<1
rows=[];checks=[]
for k,cp in [(100.,1000.),(100.,1800.),(50.,1800.),(25.,1800.)]:
  for sigma in [1000.,2000.,5000.]:
    L=sigma/RHO;alpha=k/(RHO*cp);dstep=front_step(edges,L,k,cp,'D');nstep=front_step(edges,L,k,cp,'N');gD=np.diff(dstep);gN=np.diff(nstep)
    assert min(gD)>-1e-12 and min(gN)>-1e-12
    peaks={};back={}
    for name,s in stages.items():
      flash=s['prompt_peak_reference_K']*math.sqrt(100*1000/(k*cp))
      # Prompt is over before the first slow bin. Concentrate it at its final
      # instant for a conservative tail; never add two separated peak values.
      age=np.maximum(times-s['prompt_duration_s'],1e-10)
      ptail=s['prompt_fluence']/np.sqrt(np.pi*k*RHO*cp*age)
      slowD=fftconvolve(s['q'],gD)[:len(times)]+ptail
      slowN=fftconvolve(s['q'],gN)[:len(times)]+ptail
      flash+=float(s['q'][0])*2*math.sqrt(s['prompt_duration_s'])/math.sqrt(math.pi*k*RHO*cp)
      peaks[name]={'fixed_back_isolated_front_rise_bound_K':max(flash,float(max(slowD))),'finite_bulk_nonuniform_isolated_front_bound_K':max(flash,float(max(slowN))),'flash_bound_K':flash,'slow_material_plus_decay_peak_D_K':float(max(slowD)),'slow_peak_time_s':float(times[np.argmax(slowD)])}
      # Coarse5s energy bins are sufficient for the >=850s back-face response.
      ee=np.arange(0,CUT+5,5.);qparts=pulse_bins(s,ee,False);energy=sum(qparts[:3])*5;mid=(ee[:-1]+ee[1:])/2
      tt=np.geomspace(50.,max(1e5,4*L*L/alpha),350)
      vals=np.array([float(back_pdf(np.array([t]),L,k,cp)[0]*s['prompt_fluence']+np.dot(energy,back_pdf(t-mid,L,k,cp))) for t in tt])
      j=int(np.argmax(vals));back[name]={'isolated_fixed_back_peak_W_m2':float(vals[j]),'time_s':float(tt[j]),'ratio_to_station_mean_heat':float(vals[j]/qmean),'early_fluence_fraction':s['actual_early_energy_fraction']}
    smallP=stages['B1']['minimum_slot'];smallcut=min(CUT,smallP/2)
    residual={}
    for boundary in ['D','N']:
      oldbig=max(stages[n]['fluence'] for n in big)*history(L,k,cp,GUARD,CUT,boundary)
      oldsmall=stages['B1']['fluence']*history(L,k,cp,smallP,smallcut,boundary)
      resistance=L/k/(3 if boundary=='N' else 1)
      tail=(sum(tail_bound(stages[n],GUARD,CUT) for n in big)+tail_bound(stages['B1'],smallP,smallcut))*resistance
      residual[boundary]={'old_large_K':oldbig,'old_B1_K':oldsmall,'late_tails_K':tail,'total_K':oldbig+oldsmall+tail}
    buffer_span=sum(s['fluence'] for s in stages.values())/(sigma*cp)
    boundD=1000+max(peaks[n]['fixed_back_isolated_front_rise_bound_K'] for n in big)+peaks['B1']['fixed_back_isolated_front_rise_bound_K']+residual['D']['total_K']
    boundN=1000+buffer_span+max(peaks[n]['finite_bulk_nonuniform_isolated_front_bound_K'] for n in big)+peaks['B1']['finite_bulk_nonuniform_isolated_front_bound_K']+residual['N']['total_K']+qmean*L/(6*k)
    rows.append({'sigma_kg_m2':sigma,'thickness_m':L,'k_W_m_K':k,'cp_J_kg_K':cp,'diffusivity_m2_s':alpha,'through_diffusion_s':L*L/alpha,'fixed_back_slowest_decay_s':4*L*L/(alpha*np.pi*np.pi),'finite_bulk_nonuniform_slowest_decay_s':L*L/(alpha*np.pi*np.pi),
       'mean_front_back_gradient_K':qmean*L/k,'all_stage_bulk_swing_upper_K_at_constant_mean_extraction':buffer_span,'per_stage_uniform_energy_rise_K':{n:s['fluence']/(sigma*cp) for n,s in stages.items()},
       'isolated_peaks':peaks,'fixed_back_peak_heat_pickup':back,'background_bounds':residual,'fixed_1000K_back_front_temperature_bound_K':boundD,'finite_buffer_front_temperature_bound_K_if_each_bulk_mode_minimum_1000K':boundN})
    # Independent integral of positive exponential response gives its thermal resistance.
    tt=np.array([.001*L*L/alpha,1e4*L*L/alpha]);fs=front_step(tt,L,k,cp,'D');ns=front_step(tt,L,k,cp,'N')
    checks.append({'sigma':sigma,'k':k,'cp':cp,'D_long_resistance_relative':float(fs[-1]/(L/k)-1),'N_long_resistance_relative':float(ns[-1]/(L/(3*k))-1),'D_early_halfspace_relative':float(fs[0]/(2*math.sqrt(alpha*tt[0])/(k*math.sqrt(math.pi)))-1)})

# Local unsupported tile bending screen; the cold frame/guide supports the
# entire shell globally. These are small-load elastic requirements, not a
# passive-shell or thermoshock qualification.
modal_cut=501
modal_sum=sum(1/(m*n*(m*m+n*n)) for m in range(1,modal_cut+1,2) for n in range(1,modal_cut+1,2))
modal_tail_bound=(1.5+math.log(modal_cut))/modal_cut**2
modal_factor=2*(modal_sum+modal_tail_bound)
plate=[];E=1e10;nu=.2;imp=max(s['impulse'] for s in stages.values());ppeak=max(s['peak_pressure'] for s in stages.values())
for sigma in [1000.,2000.,5000.]:
 t=sigma/RHO;D=E*t**3/(12*(1-nu*nu))
 for span in [100.,250.,1000.]:
  om=2*math.pi**2/span**2*math.sqrt(D/sigma);amp=16/math.pi**2*imp/(sigma*om)
  stress=E*t*math.pi**2*amp/(2*(1-nu)*span**2)
  plate.append({'sigma_kg_m2':sigma,'tile_span_m':span,'flexural_rigidity_N_m':D,'first_mode_period_s':2*math.pi/om,'uniform_peak_static_centre_deflection_m':.00406*ppeak*span**4/D,'single_mode_impulsive_displacement_upper_m':amp,'all_mode_absolute_impulsive_displacement_upper_m':amp*modal_factor,'single_mode_impulsive_bending_stress_Pa':stress,'load_scope':'Uniform normal pulse, simply-supported square plate. The impulsive bound ignores finite pulse-duration reduction; higher modes, joints, thermal stress and the separate mean self-gravity load remain outside it.'})

out={'scope':__doc__,'input_sha256':hashes,'code_sha256':mc.hashfile(__file__),'reference':{'radius_m':R,'helium_per_s':plant['helium_per_s'],'fixed_original_plant_mass_kg':plant['mass_kg'],'fixed_original_receiver_kg_m2':5000.,'selected_cases':plant['selected_cases'],'total_received_heat_mean_W_m2_assuming_full_retention':qmean,'large_pulse_guard_s':GUARD,'large_guard_capacity_used':big_rate*GUARD,'early_history_cut_s':CUT},
 'pulse_energy_fluence_J_m2':{n:{'kinetic':s['K_fluence'],'exact_nodal_kinetic_separate':s['K_nodal_fluence'],'advected_internal':s['U_fluence'],'prompt':s['prompt_fluence'],'delayed':s['decay_fluence'],'total':s['fluence'],'interval_s':s['period'],'early_slow_energy_fraction':s['actual_early_energy_fraction']} for n,s in stages.items()},
 'thermal_rows':rows,'local_plate_rows':plate,'plate_absolute_modal_factor':modal_factor,'checks':checks,'source_energy_convention':{'mass_conserving_smooth_mean_W_m2':qmean,'exact_plant_ledger_full_source_W_m2':plant['captured_pulse_W']/.99/AREA,'relative_difference_from_velocity_interpolation':qmean/(plant['captured_pulse_W']/.99/AREA)-1},
 'limitations':['Thermal calculation assigns full retained energy to the front surface; actual gamma containment/depth is separate and may reduce these temperature bounds.',
 'Constant thermal properties and linear slab geometry; no vapor, cracking, high-temperature material degradation, 2D joints or radiative recirculation model.',
 'Fixed-back calculation requires a heat sink able to accept its calculated transient flux. Mean-rated hardware is not credited with that peak for free.',
 'Finite-buffer case extracts constant mean heat. Each stage bulk-energy mode is shifted so its minimum is1000K; their peak-to-peak sum bounds the variable bulk temperature without an infinite heat bath. It is a conservative all-phase bound, not a solved dispatch schedule.',
 'The incoming prompt flash is evaluated separately from late material/decay heat and only its decaying tail is added at later times.',
 'Delayed heating uses actual parent counts and existing non-neutrino decay ledger. The timing envelope permits a2R/c photon travel spread; charged-product transport is a separate source-history assumption.',
 'No W/kg increase is assigned to a thinner receiver until actual photon containment, temperature-dependent engine performance, supports and process interfaces are integrated.'],
 'sources':{'graphite_1000K_heat_capacity':'https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=921062','high_temperature_graphite_grade_measurements':'https://doi.org/10.1016/j.jnucmat.2020.152176','room_temperature_graphite_modulus_scale':'https://www.sglcarbon.com/pdf/SGL-Datasheet-SIGRAFINE-R6650-EN.pdf'}}
Path(__file__).with_name('receiver-thickness.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'mean_heat':qmean,'pulse_fluence':out['pulse_energy_fluence_J_m2'],'thermal_summary':[{k:r[k] for k in ['sigma_kg_m2','k_W_m_K','cp_J_kg_K','all_stage_bulk_swing_upper_K_at_constant_mean_extraction','fixed_1000K_back_front_temperature_bound_K','finite_buffer_front_temperature_bound_K_if_each_bulk_mode_minimum_1000K']} for r in rows]},indent=2))
