"""Conditional He trap kinetics along archived chronological receiver T(t).
Run with numpy/scipy: python thermal_trapping.py
No measured trapped fraction, trap-energy distribution or repair work is inferred.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
from scipy.special import erfc
P=Path(__file__).resolve().parent
KBEV=8.617333262145e-5
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def cumulative(t,rate):
    return np.r_[0.,cumulative_trapezoid(rate,t)]

def main():
    histpath=P.parent/'geometry/temperature-histories.json'
    spectra_path=P.parent/'stopping/surface-spectra.json'
    stopping_path=P.parent/'stopping/surface-stopping.json'
    hist=json.loads(histpath.read_text()); spectra=json.loads(spectra_path.read_text()); stopping=json.loads(stopping_path.read_text())
    checks=[]; rows=[]; old=[]; sources={str(p):sha(p) for p in [histpath,spectra_path,stopping_path]}
    def check(name,ok):
        checks.append(dict(name=name,pass_=bool(ok)))
        assert ok,name
    for profile in hist['profiles']:
        path=Path(profile['csv']); promptpath=Path(profile['prompt_interval_upper']['short_time_upper_csv'])
        check('profile hash '+path.name,sha(path)==profile['sha256'])
        check('prompt hash '+promptpath.name,sha(promptpath)==profile['prompt_interval_upper']['short_time_upper_sha256'])
        sources[str(path)]=sha(path);sources[str(promptpath)]=sha(promptpath)
        a=np.genfromtxt(path,delimiter=',',names=True); b=np.genfromtxt(promptpath,delimiter=',',names=True)
        tc=a['relative_to_prompt_arrival_s']; tp=b['relative_to_prompt_arrival_s']
        # Replace the unresolved first quarter second, rather than adding two prompts.
        before=tc<0;after=tc>.25
        t=np.r_[tc[before],tp,tc[after]]
        check('increasing time '+path.name,np.all(np.diff(t)>0))
        stage=profile['stage']; sd=spectra['stages'][stage]
        arrival=np.array(sd['packet_arrival_after_ignition_s'])-(profile['prompt_arrival_s']-profile['focus_s'])
        weight=np.array(sd['species']['he4']['packet_nuclei_m2']);weight/=sum(weight)
        depth=np.array(stopping['stages'][stage]['species']['he4']['packet_mean_projected_depth_um'])*1e-6
        for key,col,pcol in [('bulk_k50','front_bulk_k50_K','bulk_front_upper_K'),('100um_skin_k4','front_with100um_skin_k4_K','skin4_front_upper_K')]:
            # +0.1K is the archived omitted-old-tail allowance, not measured scatter.
            T=np.r_[a[col][before],b[pcol],a[col][after]]+.1
            check('positive T '+path.name+key,np.all(T>0))
            D=5e-5*np.exp(-1.17/(KBEV*T));Dc=5e-5*np.exp(-1.17/(KBEV*np.minimum(T,1050.)))
            ad=cumulative(t,D);adc=cumulative(t,Dc)
            for deadline in [100.,1000.,4000.]:
                end=arrival+deadline;covered=(arrival>=t[0])&(end<=t[-1]);cw=float(sum(weight[covered]))
                start_clip=np.clip(arrival,t[0],t[-1]);end_clip=np.clip(end,t[0],t[-1])
                observed=np.maximum(end_clip-start_clip,0);missing=deadline-observed
                req=[]
                for nu in [1e12,1e13]:
                    def survival(E,missing_T,conditional=False):
                        c=cumulative(t,nu*np.exp(-E/(KBEV*T)))
                        hazard=np.interp(end_clip,t,c)-np.interp(start_clip,t,c)
                        hazard+=missing*nu*np.exp(-E/(KBEV*missing_T))
                        if conditional:return float(weight[covered]@np.exp(-hazard[covered])/cw)
                        return float(weight@np.exp(-hazard))
                    exact=brentq(lambda E:survival(E,1000.,True)-.95,1.,10.) if cw>0 else None
                    lower=brentq(lambda E:survival(E,1000.)-.95,1.,10.)
                    upper=brentq(lambda E:survival(E,1668.)-.95,1.,10.)
                    req.append(dict(attempt_frequency_s=nu,E95_covered_cohort_eV=exact,E95_all_cohort_eV_uncovered_at1000K=lower,E95_all_cohort_eV_uncovered_at1668K=upper))
                    check('ordered missing-temperature bound '+str((path.name,key,deadline,nu)),lower<=upper+1e-10)
                    if exact is not None: check('root survival '+str((path.name,key,deadline,nu)),abs(survival(exact,1000.,True)-.95)<1e-10)
                mob={}
                if cw>0:
                    for policy,integ in [('Jung_extrapolate',ad),('Jung_hold1050',adc)]:
                        A=np.interp(end[covered],t,integ)-np.interp(arrival[covered],t,integ)
                        mob[policy]=float(weight[covered]@erfc(depth[covered]/(2*np.sqrt(A)))/cw)
                rows.append(dict(stage=stage,selection=profile['selection'],temperature_model=key,repair_after_local_implantation_s=deadline,fully_covered_He_number_fraction=cw,requirements=req,conditional_mobile_front_release_fraction_covered=mob))
            # An existing trapped cohort sees the prompt and the later material flash.
            for deadline in [100.,1000.,4000.]:
                ent=[]
                for nu in [1e12,1e13]:
                    def hz(E):
                        c=cumulative(t,nu*np.exp(-E/(KBEV*T)))
                        return float(np.interp(deadline,t,c)-np.interp(0,t,c))
                    E=brentq(lambda E:hz(E)+math.log(.95),1,10)
                    # Isolate first0.25s contribution at the solved barrier.
                    c=cumulative(t,nu*np.exp(-E/(KBEV*T)))
                    hshort=float(np.interp(.25,t,c)-np.interp(0,t,c))
                    ent.append(dict(attempt_frequency_s=nu,E95_eV=E,first_quarter_second_hazard_fraction=hshort/(-math.log(.95))))
                old.append(dict(stage=stage,selection=profile['selection'],temperature_model=key,hold_after_prompt_s=deadline,requirements=ent))
    # Independent constant-temperature identity validates the rate/time convention.
    for T in [1000.,1668.]:
        E=KBEV*T*math.log(1e13*4000/(-math.log(.95)))
        check('constant T identity '+str(T),abs(math.exp(-1e13*4000*math.exp(-E/KBEV/T))-.95)<1e-14)
    out=dict(scope=__doc__,code_sha256=sha(__file__),inputs_sha256=sources,
      model={'trap_law':'S=exp[-integral nu exp(-E/kT) dt]; single first-order trap, no retrapping or irradiation detrapping',
       'acceptance':'95% survival of the He number-weighted cohort, assuming all selected atoms initially occupy the same trap type',
       'attempt_frequency_s':[1e12,1e13], 'covered_cohort':'Only arrivals whose complete postimplant hold interval fits the archived -600..4000s profile. No omitted-tail renormalization is used in the all-cohort bounds.',
       'missing_interval':'Bracket unrecorded durations at1000K and1668K. These are explicit conditional temperature bounds, not a reconstructed trajectory or credible interval.',
       'surface_temperature':'Front-surface T(t) used for micron implants, plus0.1K old-tail allowance. Prompt upper curve replaces the first0.25s; newly implanted atoms do not encounter their own earlier prompt.',
       'mobile_kernel':'erfc(depth/(2sqrt(integral Ddt))) in an absorbing-front half-space. Jung mobile/orientation D reference; mean projected depths omit straggling. Neither a trapped fraction nor real total leakage is inferred.',
       'energy':'E is an activation barrier, not the input work per recovered atom or defect.'},
      postimplant_rows=rows,preexisting_trap_rows=old,checks=checks)
    (P/'thermal-trapping.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'checks':len(checks),'all_pass':all(x['pass_'] for x in checks),'B1': [r for r in rows if r['stage']=='B1' and r['temperature_model']=='100um_skin_k4']},indent=2))
if __name__=='__main__':main()
