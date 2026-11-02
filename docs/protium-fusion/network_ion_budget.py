"""Explicit ion-energy ledger for the TENDL recycling network.

Prompt CM neutron/photon moments set the thermal/direct partition. Remaining
prompt energy is provisionally collected as charged kinetic energy; evaluated
energy inconsistencies are reported. Inelastic same-isotope events remain in
the plasma and require replacing their photon energy. This is a source ledger,
not a finite driven reactor or a qualified energy converter.
"""
from pathlib import Path
import gzip,json,argparse
from math import pi,sqrt
import numpy as np
from catalyst_network import load_data,audit,stage,HOT,label
from catalyst_catalog import nucmass,atomicmass,thermal
from recycling_spectra import weighted_moment
from nonneutral_screen import MEV,C,Q,DD_CHARGED_PER_NEW_N,DD_THERMAL_PER_NEW_N,EPS0,MU0,E_CHARGE
from catalyst_penning import equilibrium,MP
ROOT=Path(__file__).parent


def stage_energy(za,d,spec,T):
    s=stage(za,d,T);z,a=za;sv=s['sv'];out=[]
    for e in s['events']:
        mt=e['mt'];rate=e['rate']; subs=spec['channels'][str(mt)];tabs=d['MF3'][mt]
        neut=0.;phot=0.;lightcharged=0.
        for sub in subs:
            v=weighted_moment(tabs,sub,T,z,a)/sv
            if sub['ZA']==0:phot+=v
            elif sub['ZA']==1:neut+=v
            elif sub['ZA']<=2004:lightcharged+=v
        mi=nucmass(1,1)+nucmass(z,a)
        nn=e['light'].get((0,1),0)
        nthermal=rate*nn*nucmass(0,1)/mi*1.5*T
        if mt==103:
            prep=phot;charged=0.;heat=phot;com=0.
        else:
            prep=rate*(e['Ecm_MeV']+1.5*T)
            heat=neut+nthermal+phot
            charged=prep+rate*e['Q_prompt_MeV']-heat
            com=1.5*T*rate
        # CM remainder after photons/light products is a diagnostic; retained
        # nuclear excitation and evaluated rounding may contribute to it.
        residual_cm=rate*(e['Ecm_MeV']+e['Q_prompt_MeV'])-neut-phot-lightcharged
        out.append({'mt':mt,'rate_per_pn':rate,'mean_reacting_Ecm_MeV':e['Ecm_MeV'],
          'input_thermal_MeV_per_pn':prep,'charged_prompt_MeV_per_pn':charged,
          'neutron_CM_MeV_per_pn':neut,'neutron_COM_thermal_MeV_per_pn':nthermal,
          'prompt_photon_MeV_per_pn':phot,'heat_prompt_MeV_per_pn':heat,
          'prompt_Q_MeV_per_pn':rate*e['Q_prompt_MeV'],
          'diagnostic_heavy_CM_remainder_MeV_per_pn':residual_cm,
          'neutrons_per_pn':rate*nn})
    return {'ZA':za,'sv':sv,'events':out,
      'extraction_events_per_pn':sum(e['rate_per_pn'] for e in out if e['mt']!=103),
      'input_thermal_MeV_per_pn':sum(e['input_thermal_MeV_per_pn'] for e in out),
      'charged_prompt_MeV_per_pn':sum(e['charged_prompt_MeV_per_pn'] for e in out),
      'neutron_MeV_per_pn':sum(e['neutron_CM_MeV_per_pn']+e['neutron_COM_thermal_MeV_per_pn'] for e in out),
      'photon_MeV_per_pn':sum(e['prompt_photon_MeV_per_pn'] for e in out),
      'inelastic_photon_MeV_per_pn':sum(e['prompt_photon_MeV_per_pn'] for e in out if e['mt']==103),
      'neutrons_per_pn':sum(e['neutrons_per_pn'] for e in out)}


def nuclear_energy(data,spectra,T=1.6,heater=.95,direct=.95,heat=.7):
    net=audit(data,T);stages=[]
    for za in HOT:
        s=stage_energy(za,data[za],spectra['targets'][label(za)],T)
        s['pn_per_primary']=net['processor_pn_rates'][label(za)];stages.append(s)
    def total(key):return sum(s['pn_per_primary']*s[key] for s in stages)
    n=net['net_neutrons_available_for_light_fuels']; light=net['light_net_before_D_completion']
    newD=n-light.get('2003',0)+light.get('1002',0)
    dcharged=newD*DD_CHARGED_PER_NEW_N
    # Trace T and He3 completion is credited entirely as heat, conservatively.
    charged=total('charged_prompt_MeV_per_pn')+dcharged
    prep=total('input_thermal_MeV_per_pn')
    q=net['heat_after_excluding_positive_makeup_transmutation_MeV']
    heatout=prep+q-charged
    g=direct+(1-direct)*heat
    electricity=g*charged+heat*heatout-prep/heater-net['discard_electricity_reserve_1000MeV_per_event_MeV']
    return {'T_MeV':T,'network':net,'stages':stages,
      'input_thermal_MeV_per_primary':prep,'charged_MeV_per_primary':charged,
      'thermal_MeV_per_primary':heatout,'nuclear_heat_MeV_per_primary':q,
      'energy_ledger_error_MeV':charged+heatout-prep-q,
      'bare_electricity_MeV_per_primary':electricity,
      'neutron_MeV_per_primary':total('neutron_MeV_per_pn'),
      'prompt_photon_MeV_per_primary':total('photon_MeV_per_pn'),
      'inelastic_photon_MeV_per_primary':total('inelastic_photon_MeV_per_pn'),
      'extraction_events_per_primary':total('extraction_events_per_pn'),
      'D_charged_MeV_per_primary':dcharged,
      'D_new_fuel_per_primary':newD,
      'assumptions':{'heater':heater,'direct':direct,'heat':heat,'trace_T_and_He3_completion_all_heat':True}}



def geometry(eq,xwall,stage,T=1.6,B=140,heater=.95,direct=.95,heat=.7,
             potential_in=.999,potential_out=.999,end_barrier_kT=25.,end_in=.999,end_out=.999,
             extraction_MeV=.2,atomic_MeV=.1):
    def v(k):return float(np.interp(xwall,eq['x'],eq[k]))
    z,a=eq['Z'],eq['A'];mc=MP*eq['M'];tj=T*MEV;sv=stage['sv']
    l=sqrt(tj*MP)/(E_CHARGE*B);ns=EPS0*B*B/MP;om=eq['xi']*z*E_CHARGE*B/mc;R=xwall*l
    rate=2*pi*l*l*ns*ns*sv*v('ihc')
    nh=2*pi*l*l*ns*v('ih');nc=2*pi*l*l*ns*v('ic')
    r2=l*l*v('ihcx2')/v('ihc');pairrot=.5*(MP+mc)*om**2*r2/MEV
    destr=stage['extraction_events_per_pn'];nn=stage['neutrons_per_pn'];mi=nucmass(1,1)+nucmass(z,a)
    neutronrot=pairrot*nn*nucmass(0,1)/mi;rotinput=destr*pairrot
    g=direct+(1-direct)*heat
    rotloss=rotinput/heater-g*(rotinput-neutronrot)-heat*neutronrot
    fh=2*pi*R*ns*v('h')*sqrt(tj/(2*pi*MP));fc=2*pi*R*ns*v('c')*sqrt(tj/(2*pi*mc))
    wallcirc=(fh*(2*T+.5*MP*om**2*R*R/MEV)+fc*(2*T+.5*mc*om**2*R*R/MEV))/rate
    wallloss=wallcirc*(1/heater-g)
    inside=eq['x']<xwall;xx=np.r_[eq['x'][inside],xwall]
    hh=np.interp(xx,eq['x'],eq['h']);cc=np.interp(xx,eq['x'],eq['c']);pp=np.interp(xx,eq['x'],eq['psi']);dd=np.interp(xx,eq['x'],eq['dpsi'])
    meanpsi=float(np.trapezoid(hh*cc*pp*xx,xx)/v('ihc'))
    potwork=destr*(z+1)*T*(meanpsi-v('psi'))
    potloss=potwork*(1/potential_in-potential_out-(1-potential_out)*heat)
    er=-tj/(E_CHARGE*l)*dd;rad=0.
    for den,charge,m in [(hh,1,MP),(cc,z,mc)]:
        force2=(er-om*xx*l*B)**2+2*B*B*tj/m
        pone=(charge*E_CHARGE)**4/(6*pi*EPS0*m*m*C**3)*force2
        rad+=2*pi*l*l*ns*np.trapezoid(den*pone*xx,xx)/rate/MEV
    mu=MP*mc/(MP+mc);qe=E_CHARGE*(mc-z*MP)/(MP+mc);K=z*E_CHARGE**2/(4*pi*EPS0)
    brems=(2/3)*qe*qe*K*tj/(EPS0*mu*mu*C**3)/sv/MEV
    radloss=(rad+brems)*(1/heater-heat)
    endcirc=destr*(z+1)*T*end_barrier_kT
    endloss=endcirc*(1/end_in-end_out-(1-end_out)*heat)
    extraction=destr*extraction_MeV;atomic=destr*atomic_MeV
    loss=rotloss+wallloss+potloss+radloss+endloss+extraction+atomic
    area=pi*R*R/rate
    return {'ZA':[z,a],'T_MeV':T,'B_T':B,'xi':eq['xi'],'h0':eq['h0'],'c0':eq['c0'],'xwall':xwall,
      'radius_m':R,'reaction_rate_per_m_s':rate,'area_per_pn_rate_m3_s':area,
      'proton_inventory_per_m':nh,'catalyst_inventory_per_m':nc,
      'pair_rotation_MeV':pairrot,'rotation_input_MeV_per_pn':rotinput,'neutron_rotation_MeV_per_pn':neutronrot,
      'rotation_loss_MeV_per_pn':rotloss,'wall_circulation_MeV_per_pn':wallcirc,'wall_loss_MeV_per_pn':wallloss,
      'potential_circulation_MeV_per_pn':potwork,'potential_loss_MeV_per_pn':potloss,
      'Larmor_MeV_per_pn':float(rad),'ion_brems_MeV_per_pn':brems,'radiation_loss_MeV_per_pn':float(radloss),
      'end_circulation_MeV_per_pn':endcirc,'end_loss_MeV_per_pn':endloss,
      'extraction_loss_MeV_per_pn':extraction,'atomic_preparation_reserve_MeV_per_pn':atomic,
      'total_loss_MeV_per_pn':float(loss),'axis_to_wall_V':-T*1e6*v('psi'),
      'wall_E_V_m':-tj/(E_CHARGE*l)*v('dpsi'),'rotation_speed_wall_over_c':om*R/C}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--spectra',type=Path,default=ROOT/'recycling-spectrum-moments.json.gz');ap.add_argument('--out',type=Path,default=ROOT/'network-ion-budget.json');args=ap.parse_args()
    data=load_data();spec=json.loads(gzip.decompress(args.spectra.read_bytes()));d=nuclear_energy(data,spec)
    args.out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
    print({k:v for k,v in d.items() if k not in ['network','stages']})
    for s in d['stages']:
        if s['pn_per_primary']>1e-6:print(s['ZA'],'flow',s['pn_per_primary'],'input',s['input_thermal_MeV_per_pn'],'charged',s['charged_prompt_MeV_per_pn'],'n',s['neutron_MeV_per_pn'],'gamma',s['photon_MeV_per_pn'],'inel',s['inelastic_photon_MeV_per_pn'])

