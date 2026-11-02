#!/usr/bin/env python3
"""Directional thermal environment screen at5AU; direct sunlight fixed tozero.
Run python environment.py --root /path/to/website [--output environment.json].
NumPy required only for angular quadrature. No assumed actual Dyson radius/temperature.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
from functools import lru_cache
SIGMA=5.670374419e-8;AU=149597870700.;LSUN=3.828e26;D=5*AU;TCMB=2.725
TEMPS=[197,200,250,300,450];EPS=.9

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
@lru_cache(None)
def gauss(n): return np.polynomial.legendre.leggauss(n)

def view(beta,theta,n=256):
    """Cosine-weighted view of a uniform bright disk with angular radius theta.
    beta is angle between panel normal and direction to center of disk.
    Analytic azimuth integral, Gaussian polar quadrature.
    """
    x,w=gauss(n);g=(x+1)*theta/2
    A=np.cos(beta)*np.cos(g);B=np.sin(beta)*np.sin(g);az=np.zeros_like(g)
    if abs(math.sin(beta))<1e-14:az=2*np.pi*np.maximum(A,0)
    else:
        full=A>=B;part=(A<B)&(A>-B);az[full]=2*np.pi*A[full]
        phi=np.arccos(np.clip(-A[part]/B[part],-1,1))
        az[part]=2*(A[part]*phi+B[part]*np.sin(phi))
    return float(np.sum(w*az*np.sin(g))*theta/(2*np.pi))
def qgray(Tr,eps,fself,env):
    """env=[(viewfactor, fixedexternalradiosityW/m²)]. Rest ofsky isCMB.
    Selfview is identical-temperature gray surface. Exact scalar radiosity closure.
    """
    fspace=1-fself-sum(f for f,J in env)
    assert fspace>=-1e-12
    H=sum(f*J for f,J in env)+max(0,fspace)*SIGMA*TCMB**4
    den=1-(1-eps)*fself
    J=(eps*SIGMA*Tr**4+(1-eps)*H)/den
    return {'net_rejection_W_m2':(1-fself)*J-H,'external_incident_W_m2':H,
        'first_pass_absorbed_external_W_m2':eps*H,'sky_view_factor':max(0,fspace),
        'self_view_factor':fself,'radiosity_W_m2':J}
def main(root,out):
    paths=[root/'docs/slop-requests/articles/solar-dyson-swarm.json',root/'docs/slop-requests/briefs/solar-dyson-swarm.json',root/'docs/solar-dyson-swarm/engineering-notes.json']
    flux=LSUN/(4*math.pi*D*D);cases=[];checks=[]
    for T in TEMPS:
        gross=EPS*SIGMA*T**4
        for ab in [0,.1,1,5,10,12.25,25,50,100,400]:
            net=gross-ab-EPS*SIGMA*TCMB**4
            cases.append({'radiator_K':T,'epsilon':EPS,'net_absorbed_environment_W_m2':ab,
                'net_rejection_W_m2':net,'emitting_area_per_W_m2_W':1/net if net>0 else None,
                'area_multiplier_vs_zero_IR':gross/net if net>0 else None})
    external=[]
    for f in [.5,.9,1.]:
        for rAU in [.4,1,2,4]:
            u=rAU/5;th=math.asin(u);T=(f*LSUN/(4*math.pi*(rAU*AU)**2*SIGMA))**.25
            J=SIGMA*T**4
            sphereF=(1-math.cos(th))/2
            directions=[]
            for deg in [0,30,60,90,120,150,180]:
                F=view(math.radians(deg),th)
                directions.append({'normal_to_source_angle_deg':deg,'view_factor':F,
                    'incident_W_m2':F*J,'absorbed_W_m2_alpha09':EPS*F*J,
                    'net_rejection_by_temperature_W_m2':{str(Tr):EPS*(SIGMA*Tr**4-F*J-(1-F)*SIGMA*TCMB**4) for Tr in TEMPS}})
            checks.extend([abs(directions[0]['view_factor']/u**2-1)<1e-12,
                abs(directions[3]['view_factor']-(th-math.sin(th)*math.cos(th))/math.pi)<1e-12,
                abs(directions[-1]['view_factor'])<1e-12])
            external.append({'IR_luminosity_fraction_Lsun':f,'source_sphere_radius_AU':rAU,
              'effective_black_outward_temperature_K':T,'angular_radius_deg':math.degrees(th),
              'normal_incident_W_m2':f*flux,'spherical_receiver_mean_incident_W_m2_exact':J*sphereF,
              'small_source_spherical_mean_incident_W_m2':f*flux/4,
              'spherical_receiver_mean_absorbed_W_m2_alpha09':EPS*J*sphereF,'directions':directions})
    geometries=[('clear_open_sky',0.,[]),('self_view_10pct',.1,[]),('self_view_50pct',.5,[]),
      ('300K_patch_view001',0.,[(.01,SIGMA*300**4)]),('300K_patch_view01',0.,[(.1,SIGMA*300**4)]),
      ('450K_shield_view001',0.,[(.01,SIGMA*450**4)]),('450K_shield_view01',0.,[(.1,SIGMA*450**4)]),
      ('1000K_collector_view0001',0.,[(.001,SIGMA*1000**4)]),
      ('1000K_collector_view001',0.,[(.01,SIGMA*1000**4)])]
    for Tb in [150,197,200,250,300,450]:geometries.append((f'complete_black_bath_{Tb}K',0.,[(1.,SIGMA*Tb**4)]))
    envrows=[]
    for label,fs,env in geometries:
        envrows.append({'geometry':label,'fixed_external_patches':[{'F':f,'J_W_m2':j} for f,j in env],
            'rows':[{'radiator_K':T,**qgray(T,EPS,fs,env)} for T in TEMPS]})
    enclosure=[]
    for rAU in [5,6,10]:
        T=(.9*LSUN/(4*math.pi*(rAU*AU)**2*SIGMA))**.25
        enclosure.append({'enclosing_shell_radius_AU':rAU,'luminosity_fraction':.9,
            'equal_inner_outer_black_temperature_K':T,'absorbed_W_m2_alpha09':EPS*SIGMA*T**4,
            'net_rejection_by_temperature_W_m2':{str(Tr):EPS*SIGMA*(Tr**4-T**4) for Tr in TEMPS}})
    warm=[]
    for Tr in TEMPS:
        for Tw in [200,250,300,450,1000]:
            # Equalgrayabsorption/emission, otherdirectionscold, no selfview.
            warm.append({'radiator_K':Tr,'warm_black_patch_K':Tw,
              'F_max_for_positive_net':min(1.,(Tr/Tw)**4),
              'F_max_for_80pct_clear_sky_capacity':min(1.,.2*(Tr/Tw)**4)})
    checks.extend([abs(qgray(200,EPS,0,[(1,SIGMA*200**4)])['net_rejection_W_m2'])<1e-10,
                   qgray(197,EPS,0,[(1,SIGMA*300**4)])['net_rejection_W_m2']<0])
    outdata={'scope':__doc__,'constants':{'sigma_W_m2_K4':SIGMA,'AU_m':AU,'Lsun_W':LSUN,'distance_AU':5,'CMB_K':TCMB},
       'accepted_direct_solar_absorbed_W_m2':0.,'epsilon_reference':EPS,
       'canonical_inputs_sha256':{str(p.relative_to(root)):sha(p) for p in paths},'code_sha256':sha(__file__),
       'canon_scope':'Layered independentlyorbitingcollectors;90%capture;directionalwasteheat/coolingrights. Exactorbitalradii/temperatures unspecified.0.4AU is archivedillustration,notadoptedradius.',
       'compact_external_benchmark':{'normal_incident_W_m2_per_fIR':flux,'spherical_mean_incident_W_m2_per_fIR':flux/4,
          'at_fIR09_spherical_incident_W_m2':.9*flux/4,'at_fIR09_alpha09_spherical_absorbed_W_m2':.9*.9*flux/4},
       'gross_emission_W_m2':{str(T):EPS*SIGMA*T**4 for T in TEMPS},'CMB_absorbed_W_m2_alpha09':EPS*SIGMA*TCMB**4,
       'absorbed_heat_cases':cases,'external_uniform_source_spheres':external,'local_view_geometries':envrows,
       'energy_balanced_enclosing_black_shells':enclosure,'warm_view_limits':warm,
       'area_accounting':{'emitting_area':'Sumofphysicalradiatingfaces. A_two_sided_panel=2A_planform onlyifbothfacescountedwiththeirviewfactors.',
         'plant_balance':'Q_cold_internal+Q_IR_abs_oncollector+Q_conduction+Q_other=Σ A_emit_i q_net_i; do notcountthe same externalIR twice.',
         'collector_routing':'Q_IR_collector=q_abs_collector*A_collector. Its contribution peremittingarea isQ_IR_collector/A_emit, notq_abs_collector repeatedonA_emit.',
         'projection':'Small centralisotropic source heats a fullysphericalreceiver byprojectedareaπR²; outwardemittingarea4πR² givesfactor1/4. A facing flatpanel usesitsprojectedarea,notfactor1/4.'},
       'assumptions':['Fixedblackexternalradiosity or statedisotropicluminositybenchmarks;actualswarmhasdirectionaland spectral structure.',
        'GrayalphaIR=epsilon=.9 fordirectionalcases. A spectralmodelmay use differentPlanckmeanabsorptionforwarmandcold spectra whileobeyingKirchhoffat each wavelength/direction.',
        'No uniform300K enclosingshell at5AU isimpliedbycanon. Itsenergybudget requiresadditionalpower or differentinner/outersurfaces/escapeemissivity.',
        'Opposite-facing lowIRpanels requireactualspaceview andthermalrouting. Shieldtemperature, occlusionandconduction remainengineeringinputs.',
        'Same-temperature selfview solvedas grayradiosity;external fixedsurfaces do notreactthermally in thisscreen. A finite shieldmustgets itsownbalance.',
        'No probability distribution overunknowngeometry isassigned. Tablesareconditionalcases.'],
       'sources':{'NASA_thermal_balance':'https://www.nasa.gov/smallsat-institute/sst-soa/thermal-control/',
         'Dyson_radiative_thermodynamics':'https://arxiv.org/abs/2309.06564'},'checks':{'passed':all(checks),'count':len(checks)}}
    out.write_text(json.dumps(outdata,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:outdata[k] for k in ['compact_external_benchmark','gross_emission_W_m2','energy_balanced_enclosing_black_shells','checks']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[3]);p.add_argument('--output',type=Path);a=p.parse_args()
 main(a.root,a.output or Path(__file__).with_suffix('.json'))
