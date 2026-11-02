"""Conditional thermal-neutron competition in the saved 140 T plant.

The neutron spectrum, moderation and spatial flux are not solved. All included
materials see a common 2200 m/s equivalent scalar flux. Desired catalyst-bank
capture rates remain fixed by the prior network, so the residual light-fuel
neutrons compete between H and the added structural/impurity absorbers. No
useful energy is credited to parasitic captures. A coating is placed at the
unchanged plasma radius; the carbon backing is moved outward and thickened to
preserve its original wall-thickness/radius buckling ratio. It displaces water.
"""
from pathlib import Path
from math import pi,sqrt
import hashlib,json
from copy import deepcopy
import numpy as np
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
YEAR=365.25*86400
U=1.66053906660e-27
MEV=1.602176634e-13
SB=5.670374419e-8
# NIST NCNR compilation, absorption for 2200 m/s neutrons, natural elements.
# Atomic weights are explicit conventional inputs, not precision determinations.
# Densities are rounded engineering scenarios; alloy uses additive volumes.
ELEMENTS={
 'H':{'A':1.008,'abs_barn':.3326,'rho_kg_m3':None},
 'O':{'A':15.999,'abs_barn':.00019,'rho_kg_m3':None},
 'C':{'A':12.011,'abs_barn':.0035,'rho_kg_m3':2200},
 'B':{'A':10.81,'abs_barn':767.,'rho_kg_m3':2340},
 'N':{'A':14.007,'abs_barn':1.9,'rho_kg_m3':None},
 'Al':{'A':26.9815385,'abs_barn':.231,'rho_kg_m3':2700},
 'Ti':{'A':47.867,'abs_barn':6.09,'rho_kg_m3':4500},
 'V':{'A':50.9415,'abs_barn':5.08,'rho_kg_m3':6100},
 'Zr':{'A':91.224,'abs_barn':.185,'rho_kg_m3':6500},
 'Fe':{'A':55.845,'abs_barn':2.56,'rho_kg_m3':7870},
 'Ni':{'A':58.6934,'abs_barn':4.49,'rho_kg_m3':8900},
 'Cu':{'A':63.546,'abs_barn':3.78,'rho_kg_m3':8960},
}
COATINGS={e:{e:1.} for e in ['C','Al','Ti','V','Zr','Fe','Ni','Cu']}
COATINGS['Ti30Zr30V40']={'Ti':.3,'Zr':.3,'V':.4}

def read_reference():
 p=json.loads((ROOT/'network-plant.json').read_text())
 return next(r for r in p['cases'] if r['assumptions']['B']==140 and r['assumptions']['heat']==.7 and r['assumptions']['end_in']==.999)

def mixture(name):
 f=COATINGS[name]
 A=sum(v*ELEMENTS[e]['A'] for e,v in f.items())
 massf={e:v*ELEMENTS[e]['A']/A for e,v in f.items()}
 rho=1/sum(v/ELEMENTS[e]['rho_kg_m3'] for e,v in massf.items())
 return {'atomic_fractions':f,'mass_fractions':massf,'atomic_weight':A,
         'density_kg_m3':rho,'abs_barn':sum(v*ELEMENTS[e]['abs_barn'] for e,v in f.items())}

def geometry(r,thickness_m=0.,coating='C'):
 a=r['assumptions'];L=r['active_length_m'];coat=0.;carbon_extra=0.;displace=0.
 for c in r['channel_allocations']:
  rr,t,n=c['radius_m'],c['wall_m'],c['channels']
  dv=pi*(2*rr*thickness_m+thickness_m**2)*n*L
  factor=(1+t/rr)**2
  coat+=dv;carbon_extra+=dv*(factor-1);displace+=dv*factor
 water=r['water_volume_m3']-displace
 if water<=0:raise ValueError('coating exceeds available water volume')
 carbon=r['mass_breakdown_kg']['channel_walls']+carbon_extra*a['wall_density']
 mix=mixture(coating)
 return {'coating':coating,'thickness_m':thickness_m,'coating_volume_m3':coat,
         'coating_mass_kg':coat*mix['density_kg_m3'],'coating_material':mix,
         'carbon_wall_mass_kg':carbon,'water_volume_m3':water,
         'displaced_water_volume_m3':displace,'additional_carbon_mass_kg':carbon-r['mass_breakdown_kg']['channel_walls']}

def competition(r,thickness_m=0.,coating='C',water_boron_ppm=0.,wall_boron_ppm=0.,wall_nitrogen_ppm=0.,other_neutron_loss=0.,carbon_flux_ratio=1.,coating_flux_ratio=1.):
 g=geometry(r,thickness_m,coating);a=r['assumptions'];wm=g['water_volume_m3']*a['warm_water_density']
 wb=water_boron_ppm*1e-6;cb=wall_boron_ppm*1e-6;cn=wall_nitrogen_ppm*1e-6
 if max(wb,cb+cn)>=1:raise ValueError('invalid impurity mass fraction')
 # Match the original water number density, computed with H1/O16 atomic masses.
 h=r['H_macroscopic_capture_m_inv']*g['water_volume_m3']*(1-wb)
 o=h*ELEMENTS['O']['abs_barn']/(2*ELEMENTS['H']['abs_barn'])
 c=g['carbon_wall_mass_kg']*(1-cb-cn)/(ELEMENTS['C']['A']*U)*ELEMENTS['C']['abs_barn']*1e-28*carbon_flux_ratio
 bc=g['carbon_wall_mass_kg']*cb/(ELEMENTS['B']['A']*U)*ELEMENTS['B']['abs_barn']*1e-28*carbon_flux_ratio
 nc=g['carbon_wall_mass_kg']*cn/(ELEMENTS['N']['A']*U)*ELEMENTS['N']['abs_barn']*1e-28*carbon_flux_ratio
 bw=wm*wb/(ELEMENTS['B']['A']*U)*ELEMENTS['B']['abs_barn']*1e-28
 mix=g['coating_material'];film=g['coating_mass_kg']/(mix['atomic_weight']*U)*mix['abs_barn']*1e-28*coating_flux_ratio
 channels={'hydrogen':h,'oxygen':o,'carbon_wall':c,'boron_wall':bc,'nitrogen_wall':nc,'boron_water':bw,'coating':film}
 total=sum(channels.values());fractions={k:v/total for k,v in channels.items()}
 eta=(1-other_neutron_loss)*fractions['hydrogen']
 return {'scope':__doc__,'geometry':g,'water_boron_ppm':water_boron_ppm,'wall_boron_ppm':wall_boron_ppm,'wall_nitrogen_ppm':wall_nitrogen_ppm,
         'other_neutron_loss':other_neutron_loss,'carbon_flux_ratio':carbon_flux_ratio,'coating_flux_ratio':coating_flux_ratio,
         'absorption_inventory_m2':channels,'absorption_fractions':fractions,
         'hydrogen_capture_fraction_after_other_losses':eta,'capture_loss_fraction':1-eta}

def allocated_plant(r,c):
 """Reclose saved fixed-geometry plant with capture, material and maintenance changes.
 Parasitic-capture heat receives zero credit, consistent with the original
 conservative treatment. No fixed-source-profile optimization is done here.
 """
 a=r['assumptions'];g=c['geometry'];eta=c['hydrogen_capture_fraction_after_other_losses'];oldeta=a['neutron_capture']
 # Existing penalties are linear in missing light-fuel neutron captures.
 capture_penalty=r['capture_penalty_W']*(1-eta)/(1-oldeta)
 nuclear_base=json.loads((ROOT/'network-ion-budget.json').read_text())
 net=nuclear_base['network'];light_n=net['net_neutrons_available_for_light_fuels']-net['light_net_before_D_completion'].get('2003',0.)
 from nonneutral_screen import Q
 qnp=Q['np_capture']
 old_qnuc=r['nuclear_after_capture_W'];full_qnuc=nuclear_base['nuclear_heat_MeV_per_primary']*r['primary_pn_rate_s']*MEV
 qnuc=full_qnuc-(full_qnuc-old_qnuc)*(1-eta)/(1-oldeta)
 dgamma=(eta-oldeta)*light_n*qnp*r['primary_pn_rate_s']*MEV*a['gamma_leak']
 cold_radiation=r['cold_radiation_W']+dgamma
 cold=r['total_cold_load_W']+dgamma
 fridge=cold*(a['radiator_K']/a['coil_cold_K']-1)/a['cryo_Carnot_fraction']
 flux=r['primary_pn_rate_s']*light_n*eta/c['absorption_inventory_m2']['hydrogen']
 stock_ratio=r['capture_bank']['equivalent_2200m_s_flux_m2_s']/flux
 stock=r['mass_breakdown_kg']['capture_isotope_stock_lower_bound']*stock_ratio
 stock_power=r['capture_stock_processing_W']*stock_ratio
 # Recover the original preaux from bus, then apply changed accounted loads.
 preaux=r['plant_bus_net_W']/(1-a['other_fraction'])
 preaux+=r['capture_penalty_W']-capture_penalty
 preaux-=a['heat']*dgamma
 preaux+=r['refrigeration_W']-fridge
 preaux+=r['capture_stock_processing_W']-stock_power
 bus=preaux-a['other_fraction']*max(0.,preaux)
 masses={k:v for k,v in r['mass_breakdown_kg'].items() if k not in ['radiators','replacement_manufacturing_equipment','replacement_inventory']}
 masses['water_moderator_and_shield']=g['water_volume_m3']*a['warm_water_density']
 masses['channel_walls']=g['carbon_wall_mass_kg'];masses['conducting_or_getter_coating']=g['coating_mass_kg']
 masses['capture_isotope_stock_lower_bound']=stock
 masses['capture_bank_processing_equipment']=a['manufacturing_kg_W']*stock_power
 # The existing D-plant rating changes with successful H captures.
 masses['closed_D_plant']+=a['Dplant_kg_W']*r['primary_pn_rate_s']*light_n*(eta-oldeta)*Q['net_D']/2*MEV
 fixed=sum(masses.values())
 radflux=a['radiator_emissivity']*SB*(a['radiator_K']**4-a['ambient_K']**4)
 radkgW=a['radiator_areal_mass']/radflux;k=a['annual_replacement_fraction']*a['embodied_J_kg']/YEAR
 spare=a['annual_replacement_fraction']*a['replacement_inventory_years']
 M=(fixed+(qnuc-bus)*radkgW)/(1-k*(radkgW+a['manufacturing_kg_W'])-spare)
 maintenance=k*M;P=bus-maintenance;waste=qnuc-P
 masses.update(radiators=waste*radkgW,replacement_manufacturing_equipment=a['manufacturing_kg_W']*maintenance,replacement_inventory=spare*M)
 return {'net_W':P,'total_mass_kg':M,'net_W_kg':P/M,'plant_bus_net_W':bus,'capture_penalty_W':capture_penalty,
         'nuclear_after_capture_W':qnuc,'waste_heat_W':waste,'cold_radiation_W':cold_radiation,'refrigeration_W':fridge,
         'hydrogen_capture_fraction':eta,'thermal_flux_m2_s':flux,'capture_stock_ratio':stock_ratio,
         'capture_stock_kg':stock,'stock_processing_W':stock_power,'replacement_manufacturing_W':maintenance,
         'mass_breakdown_kg':masses,'mass_ledger_relative_error':sum(masses.values())/M-1,
         'energy_ledger_relative_error':(P+waste)/qnuc-1}

def main():
 r=read_reference();cases=[]
 hex_limit=min(c['radius_m']*(sqrt(pi/(2*sqrt(3))/r['assumptions']['packing'])/(1+c['wall_m']/c['radius_m'])-1) for c in r['channel_allocations'])
 for coat in COATINGS:
  for t in [0.,1e-7,1e-6,1e-5,1e-4]:
   c=competition(r,t,coat);c['plant']=allocated_plant(r,c);cases.append(c)
 impurities=[]
 for kind in ['water_boron_ppm','wall_boron_ppm','wall_nitrogen_ppm']:
  for value in [0.,.1,1.,10.,100.,1000.]:
   c=competition(r,**{kind:value});c['plant']=allocated_plant(r,c);impurities.append(c)
 tolerances=[]
 for coat in COATINGS:
  d={'coating':coat}
  for label,target in [('capture_99_5',.995),('specific_power_1mW',.001)]:
   def f(t):
    c=competition(r,t,coat)
    return (c['hydrogen_capture_fraction_after_other_losses'] if label.startswith('capture') else allocated_plant(r,c)['net_W_kg'])-target
   lo,hi=0.,hex_limit
   if f(lo)*f(hi)<0:d[label+'_thickness_m']=brentq(f,lo,hi,xtol=1e-15)
   else:d[label+'_thickness_m']=None
  tolerances.append(d)
 impurity_limits=[]
 for kind in ['water_boron_ppm','wall_boron_ppm','wall_nitrogen_ppm']:
  d={'impurity':kind}
  for label,target in [('capture_99_5',.995),('specific_power_1mW',.001)]:
   def f(ppm):
    c=competition(r,**{kind:ppm})
    return (c['hydrogen_capture_fraction_after_other_losses'] if label.startswith('capture') else allocated_plant(r,c)['net_W_kg'])-target
   d[label+'_ppm']=brentq(f,0.,1e5,xtol=1e-9)
  impurity_limits.append(d)
 flux_sensitivity=[]
 for coat,t in [('C',0.),('Ti30Zr30V40',1e-6),('Cu',1e-5)]:
  for ratio in [.5,1.,2.]:
   c=competition(r,t,coat,carbon_flux_ratio=ratio,coating_flux_ratio=ratio);c['plant']=allocated_plant(r,c);flux_sensitivity.append(c)
 area=sum(2*pi*c['radius_m']*c['channels']*r['active_length_m'] for c in r['channel_allocations'])
 inputs=['network-plant.json','network-ion-budget.json','network_plant.py','capture_bank_screen.py','nonneutral_screen.py']
 result={'scope':__doc__,'reference':{'net_W':r['net_W'],'total_mass_kg':r['total_mass_kg'],'net_W_kg':r['net_W_kg'],
  'segregated_equal_radius_hex_lattice_max_coating_m':hex_limit,'channel_inner_surface_area_m2':area,'wall_mass_kg':r['mass_breakdown_kg']['channel_walls'],'water_volume_m3':r['water_volume_m3']},
  'sources':[{'url':'https://www.ncnr.nist.gov/resources/n-lengths/list.html','use':'2200 m/s natural-element absorption cross sections; table attributes compilation to Neutron News 3(3), 1992, 29–37.'}],
  'element_inputs':ELEMENTS,'coating_definitions':{k:mixture(k) for k in COATINGS},'cases':cases,'impurity_cases':impurities,
  'coating_tolerances':tolerances,'impurity_limits':impurity_limits,'flux_ratio_sensitivity':flux_sensitivity,
  'input_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs}}
 (HERE/'capture-competition.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print('reference',result['reference'])
 for d in tolerances:print('coating limit',d)
 for d in impurity_limits:print('impurity limit',d)
 for coat in ['C','Ti','Zr','V','Ti30Zr30V40','Cu']:
  for c in cases:
   if c['geometry']['coating']==coat and c['geometry']['thickness_m'] in [1e-6,1e-5,1e-4]:
    print(coat,c['geometry']['thickness_m'],'loss',c['capture_loss_fraction'],'mW/kg',c['plant']['net_W_kg']*1000)
if __name__=='__main__':main()
