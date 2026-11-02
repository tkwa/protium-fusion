#!/usr/bin/env python3
"""Uniform warm-cavity gas inventory and cold-packet acceptance; no gas loss inferred.
Source strengths are explicit root-ledger scenarios, not a new reaction-network result.
"""
import argparse,hashlib,json,math
from pathlib import Path
KB=1.380649e-23;MU=1.66053906660e-27;R=1.1e9;AREA=4*math.pi*R**2
RATE=3.78224727e32;MP=1.67262192369e-27;MHE=6.64647699e-27

def gas(Q,m,T):
 v=math.sqrt(8*KB*T/(math.pi*m));tau=4*R/(3*v);sig=4*Q*R/(AREA*v);n=sig/(R*m)
 return {'mass_source_kg_s':Q,'mass_kg':m,'temperature_K':T,'mean_thermal_speed_m_s':v,'unit_capture_residence_s':tau,'unit_capture_inventory_kg':Q*tau,'unit_capture_radial_column_kg_m2':sig,'unit_capture_pressure_Pa':n*KB*T,'unit_capture_number_density_m3':n}
def main():
 scenarios={};checks=[]
 for name,Hp,rate in [('closed_reference',1.85904844964,RATE),('A_H2_B1_H2_candidate',3.817,RATE),('provisional_A3_B1H2_B2H2p5_C13H4',3.997e6/(4.473e32*MP),4.473e32)]:
  rows={}
  for T in [300.,1000.,2000.]:
   H=gas(rate*Hp*MP,MP,T);H2=gas(rate*Hp*MP,2*MP,T);He=gas(rate*MHE,MHE,T)
   combinations={}
   for mode,gases in [('atomic_H_plus_He',[H,He]),('H2_plus_He',[H2,He])]:
    sig=sum(g['unit_capture_radial_column_kg_m2'] for g in gases);p=sum(g['unit_capture_pressure_Pa'] for g in gases);q=sum(g['mass_source_kg_s'] for g in gases)
    tests=[]
    for a,budget in [(a,b) for a in [.001,.01] for b in [1e3,1e4,1e5]]:
     limit=8*1000*a*budget/(3*(2.4e6)**2)
     tests.append({'grain_radius_m':a,'grain_density_kg_m3':1000,'grain_velocity_m_s':2.4e6,'allowed_gas_heating_J_kg':budget,'full_kinetic_retention':True,'maximum_radial_column_kg_m2':limit,'full_injection_unit_capture_heating_J_kg':budget*sig/limit,'maximum_common_cavity_injection_fraction_over_effective_capture_fraction':limit/sig,'maximum_common_actual_residence_s_if_all_gas_injected':limit*AREA/(3*q),'maximum_local_escape_fraction_by_global_capture':{str(g):min(1.,g*limit/sig) for g in [1.,.1,.01,.001]}})
    combinations[mode]={'radial_column_kg_m2':sig,'pressure_Pa':p,'total_gas_source_kg_s':q,'acceptance':tests,'thermal_gas_collision_mfp_m_for_cross_sections_m2':{str(s):1/(sum(g['unit_capture_number_density_m3'] for g in gases)*s) for s in [1e-20,1e-19,1e-18]}}
   rows[str(T)]={'H':H,'H2':H2,'He':He,'combinations':combinations}
   for g in [H,H2,He]:checks.append(abs(3*g['unit_capture_inventory_kg']/AREA/g['unit_capture_radial_column_kg_m2']-1)<1e-12)
  scenarios[name]={'He_rate_s':rate,'unburned_H_atoms_per_He':Hp,'rows':rows}
 return {'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'geometry':{'R_m':R,'area_m2':AREA,'He_rate_s':RATE},'scenarios':scenarios,'equations':{'Sigma':'4 Qmass R/(A vbar) * f/g','residence':'4 R/(3 vbar g)','pressure':'Sigma kT/(R m)','packet_heating':'3 Sigma vpacket^2/(8 rho_grain radius), full kinetic retention','eventual_port_escape':'fport/[fport+s(1-fport)], before downstream recovery'},'checks':{'count':len(checks),'passed':sum(checks),'all_passed':all(checks)},'assumptions':['Uniform steady injection over a sphere; gas thermalizes/reemits at selected wall temperature. g is area- and speed-weighted effective removal into the return system; f is fraction of total gas source injected into the central cavity.','A uniform Maxwellian is also a stationary collisional equilibrium: gas-gas collisions do not automatically alter these steady boundary balances. Heterogeneous sources, directional plumes, temperature gradients and selective traps require transport.','The mean collision cross-section values are diagnostic parameters, not measured species values or a credible interval.','He product also contributes cold-grain drag and must enter a separate return/product line. It is not a catalyst or fresh-H loss.','Local closed backplenums can contain long-lived gas inventory without exposing incoming packets; no global cryocondensation requirement is inferred.','Gas heating allocations 1, 10 and 100 kJ/kg and complete retention are chosen conservative tests. Material stopping/ablation may set stricter surface criteria.']}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default=str(Path(__file__).with_name('return-gas.json')));a=p.parse_args();o=main();Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:v['rows']['1000.0']['combinations'] for k,v in o['scenarios'].items()},indent=2))
