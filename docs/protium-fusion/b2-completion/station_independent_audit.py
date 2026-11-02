"""Independent nonlinear mass closure and global-energy audit of station cases.

Computes the mass root numerically from separately reconstructed components;
it does not use station.py's quadratic expression for its independent result.
Ratings and service allowances remain selected engineering requirements.
"""
import argparse,json,hashlib,sys,math
from pathlib import Path
from scipy.optimize import brentq
BASE=Path(__file__).resolve().parent.parent
MEV=1.602176634e-13;G=6.67430e-11;SB=5.670374419e-8;YEAR=365.25*86400

def independent(point,inputs):
 c=inputs['cycle'];r=inputs['receiver'];R=point['radius_m'];F=point['helium_per_s'];A=4*math.pi*R*R
 eta_d=.6;eta_h=.4;capture=.99;aux=.1;repfrac=1e-4;repcost=1e8;reserve=.2;vstream=1e4
 totalC=.007;Tr=point['radiator_K'];Tf=point['receiver_base_K'];sigma=point['absorber_sigma_kg_m2'];service=point['service_sigma_kg_m2'];rating=point['driver_peak_W_kg'];timing=point['timing_basis']
 H=c['H_nonnu_MeV_per_He']*MEV;U=c['U_prepared_MeV_per_He']*MEV
 captured_heat=F*capture*(H+U);gross=eta_h*captured_heat;driver=F*U/eta_d;auxiliary=aux*(gross-driver)
 source_peak=r['peak_source_electric_W_at_R1e9']*1e9/R
 storage_E=2*max(p['initial_energy_J'] for p in c['stage_pulse_data'].values())/eta_d
 fixed={'receiver':sigma*A,'insulation':100*A,'processing_structure_thermal_service':service*A,
   'peak_driver_hardware':source_peak/rating,'two_pulse_storage':storage_E/1e5,
   'average_driver_conditioning':driver*.1,'gross_electric_conversion':gross*.1,
   'working_inventory':c['gross_fuel_kg_per_He']*F*point.get('cold_process_buffer_s',86400.),
   'inflight_fuel_inventory':sum(p['gross_fuel_kg_per_He']*F*r['stages'][k]['timing_at_R1e9'][timing]['occupied_s']*R/1e9 for k,p in c['stage_pulse_data'].items()),
   'fresh_feed_buffer':(1+totalC)*c['fresh_hydrogen_kg_per_He']*F*30*86400}
 solar=A*1361/(25*4);emissive=.9*SB*Tr**4
 def physical(M):
  replacement=repfrac*M*repcost/YEAR
  # Independently start from stream kinetic energy, then charge0.1% per circuit.
  # Maximum radial Hill-frame tide is3GM_sunR/a^3. Add its
  # absolute envelope to the inward self-load, then divide by v^2/R.
  gself=G*M/(2*R*R)
  gtide=3*G*1.98847e30*R/(5*149597870700.)**3
  stream_mass=M*(gself+gtide)*R/(vstream*vstream)
  stream_energy=.5*stream_mass*vstream*vstream
  stream_loss=.001*stream_energy/(2*math.pi*R/vstream)
  P=gross-driver-auxiliary-replacement-stream_loss
  # Captured original nuclear energy less escaped preparation, then export.
  waste=captured_heat-F*U-P
  radarea=(waste+solar)/emissive
  components=dict(fixed, radiator_panels=20*radarea,replacement_manufacturing=replacement,
    ten_year_spares=repfrac*10*M,support_streams=stream_mass,installed_reserve=reserve*M)
  return components,P,waste,radarea,replacement,stream_loss
 scale=1e23
 def root(x):return sum(physical(x*scale)[0].values())/scale-x
 hi=1.
 while root(hi)>0 and hi<1e9:hi*=2
 if root(hi)>0:raise ValueError('No lower mass root bracket found')
 M=brentq(root,0,hi,xtol=1e-13)*scale
 comp,P,waste,radarea,replacement,support=physical(M)
 alpha=1-eta_h/.8
 def derivative(x):return 1/(Tf-x)**2-(1-eta_h)*alpha/(alpha*x-Tr)**2
 enginehot=brentq(derivative,Tr/alpha+1e-6,Tf-1e-6,xtol=1e-10);enginecold=alpha*enginehot
 UA=captured_heat*(1/(Tf-enginehot)+(1-eta_h)/(enginecold-Tr))
 occupancy=F*sum(c['stage_pulse_data'][k]['pulses_per_He']*s['timing_at_R1e9'][timing]['occupied_s'] for k,s in r['stages'].items())*R/1e9
 feedh=c['fresh_hydrogen_kg_per_He']*F;feedc=c['fresh_carbon_kg_per_He']*F;structC=(totalC-c['fresh_C_H_mass_ratio'])*feedh
 data={'mass_kg':M,'net_W':P,'net_W_kg':P/M,'gross_W':gross,'driver_W':driver,'auxiliary_W':auxiliary,
   'replacement_W':replacement,'support_W':support,'waste_W':waste,'radiator_area_fraction':radarea/A,
   'source_peak_W':source_peak,'stored_electric_J':storage_E,'absorbed_solar_W':solar,
   'engine_internal_hot_K':enginehot,'engine_internal_cold_K':enginecold,'heat_exchanger_conductance_W_K':UA,
   'sequential_occupancy':occupancy,'feed_H_kg_s':feedh,'feed_C_kg_s':feedc,'net_material_C_allowance_kg_s':structC,
   'net_material_C_allowance_kg_m2_year':structC/A*YEAR,
   'receiver_refurbishment_C_return_requirement':1-structC/(repfrac*sigma*A/YEAR),
   'mass_breakdown_kg':comp}
 relative={k:(v-point[k])/max(abs(v),abs(point[k]),1e-300) for k,v in data.items() if k!='mass_breakdown_kg'}
 components={k:(v-point['mass_breakdown_kg'][k])/max(abs(v),abs(point['mass_breakdown_kg'][k]),1e-300) for k,v in comp.items()}
 energy_residual=(F*H-P-waste-(1-capture)*(H+U)*F)/(F*H)
 mass_residual=sum(comp.values())/M-1
 # A simple through-thickness steady gradient, separate from pulse response.
 mean_flux=captured_heat/A;thickness=sigma/2000;gradient=mean_flux*thickness/100
 interpretations={'mean_captured_hot_heat_W_m2':mean_flux,'full_thickness_carbon_gradient_K_at_k100':gradient,
   'hot_exchanger_approach_K':Tf-enginehot,'cold_exchanger_approach_K':enginecold-Tr,
   'required_UA_per_service_kg_W_K_kg':UA/comp['processing_structure_thermal_service'],
   'gross_receiver_refurbishment_kg_s':repfrac*sigma*A/YEAR,
   'net_structural_C_atoms_per_He_approx':structC/F/(12*1.66053906660e-27),
   'source_peak_to_average_driver_ratio':source_peak/driver,
   'maximum_solar_tide_m_s2':3*G*1.98847e30*R/(5*149597870700.)**3,
   'self_gravity_load_m_s2':G*M/(2*R*R),
   'maximum_solar_tide_to_self_load_ratio':6*1.98847e30*(R/(5*149597870700.))**3/M}
 checks={'scalar_fields_match':max(abs(x) for x in relative.values())<1e-9,
   'mass_components_match':max(abs(x) for x in components.values())<1e-9,
   'nonlinear_mass_root_closes':abs(mass_residual)<1e-12,'global_energy_closes':abs(energy_residual)<1e-12,
   'finite_engine_approaches':Tf>enginehot>enginecold>Tr,'source_peak_exceeds_mean':source_peak>=driver}
 return {'independent':data,'relative_differences':relative,'component_relative_differences':components,
   'mass_residual':mass_residual,'global_energy_residual':energy_residual,'interpretation':interpretations,
   'checks':checks,'published_point':point}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--station-report');ap.add_argument('--output',required=True);a=ap.parse_args();inputs=json.loads(Path(a.inputs).read_text())
 sys.path.insert(0,str(BASE/'open-loop-redesign'));import station
 cycle=inputs['cycle'];receiver=inputs['receiver'];R=1.1e9;F=.8/(receiver['occupied_s_per_He_at_R1e9']['mass999']*R/1e9)
 points={}
 for name,rating,service in [('rounded_10Wkg_5000service',10.,5000.),('aggressive_100Wkg_1000service',100.,1000.)]:
  points[name]=station.evaluate(F,R,cycle,receiver,radiator_T=450,hot_T=1000,absorber_sigma=5000,driver_peak_W_kg=rating,service_sigma=service)
 if a.station_report:
  st=json.loads(Path(a.station_report).read_text());points['best_tested_family']=st['best_within_tested_family']
 rows={k:independent(p,inputs) for k,p in points.items() if p is not None}
 out={'scope':'Independent numeric nonlinear closure for stated engineering allocations and supplied feed. Does not establish hardware attainability, validated reactor power, or a global optimum.',
  'input_sha256':hashlib.sha256(Path(a.inputs).read_bytes()).hexdigest(),'station_code_sha256':hashlib.sha256((BASE/'open-loop-redesign/station.py').read_bytes()).hexdigest(),
  'cases':rows,'checks_passed':sum(sum(x['checks'].values()) for x in rows.values()),'checks_failed':sum(sum(not v for v in x['checks'].values()) for x in rows.values()),
  'limits':['1kg/W maintenance fabrication,10W/kg conversion and conditioning,100kJ/kg buffer and specified source rating are selected requirements.',
    'All remaining structural carbon assigned to receiver refurbishment gives a best-case minimum return; other carbon losses tighten it.',
    'Protection/recovery of slow tails, implantation/nuclear survival and engineering of the stated conductor/coolant paths remain separate.',
    'Full-thickness gradient is only an illustrative steady conduction bound at assumed k100,rho2000; periodic pulses are audited separately.',
    'Solar tidal support uses the maximum rotating Hill-frame acceleration as an absolute scalar envelope added to self-load; guide topology and orbit stability are not solved.']}
 Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'checks_passed':out['checks_passed'],'checks_failed':out['checks_failed'],'cases':{k:{'feasible':x['published_point']['feasible'],**{j:x['independent'][j] for j in ['mass_kg','net_W','net_W_kg','radiator_area_fraction']},'interpretation':x['interpretation']} for k,x in rows.items()}},indent=2))
if __name__=='__main__':main()
