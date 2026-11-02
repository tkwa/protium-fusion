"""Independent arithmetic audit of the compact opaque station's selected row.

A passing budget is conditional on the stated future component performance.
This checks arithmetic, rating, radiative view area and finite launch tails.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
YEAR=365.25*86400
MEV=1.602176634e-13
G=6.67430e-11
SB=5.670374419e-8
checks=[]
def near(label,a,b,tol=2e-10):
 assert abs(a-b)<=tol*max(abs(a),abs(b),1),(label,a,b)
 checks.append(label)

def main():
 packet=json.loads((HERE/'opaque-station.json').read_text());s=packet['cases'][0];a=s['parameters'];p=s['power'];o=s['operations']
 variants=json.loads((HERE/'small-stages/cycle-He3-check.json').read_text())['variants']
 d=next(x for x in variants[0]['cases'] if x['collection_loss_per_stage']==1e-5)
 H=d['nuclear_Q_MeV_per_helium']-d['neutrino_MeV_per_helium'];U=d['fuel_preparation_MeV_per_helium']
 near('nuclear heat includes all nonneutrino output',H,s['H_MeV_He'])
 near('prepared input',U,s['U_MeV_He'])
 M=a['mass_kg'];R=a['collector_radius_m'];A=4*math.pi*R**2
 support=.001*G*M*M/(4*R)/(2*math.pi*R/1e4)
 maintenance=1e-4*M*1e8/YEAR
 near('support',support,p['support_W']);near('maintenance',maintenance,p['replacement_W'])
 specific_net=a['thermal_efficiency']*a['capture']*(H+U)-U/a['driver_efficiency']
 rate=(p['useful_W']+support+maintenance)/(.9*specific_net*MEV)
 near('independent throughput',rate,o['helium_atoms_per_s'])
 escaped=(1-a['capture'])*(H+U)*rate*MEV
 near('global nuclear energy closure',H*rate*MEV,p['waste_W']+p['useful_W']+escaped)
 near('electric conversion closure',p['gross_generation_W'],p['driver_W']+p['other_auxiliary_W']+support+maintenance+p['useful_W'])
 near('mass closure',sum(s['mass_budget_kg'].values()),M)
 assert s['mass_budget_kg']['unassigned_installed_design_reserve']>0;checks.append('positive reserve')
 flux=a['radiator_emissivity']*SB*a['radiator_temperature_K']**4
 near('radiator area',p['waste_W']/flux,o['radiator_area_m2'])
 sky=1361/packet['logistics']['destination_AU']**2/4
 radiator_only_area=p['waste_W']/(flux-sky)
 # Conservative whole-project sky budget: all outer surfaces absorb, while
 # only the selected panels radiate at330 K. Extra panels use fixed reserve.
 sky_absorbed=sky*A
 sky_area=(p['waste_W']+sky_absorbed)/flux
 assert sky_area<A;checks.append('radiator fits including full solar luminosity')
 receivers=json.loads((HERE/'opaque-receiver-audit.json').read_text())
 rating={z['stem']:z['source_rating'] for z in receivers['cases']}
 max_energy=max(z['initial_electric_energy_J'] for z in rating.values())
 near('max pulse electric energy',max_energy,o['maximum_electric_pulse_J'])
 near('two pulse energy storage',2*max_energy/1e5,s['mass_budget_kg']['pulse_storage_two_max_pulses_at100kJ_kg'])
 peak=max(z['peak_kinetic_launch_electric_W'] for z in rating.values())
 # Tight/no-rad B2 are alternatives, not two simultaneously installed stages.
 peaks=[v['peak_kinetic_launch_electric_W'] for k,v in rating.items() if not k.endswith('-norad')]
 b1=HERE/'small-stages/scaled-hot-B1-64-R30-initial.json';initial=json.loads(b1.read_text())
 r=np.asarray(initial['outer_radius_m']);v=-np.asarray(initial['velocity_m_s']);w=np.asarray(initial['mass_fraction_of_total']);ri=np.r_[0,r[:-1]];vi=np.r_[0,v[:-1]]
 rc=.75*(r**4-ri**4)/(r**3-ri**3);vc=vi+(v-vi)*(rc-ri)/(r-ri);flight=(R-rc)/vc
 idx=np.argsort(flight);quant=np.interp([.5,.9,.95,.99,.999],np.cumsum(w[idx])-.5*w[idx],flight[idx])
 receiver=next(z['receiver'] for z in receivers['cases'] if z['stem']=='scaled-hot-B1-64-R30')
 clear=receiver['arrival_1_50_99_percent_energy_s'][2];period=o['pulse_intervals_s']['n15'];gap=period-clear
 out=dict(scope=__doc__,baseline=dict(net_W=p['useful_W'],mass_kg=M,net_mW_kg=p['useful_W']/M*1e3,
   nuclear_heat_MeV_He=H,prepared_MeV_He=U,preaux_electric_MeV_He=specific_net,
   throughput_He_s=rate,mean_collector_flux_W_m2=p['captured_heat_W']/A,
   waste_W=p['waste_W'],replacement_W=maintenance,support_W=support,
   radiator_fraction_3K_only=o['radiator_area_fraction_of_outer_sphere'],
   radiator_fraction_with_full_site_luminosity=sky_area/A,
   destination_AU=packet['logistics']['destination_AU'],external_flux_mean_W_m2=sky,
   radiator_fraction_if_only_panels_absorb=radiator_only_area/A,
   additional_solar_absorbed_W=sky_absorbed,
   additional_panel_mass_from_reserve_kg=20*(sky_area-o['radiator_area_m2']),
   storage_max_pulse_J=max_energy,storage_two_pulse_kg=2*max_energy/1e5,
   peak_source_kinetic_electric_W=peak,sum_all_stage_peak_kinetic_electric_W=sum(peaks),
   peak_W_per_fixed_driver_kg=peak/s['mass_budget_kg']['driver_guidance_and_peak_power_hardware_allowance'],
   summed_peak_W_per_fixed_driver_kg=sum(peaks)/s['mass_budget_kg']['driver_guidance_and_peak_power_hardware_allowance'],
   unassigned_installed_reserve_kg=s['mass_budget_kg']['unassigned_installed_design_reserve'],
   reserve_after_external_load_panels_kg=s['mass_budget_kg']['unassigned_installed_design_reserve']-20*(sky_area-o['radiator_area_m2'])),
   B1_timing=dict(period_s=period,outgoing_99percent_energy_arrival_s=clear,
     available_incoming_window_after_previous_outgoing_99percent_energy_s=gap,
     flight_mass_quantiles_50_90_95_99_99point9_s=quant.tolist(),
     mass_fraction_flight_longer_than_window=float(w[flight>gap].sum()),
     mass_fraction_flight_longer_than_period=float(w[flight>period].sum()),
     note='Mass-bin screen, not a collision-free schedule. Energy-arrival percentile does not bound residual gas or mass clearance.'),
   checks=dict(count=len(checks),passed=checks),
   input_sha256={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'opaque-station.json',HERE/'opaque-receiver-audit.json',HERE/'small-stages/cycle-He3-check.json',b1,Path(__file__)]})
 (HERE/'opaque-station-audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
