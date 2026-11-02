"""Homogeneous 2200 m/s equivalent thermal-flux catalyst inventory screen; not neutron transport.

Most tabulated cross sections are NIST total absorption values, so assigning
all absorption to the desired capture gives a lower bound on required stock.
V50 uses the 1991 measured radiative-capture value, 21(+4,-2) barn.
"""
from pathlib import Path
import json,sys
from math import log
from catalyst_catalog import atomicmass
from nonneutral_screen import MEV,C
ROOT=Path(__file__).parent
YEAR=365.25*86400
# NIST NCNR table, 2200 m/s absorption cross sections (barn).
ABS={'21045':27.5,'22046':.59,'22047':1.7,'22048':7.84,'22050':.179,
     '24050':15.8,'26054':2.25,'28058':4.6}
CAP={'23050':21.}

def screen(network,source_rate=1e24,flux=1e13):
    rows=[];uncovered=[]
    for key,rate in network['capture_steps_per_primary_pn'].items():
        if key not in ABS and key not in CAP:
            if rate:uncovered.append({'target':key,'capture_rate_per_primary_pn':rate})
            continue
        sig=CAP.get(key,ABS.get(key));z,a=divmod(int(key),1000)
        tau=1/(flux*sig*1e-28);mass=rate*source_rate*tau*atomicmass(z,a)*MEV/C**2
        rows.append({'target':key,'capture_per_primary_pn':rate,'cross_section_barn':sig,
          'cross_section_kind':'measured radiative capture' if key in CAP else 'total absorption used as capture upper bound',
          'mean_capture_time_years':tau/YEAR,'inventory_kg':mass})
    total=sum(r['inventory_kg'] for r in rows)
    return {'scope':__doc__,'source_rate_s':source_rate,'equivalent_2200m_s_flux_m2_s':flux,
      'inventory_rows':rows,'total_listed_inventory_lower_bound_kg':total,
      'uncovered_capture_steps':uncovered,'uncovered_captures_per_primary_pn':sum(r['capture_rate_per_primary_pn'] for r in uncovered),
      'Ti49_makeup_kg_year':network['Ti49_feed_for_boundary_per_primary_pn']*source_rate*atomicmass(22,49)*MEV/C**2*YEAR,
      'processing_sensitivities':[{'mean_inventory_processing_years':t,'specific_electricity_J_kg':e,
        'isotope_inventory_throughput_kg_s':total/(t*YEAR),'power_W':total*e/(t*YEAR),
        'Ti49_product_recapture_probability_before_extraction':(flux*2.2e-28)/(flux*2.2e-28+1/(t*YEAR))}
       for t in [100,1000,10000] for e in [1e7,1e8,1e9]]}

if __name__=='__main__':
    network=next(d for d in json.loads((ROOT/'catalyst-network-audit.json').read_text())['cases'] if d['T_MeV']==1.6)
    cases=[screen(network)]
    baseline=next(d for d in json.loads((ROOT/'end-transfer-sensitivity.json').read_text())['cases'] if d['assumptions']['nominal_B_T']==140 and d['assumptions']['heat_efficiency']==.7 and d['assumptions']['end_transfer_in']==.999)
    vol=baseline['mass_breakdown_kg']['water_moderator_and_shield']/1000
    nh=2000/((atomicmass(8,16)+2*atomicmass(1,1))*MEV/C**2)
    sigma_h=nh*.3326e-28
    light=network['light_net_before_D_completion'];useful_n=network['net_neutrons_available_for_light_fuels']-light.get('2003',0)
    flux=baseline['reaction_rate_s']*useful_n/(sigma_h*vol)
    d=screen(network,baseline['reaction_rate_s'],flux)
    d['homogeneous_balance_assumptions']={'water_volume_m3':vol,'macroscopic_H_absorption_per_m':sigma_h,
      'H_capture_fraction_of_source_neutrons':useful_n/(useful_n+network['all_neutrons_captured_on_catalysts']),
      'reference_plant_mass_kg':baseline['total_mass_kg'],
      'inventory_fraction_of_reference_plant':d['total_listed_inventory_lower_bound_kg']/baseline['total_mass_kg'],
      'limit':'Uses the earlier plant rate and water mass solely as a scale. Auxiliary ion channels and revised transport have not been allocated into that plant.'}
    cases.append(d)
    (ROOT/'capture-bank-screen.json').write_text(json.dumps({'scope':__doc__,'cases':cases},indent=2,allow_nan=False)+'\n')
    for d in cases:print(d['equivalent_2200m_s_flux_m2_s'],d['total_listed_inventory_lower_bound_kg'],d['Ti49_makeup_kg_year'])
