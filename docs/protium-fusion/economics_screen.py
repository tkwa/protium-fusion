"""Transparent economic scales; assumptions, not future price forecasts."""
import json
from pathlib import Path

YEAR=365.25*86400
SPECIFIC_POWER=.001
P=1e9
HEAT_PER_KG_H=24.68*1.602176634e-13/(4*1.6735575e-27)
CONVERSION=.4
result={
  'threshold_net_specific_power_W_kg':SPECIFIC_POWER,
  'net_power_example_W':P,
  'installed_mass_at_threshold_kg':P/SPECIFIC_POWER,
  'annual_net_energy_per_kg_plant_J':SPECIFIC_POWER*YEAR,
  'hydrogen_nuclear_heat_J_kg':HEAT_PER_KG_H,
  'assumed_heat_to_net_efficiency':CONVERSION,
  'hydrogen_kg_per_GW_year_at_assumed_efficiency':P*YEAR/(HEAT_PER_KG_H*CONVERSION),
  'embodied_energy_payback_years':{
    str(e):e/(SPECIFIC_POWER*YEAR) for e in [1e5,1e6,1e7,1e8]},
  'illustrative_price_model':{
    'currency':'constant arbitrary monetary unit; not a future-dollar prediction',
    'annual_capital_recovery_factor':.05,'capacity_factor':.9,
    'capital_only_cost_per_kWh_per_unit_cost_per_kg':.05/(SPECIFIC_POWER*YEAR/3.6e6*.9),
    'unit_cost_per_kg_for_capital_only_0_05_per_kWh':.05*(SPECIFIC_POWER*YEAR/3.6e6*.9)/.05},
  'limitations':['threshold examples, not predicted performance',
    'no operations, replacement, financing history or infrastructure model',
    '40% net conversion here is an assumption; actual recirculation must be deducted',
    'all fuel-processing and required operating inventory must enter plant mass']}
if __name__=='__main__':
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
