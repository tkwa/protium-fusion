"""Separate capital carbon stock from operational fresh-feed ratio."""
import json,hashlib
from pathlib import Path
P=Path(__file__).resolve().parents[1];cp=P/'cold-sink-redesign/cases.json';fp=P/'open-loop-redesign/operating-cases.json'
c=json.loads(cp.read_text())['cases']['working_200K_60pct_IR11_bypass1'];f=json.loads(fp.read_text())['cases']['high_source_rating_surface'];Y=365.25*86400
C=c['mass_breakdown_kg']['receiver'];h=f['feed_H_kg_s'];CH=0.007;rows=[]
for life in [20.,100.,1000.]:
 H=h*life*Y
 rows.append({'life_years':life,'fresh_H_delivered_kg':H,'receiver_capital_C_per_lifetime_fresh_H_kg_kg':C/H,'operational_C_ceiling_kg_at_0_007':H*CH,'scope':'Comparison of stocks and lifetime flows; operational ceiling is not imposed on construction.'})
ratio=4*1.00794/12.0107;byH=C*ratio;G=4.22689768e6
j={'scope':'Frozen thermal reference resource context. Receiver mass treated as carbon; alternative catcher/absorber materials may change this. Upstream material acquisition is not included in station W/kg.','input_sha256':{'thermal':hashlib.sha256(cp.read_bytes()).hexdigest(),'operating':hashlib.sha256(fp.read_bytes()).hexdigest()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'receiver_capital_C_kg':C,'receiver_fraction_of_installed_mass':C/c['mass_kg'],'fresh_H_kg_s':h,'operational_fresh_C_H_ceiling':CH,'lifetime_rows':rows,'methane_recipe':{'stoichiometry':'CH4 -> C(graphite) + 2 H2','H_byproduct_mass_per_C_mass':ratio,'stoichiometric_H_byproduct_kg':byH,'years_of_reference_H_feed_in_byproduct':byH/(h*Y),'CH4_feed_kg':C+byH,'standard298K_deltaG_J_kg_C':G,'work_at_50pct_efficiency_J':C*G/.5,'equivalent_years_reference_net_export_for_this_chemical_step':C*G/.5/(c['net_W']*Y),'scope':'Chemical step only, complete byproduct retained; not total factory, extraction, transport, installation or thermal-preparation cost. Standard-state JANAF value from manufacturing-efficiency/energy-examples.json.'}}
Path(__file__).with_name('capital-stock.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
