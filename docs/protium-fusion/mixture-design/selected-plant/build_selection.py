"""Select the documented near-maximum receiver point and sampling sensitivity."""
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
sys.path.insert(0,str(HERE.parent))
import receiver_tradeoff as rt
import plant_compare as pc
import station_generalized as station

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    family_path=ROOT/'component-design/final-family/delivered-B1-C13H2/integration.json';trade_path=HERE/'receiver-tradeoff.json'
    family=read(family_path);trade=read(trade_path);cycle=family['precursor_cycle'];hot=family['hot_cycle']
    row=next(v for v in trade['rows'] if v['receiver_sigma_kg_m2']==350 and v['prompt_policy']=='EM20' and v['missed_hot_routing']=='cold_equipment' and v['additional_missed_hot_fraction_HplusUprecursor']==0)
    p=row['plant'];selected=family['plant_selected'];pc.cold=station
    fraction=(hot['H_nonnu_MeV_per_He']+hot['U_prepared_MeV_per_He'])/(cycle['H_nonnu_MeV_per_He']+cycle['U_prepared_MeV_per_He'])
    width=row['transport_sampling_80pct_half_width_HplusUhot']*fraction
    conditional=[]
    for delta in [-width,width]:
        x=rt.optimized(selected,cycle,350,1.,row['total_hot_fraction_HplusUprecursor']+delta)
        conditional.append(1000*x['net_W_kg'])
    stages=[]
    for s in family['hot_selected']:
        r=s['record'];cat=pc.CAT[r['stage']];energy=next(z['source_prepared_energy_J'] for z in selected if z['record']['stage']==r['stage'])
        stages.append({'stage':r['stage'],'catalyst':cat,'H_per_catalyst':s['stage_data']['input_counts']['p'],'mass_kg':r['mass_kg'],'precursor_prepared_J':energy,'hot_start_J':r['initial_energy_J'],'deposited_nuclear_over_precursor':(r['nuclear_energy_J']-r['weak_energy_escaped_upper_J'])/energy,'cold_outputs_per_initial_catalyst':s['stage_data']['cold_output'],'average_interval_s':p['pulse_interval_s'][cat]})
    checks={'positive_net':p['net_W']>0,'same_specific_power':abs(p['net_W']/p['mass_kg']/p['net_W_kg']-1)<1e-12,'full_mass_closure':abs(sum(p['mass_breakdown_kg'].values())/p['mass_kg']-1)<1e-12,'closed_material_return':cycle['carbon_return']==cycle['hydrogen_return']==1,'budget':cycle['fresh_C_H_mass_ratio']<.007,'energy_closure':max(abs(x) for x in p['closure_relative_residuals'].values())<1e-10,'sampling_brackets_central':min(conditional)<=1000*p['net_W_kg']<=max(conditional)}
    assert all(checks.values()),checks
    out={'family':family['family'],'selection':'350 kg/m² is the representative supported point on the near-tied250–350 kg/m² plateau. It has a lower calculated thermal envelope and smaller panel bow. The250 numerical score exceeds350 by0.043%, below Monte Carlo sampling uncertainty. No global optimum or complete hardware qualification is asserted.','plant':p,'cycle':cycle,'hot_cycle':hot,'stages':stages,'receiver_scenario':{k:v for k,v in row.items() if k!='plant'},'sampling_only_80pct_specific_power_mW_kg':sorted(conditional),'sampling_scope':'Propagated Monte Carlo transport sampling only, holding the source spectrum, geometry, all physical models and component assumptions fixed. This is not an80% credible interval for an engineered reactor.','checks':checks,'input_sha256':{str(q):sha(q) for q in [family_path,trade_path,Path(__file__),Path(rt.__file__),Path(pc.__file__),Path(station.__file__)]}}
    (HERE/'point.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'mW_kg':1000*p['net_W_kg'],'P_W':p['net_W'],'M_kg':p['mass_kg'],'sampling_only_80pct':conditional,'checks':checks},indent=2))
if __name__=='__main__':main()
