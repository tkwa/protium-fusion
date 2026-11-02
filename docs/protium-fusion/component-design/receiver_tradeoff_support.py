"""Apply prestress to exact relaxed-capture masses, and finite bulk heat storage.
No hydro or repeated plant scan. Full-retained pulse energy is a conservative
hot-buffer scenario even when the plant routes some radiation into cold equipment.
"""
from pathlib import Path
import ast,json,math,hashlib
import numpy as np
P=Path(__file__).parent;BASE=P.parent;code=P/'receiver_prestress.py';trade=Path('/home/tkwa/src/protium-fusion/output/current-rate128-receiver-tradeoff.json')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
tree=ast.parse(code.read_text());nodes=[]
for x in tree.body:
 if isinstance(x,ast.For) and isinstance(x.target,ast.Name) and x.target.id=='sigma_receiver':break
 nodes.append(x)
ns={'__file__':str(code)};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(code)+'::definitions','exec'),ns)
d=json.loads(trade.read_text());R=ns['R'];A=ns['area'];G=ns['G'];n2=ns['n2'];omega=ns['omega'];E2D=ns['E2D'];nu=ns['nu'];MEV=1.602176634e-13
energies={}
for entry in ns['selected']:
 st=entry['stage_data'];r=entry['receiver'];stage=entry['record']['stage'];e=r['energy_channels_J'];late=(st['Q_cold_MeV']-st['nu_cold_MeV'])*st['attempts_per_pulse']*MEV
 energies[stage]={'kinetic_smooth_J_m2':r['ballistic']['kinetic_J_smooth_mass_conserving']/A,'kinetic_nodal_J_m2':e['outgoing_kinetic_nodal']/A,'internal_J_m2':e['endpoint_internal']/A,'prompt_J_m2':e['escaped_prompt_radiation']/A,'delayed_J_m2':late/A}
 energies[stage]['total_smooth_J_m2']=sum(energies[stage][k]for k in ['kinetic_smooth_J_m2','internal_J_m2','prompt_J_m2','delayed_J_m2'])
rows=[]
for r in d['rows']:
 if r['receiver_sigma_kg_m2']not in[500,750,1000] or r['additional_missed_hot_fraction_HplusU']!=.005 or r['missed_hot_routing']!='cold_equipment':continue
 p=r['plant'];M=p['mass_kg'];sigma=M/A;a=G*M/(2*R**3);T=max(1e4,(4.4*n2-a)*sigma*R*R/2);ap=2*T/(sigma*R*R);b=np.array([(a+ap+4*n2)/2,(a+ap-2*n2)/2,(a+ap-4*n2)/2]);ss=sigma*b/omega**2;sm=float(sum(ss)*A)
 mode=[ns['modes'](a,sigma,T,days)for days in[30.,100.,300.]];assert not any(x['growing']for x in mode)
 recv=r['receiver_sigma_kg_m2'];cp=1800.;rho=2000.;k=50.;L=recv/rho;alpha=k/(rho*cp);total=sum(e['total_smooth_J_m2']for e in energies.values());mean=sum(energies[s]['total_smooth_J_m2']/p['pulse_interval_s'][cat]for s,cat in [('A','n14'),('B1','n15'),('B2','c12'),('C13','c13')])
 row={'receiver_kg_m2':recv,'plant_mass_kg':M,'hot_fraction':r['total_hot_fraction'],'heat_routing':r['missed_hot_routing'],'He_s':p['helium_per_s'],'frame_tension_N_m':T,'biaxial_prestrain':(1-nu)*T/E2D,'cold_stress_for_8kg_m2_equivalent_load_layer_Pa':T/.004,'stream_axis_kg_m2':ss.tolist(),'stream_mass_kg':sm,'stream_mass_fraction_of_existing_scalar_budget':sm/p['mass_breakdown_kg']['support_streams'],'stream_loss_W_m2_at_0p001_per_turn':.001*sm*1e8/3/(2*math.pi/omega)/A,'affine_modes':mode,
 'thermal_buffer':{'density_kg_m3':rho,'cp_J_kg_K':cp,'bulk_k_W_m_K':k,'thickness_m':L,'bulk_diffusion_s':L*L/alpha,'full_retained_all_stage_bulk_swing_upper_K_constant_mean_pickup':total/(recv*cp),'largest_single_stage_uniform_rise_K':max(e['total_smooth_J_m2']for e in energies.values())/(recv*cp),'full_retained_mean_W_m2':mean,'mean_front_back_gradient_K_if_steady_1D':mean*L/k,'B1_interval_s':p['pulse_interval_s']['n15'],'large_event_one_day_guard_capacity_fraction':86400*sum(1/p['pulse_interval_s'][c]for c in['n14','c12','c13']),'bulk_buffer_scope':'For independently repeated stages and constant mean extraction, each zero-mean source contributes an energy range no greater than its total pulse energy, hence sum(E)/(Sigma cp). This is a bulk-mode energy bound, not a bound on local front/back temperature or a solved controller. Depositing all source heat here overstates hot-buffer energy when photons enter cold equipment. A finite temperature swing must be allowed; a fixed1000K back face would demand much larger instantaneous pickup.'}}
 assert row['thermal_buffer']['large_event_one_day_guard_capacity_fraction']<1
 rows.append(row)
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':{str(p):sha(p)for p in[code,trade,BASE/'mixture-design/current-rate128-reduced-inputs.json']},'pulse_fluences':energies,'rows':rows,'limits':['No source or receiver spectrum was relabeled: actual archived nominal128 quartet used.','50 W/m/K bulk conductivity is a requirement with separately qualified damaged50–100μm skin.','This bulk energy bound does not establish maximum surface temperature; prompt spectral deposition, contact, local plate stress and late-tail routing remain separate.','The relaxed-capture scalar routes specify the energy destination; actual radiation transport through cold services must realize that geometry.']}
(P/'receiver-tradeoff-support.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(rows,indent=2))
