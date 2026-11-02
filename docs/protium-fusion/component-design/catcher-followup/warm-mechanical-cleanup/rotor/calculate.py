from pathlib import Path
import json,math,hashlib
p=Path(__file__).resolve().parent
raw=Path('/home/tkwa/src/protium-fusion/output/catcher-followup/warm-mechanical-cleanup/rotor')
area=math.pi*(.260/2)**2
out={'scope':'Primary benchmark normalization and conditional design scaling; no new pump operating point.','primary_1994':{'rotor_diameter_m':.260,'rotor_mass_kg':10,'rotor_rpm':27000,'suspension_input_W':13,'motor_input_W':35,'total_input_W':48,'rotor_disk_area_m2':area,'input_per_rotor_disk_W_m2':48/area},'conditional_target':{'idle_input_W_m2_pumpface':4,'pumpface_fraction':.245,'idle_W_m2_receiver':4*.245,'same48W_module_required_area_m2':12,'same48W_module_equivalent_radius_m':math.sqrt(12/math.pi),'rpm_at700m_s':700/math.sqrt(12/math.pi)*60/(2*math.pi)},'stress_scale':{'density_kg_m3':1800,'rim_speed_m_s':700,'thin_ring_Pa':1800*700**2,'ideal_disk_nu0p2_Pa':(3+.2)/8*1800*700**2},'thermal':{'black1000K_to300K_W_m2_viewed':5.670374419e-8*(1000**4-300**4),'scope':'Multiply by emissivity/view factor/cooled area ratio. Heat bypass is separate from electrical idle input.'},'nasa190lbf_N':190*4.4482216152605,'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in raw.glob('*.pdf')}}
assert out['primary_1994']['suspension_input_W']+out['primary_1994']['motor_input_W']==48
assert abs(out['conditional_target']['idle_W_m2_receiver']-.98)<1e-15
assert 844<out['nasa190lbf_N']<846
(p/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
