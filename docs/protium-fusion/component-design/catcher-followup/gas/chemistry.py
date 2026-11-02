"""Chemical-form thresholds for a ten-percent total residual-gas column change.

No volatile production, cavity-entry fraction or permanent material loss is
assigned. Mass added from C/N is counted even when the hydrogen flow is fixed.
"""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent
d=json.loads((P/'inputs.json').read_text());a=json.loads((P/'chemistry-inputs.json').read_text());r=json.loads((P/'chamber-dt150.0-q128-cap0.1_eta100000.0.json').read_text());m=a['atomic_mass_kg'];k=1.380649e-23;T=1000.;R=d['R_m'];A=4*math.pi*R*R;g0=.01
qH,qHe=r['cavity_source_mass_rate_atomic_kg_s'];qHfull=d['plant']['recovered_H_kg_s'];fH=qH/qHfull
def v(mass):return math.sqrt(8*k*T/(math.pi*mass))
sigmaH=4*R*qH/(A*g0*v(2*m['p']));sigmaHe=4*R*qHe/(A*g0*v(m['he4']));base=sigmaH+sigmaHe;wH=sigmaH/base
escaped=[];returned=[]
for label,nH,other,number in [('C12H4',4,'c12',1),('C13H4',4,'c13',1),('N14H3',3,'n14',1),('N15H3',3,'n15',1)]:
 mass=nH*m['p']+number*m[other];added=number*m[other]/(nH*m['p'])
 for g in [.001,.01,.1]:
  ratio=mass/(nH*m['p'])*v(2*m['p'])/v(mass)*g0/g
  threshold=.1/(wH*(ratio-1)) if ratio>1 else None
  escaped.append({'molecule':label,'H_nuclei_per_molecule':nH,'global_capture_g':g,'mass_column_per_same_escaped_H_mass_relative_H2g0p01':ratio,'pressure_per_same_escaped_H_mass_relative_H2g0p01':2/nH*v(2*m['p'])/v(mass)*g0/g,'escaped_H_nucleus_fraction_for10percent_total_column_increase':threshold,'threshold_attainable_by_mixture':threshold is not None and threshold<=1,'extra_C_or_N_mass_per_escaped_H_mass':added,'added_C_or_N_kg_s_at_threshold':qH*threshold*added if threshold is not None else None})
  for f in [.0001,.001,.01,1.]:
   factor=f/fH*ratio
   y=.1/(wH*(factor-1)) if factor>1 else None
   returned.append({'molecule':label,'global_capture_g':g,'molecular_branch_H_cavity_entry_fraction':f,'fraction_of_total_returned_H_in_molecular_branch_for10percent_total_column_increase':y,'threshold_attainable_by_mixture':y is not None and y<=1,'comparison':'A fraction of total returned H takes this chemical branch, with its specified cavity-entry fraction; that fraction replaces the same fraction of baseline H2 routing. This separates chemical production from front transport.'})
# Statistical recombination of the exact returned N isotope flux, explicitly
# selected only to retain isotope masses. Pure-isotope alternatives differ little.
N=a['mean_fully_decayed_N_output_atomic_kg_s'];atoms={sp:N[sp]/m[sp] for sp in N};totalatoms=sum(atoms.values());x=atoms['n14']/totalatoms
components=[('N14N14',x*x,2*m['n14']),('N14N15',2*x*(1-x),m['n14']+m['n15']),('N15N15',(1-x)**2,2*m['n15'])]
nitrogen=[]
for g in [.001,.01,.1,1.]:
 all_sigma=sum(4*R*(totalatoms/2*f*mass)/(A*g*v(mass)) for _,f,mass in components)
 threshold=.1*base/all_sigma
 nitrogen.append({'global_capture_g':g,'all_returned_N_as_cavity_N2_column_kg_m2':all_sigma,'fraction_of_total_returned_N_entering_cavity_as_N2_for10percent_total_column_increase':threshold,'N_mass_to_cavity_at_threshold_kg_s':sum(N.values())*threshold,'scope':'The threshold combines actual chemical liberation and cavity routing. N retained in the collection layer or isolated process loop contributes no cavity gas. The material remains recoverable station inventory.'})
out={'criterion':'Ten-percent increase of the total steady H-bearing plusHe mass column at1000K, relative to ideal100nm-cap routing andH2/Heg=.01. This is a comparison criterion, not a universal pellet failure limit.','baseline':{'H_cavity_atomic_kg_s':qH,'He_cavity_atomic_kg_s':qHe,'H_fraction_entering_cavity':fH,'H2_column_kg_m2':sigmaH,'He_column_kg_m2':sigmaHe,'total_column_kg_m2':base,'H_column_fraction':wH,'total_returned_N_atomic_kg_s':sum(N.values())},'escaped_H_chemical_form_scenarios':escaped,'total_returned_H_chemical_and_routing_scenarios':returned,'N2_thresholds':nitrogen,'N2_isotope_model':[{'isotopologue':n,'molecular_number_fraction':f,'mass_kg':mass} for n,f,mass in components],'physical_interpretation':['At fixed hydrogen nuclei, CH4/NH3 add carbon/nitrogen mass and lower thermal speed. Their columns cannot be obtained by changing only molecular count.','A molecular compound released into closed channels can be dissociated/separated before return. Cavity-emitted molecules require physical ambient cleanup, whose ionization and final transfer differ fromHe.','CH4 carbon is a transfer from retained carbon inventory, not automatically a fresh-carbon debit. The same applies to nitrogen inNH3/N2.','Chemical-form fractions, release populations and species-dependent capture are unmeasured here; the tables expose the threshold for a subsequent chemistry test.','The formulas describe a steady diffuse reservoir. They neither bound directed molecular plumes nor solve high-energy reflected CNO atoms.','Atomic molecular-binding mass defects are neglected relative to atomic masses; this is negligible for the column threshold.'], 'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['inputs.json','chemistry-inputs.json','chamber-dt150.0-q128-cap0.1_eta100000.0.json']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'chemistry.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'baseline':out['baseline'],'same_g':[x for x in escaped if x['global_capture_g']==.01],'N2':nitrogen},indent=2))
