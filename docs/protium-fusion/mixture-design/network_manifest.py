"""Archive the exact selected REACLIB network used by saved hydro runs."""
from pathlib import Path
import argparse, collections, hashlib, json, platform, sys
import numpy as np
import scipy
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(root,output):
    root=root.resolve();base=root/'docs/protium-fusion';sys.path.insert(0,str(base))
    import pulse_model as pm
    import pynucastro as p
    from pynucastro.rates.files import _find_rate_file
    libfile=Path(_find_rate_file('reaclib_default2_20250330'))
    sources=Path(_find_rate_file('rate_sources.csv'))
    temps=[.01,.03,.1,.2,.3,.4,.5,.6,.8,1.,1.5,2.]
    rates=[];max_error=0.
    for j,r in enumerate(pm.RATES):
        values=[]
        for t in temps:
            tf=np.array([1.,1/t,t**(-1/3),t**(1/3),t,t**(5/3),np.log(t)])
            own=float(sum(np.exp(np.array(s.a)@tf) for s in r.sets));ref=float(r.eval(t*1e9))
            if ref>1e-290:max_error=max(max_error,abs(own/ref-1))
            values.append(own)
        count=lambda seq:dict(collections.Counter(str(n).lower() for n in seq))
        rates.append({'index':j,'name':r.fname,'reaction':str(r),'reactants':count(r.reactants),'products':count(r.products),
                      'chapter':r.chapter,'source_label':r.src,'source':r.source,'label_properties':r.labelprops,
                      'weak_type':r.weak_type,'derived_from_inverse':bool(r.derived_from_inverse),
                      'density_exponent':int(r.dens_exp),'identical_particle_prefactor':float(r.prefactor),
                      'use_ye_weighting_in_rate_object':bool(r.use_ye_weighting),
                      'library_Q_MeV':float(r.Q),'hydro_Q_from_atomic_masses_MeV':float(pm.Q[j]),
                      'hydro_weak_Q_escape_MeV':float(pm.QNU[j]),
                      'stoichiometry':{n:float(pm.STOICH[i,j]) for i,n in enumerate(pm.SPECIES) if pm.STOICH[i,j]!=0},
                      'sets':[{'coefficients_a0_through_a6':list(s.a),'label_properties':s.labelprops,'label':s.label,
                               'resonant':bool(s.resonant),'weak':bool(s.weak),'derived_from_inverse':bool(s.derived_from_inverse)} for s in r.sets],
                      'sample_rates':values,'raw_reaclib_records':r.original_data})
    print("validation diagnostics",len(rates),pm.NN,max_error)
    assert len(rates)==194 and pm.NN==27 and max_error<2e-10
    assert np.max(abs(pm.AA@pm.STOICH))==0
    out={'scope':'Exact archival data for the selected 27-species/194-rate hydro network. A library date does not imply every constituent measurement is that recent.',
         'versions':{'python':platform.python_version(),'pynucastro':p.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
         'selection':'ReacLibLibrary().linking_nuclei(NUCS, with_reverse=True).get_rates(), in saved pulse_model.py order',
         'reaclib_library':{'filename':libfile.name,'sha256':sha(libfile)},'rate_source_database_sha256':sha(sources),
         'source_sha256':{'pulse_model.py':sha(base/'pulse_model.py'),'nuclear_screen.py':sha(base/'nuclear_screen.py')},
         'code_sha256':sha(__file__),'species':[{'name':n,'A':int(pm.AA[i]),'Z':int(pm.ZZ[i]),'atomic_mass_MeV':float(pm.MASS[i])} for i,n in enumerate(pm.SPECIES)],
         'rate_formula':'sum_sets exp[a0+a1/T9+a2*T9^(-1/3)+a3*T9^(1/3)+a4*T9+a5*T9^(5/3)+a6*ln(T9)]',
         'flow_formula':'rate * prefactor * rho_cgs^density_exponent * product(reactant molar abundances Y); no screening multiplier',
         'energy_convention':'Hydro Q is computed from species atomic masses, not the rounded library Q. All weak Q is subtracted from hot-stage deposition as an upper escape convention.',
         'sample_temperatures_GK':temps,'sample_rate_units':'s^-1 for one-body, cm^3 mol^-1 s^-1 for two-body, cm^6 mol^-2 s^-1 for three-body',
         'rates':rates,'checks':{'rate_count':len(rates),'species_count':pm.NN,'formula_vs_library_max_relative':max_error,'baryon_stoichiometry_exact':True}}
    output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out['checks']))
    for r in rates:
        if r['name'] in ['C12_p_to_N13_reaclib','C13_p_to_N14_reaclib','N14_p_to_O15_reaclib','N15_p_to_He4_C12_reaclib','N13_p_to_O14_reaclib']:
            print(r['name'],r['source_label'],r['source'],dict(zip(temps,r['sample_rates'])))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.root,a.output)
