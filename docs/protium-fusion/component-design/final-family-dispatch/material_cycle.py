"""General material/energy ledger for actual isotope mixtures in staged pulses.

Input composition must be recorded or explicitly supplied. No isotope name or
stage label determines proton count. The default recycle quartet is C12/C13/
N14/N15; arbitrary input/output species and explicit cold decay maps are allowed.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.optimize import brentq
import pynucastro as pyna
MU=1.66053906660e-27
MEV=1.602176634e-13
U_MEV=931.49410242
DEFAULT_STATES=('c12','c13','n14','n15')
DEFAULT_DECAYS={'n':{'p':1.},'n13':{'c13':1.},'o14':{'n14':1.},'o15':{'n15':1.},
 'f17':{'o17':1.},'f18':{'o18':1.},'ne18':{'o18':1.},'ne19':{'f19':1.},
 'na21':{'ne21':1.},'na22':{'ne22':1.},'mg22':{'ne22':1.},'mg23':{'na23':1.}}
DEFAULT_COLD_NU={'n13':.706,'o14':1.05,'o15':1.00}
def nuc(name):return pyna.Nucleus(name)
def atomic_kg(name):return nuc(name).mass/U_MEV*MU
def add(dst,src,factor=1.):
 for k,v in src.items():dst[k]=dst.get(k,0.)+factor*v
def baryons(counts):return sum(nuc(k).A*v for k,v in counts.items())
def mass_MeV(counts):return sum(nuc(k).mass*v for k,v in counts.items())
def hashfile(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def clean_positive(values,tolerance=1e-7):
 neg=sum(-v for v in values.values() if v<0)
 if neg>tolerance:raise ValueError(f'Negative mass fraction {neg} exceeds {tolerance}')
 norm=sum(max(v,0.) for v in values.values())
 if norm<=0:raise ValueError('Empty composition')
 return {k:max(v,0.)/norm for k,v in values.items()},neg,norm

def cold_decay(counts,decay_map=None,neutrino_MeV=None):
 """Map each supplied parent to final daughters; maps may contain chains.

 Each dictionary value is daughter multiplicity/branch-averaged multiplicity.
 Nu energy may be specified per decay; absent values conservatively remove the
 entire positive Q. Unmapped species remain unchanged, so callers must supply
 relevant decays for newly added nuclides. Atomic masses include electron/positron
 bookkeeping consistently with the original ledger.
 """
 dm=DEFAULT_DECAYS if decay_map is None else decay_map
 nd=DEFAULT_COLD_NU if neutrino_MeV is None else neutrino_MeV
 def expand(name,trail=()):
  if name not in dm:return {name:1.},0.
  if name in trail:raise ValueError('Cyclic decay map')
  daughters=dm[name]
  if abs(baryons(daughters)-nuc(name).A)>1e-9:raise ValueError('Decay baryon failure '+name)
  q=nuc(name).mass-mass_MeV(daughters)
  if q < -1e-6:raise ValueError('Endothermic spontaneous decay '+name)
  nu=nd.get(name,max(0.,q));out={}
  if not 0<=nu<=max(0.,q)+1e-6:raise ValueError('Invalid decay neutrino energy '+name)
  for daughter,mult in daughters.items():
   fin,nuchain=expand(daughter,trail+(name,));add(out,fin,mult);nu+=mult*nuchain
  return out,nu
 out={};nu=0.
 for name,number in counts.items():
  fin,n=expand(name);add(out,fin,number);nu+=number*n
 q=mass_MeV(counts)-mass_MeV(out)
 if abs(baryons(counts)-baryons(out))>1e-9*max(1.,abs(baryons(counts))):raise ValueError('Cold decay baryons')
 return out,q,nu

def load_stage(path,*,initial_Y=None,initial_counts=None,trace_path=None,states=DEFAULT_STATES,
               decay_map=None,cold_neutrino_MeV=None,require_completed=True):
 """Load one pulse using actual input composition, normalize to one catalyst.

 Precedence: one explicit initial_Y or initial_counts; recorded initial_Y;
 recorded initial_mass_fractions; mass-weighted first saved trace state.
 Initial_Y is nuclei per baryon. initial_counts can use any common number basis.
 No legacy stage-name/proton-count inference is performed.
 """
 path=Path(path).resolve();d=json.loads(path.read_text());provenance={'record':str(path),'sha256':hashfile(path)}
 if initial_Y is not None and initial_counts is not None:raise ValueError('Choose one explicit composition')
 if initial_counts is not None:
  b=baryons(initial_counts)
  if b<=0:raise ValueError('Empty initial counts')
  y={k:float(v)/b for k,v in initial_counts.items()};origin='explicit_initial_counts'
 elif initial_Y is not None:y={k:float(v) for k,v in initial_Y.items()};origin='explicit_initial_Y'
 elif 'initial_Y' in d:y={k:float(v) for k,v in d['initial_Y'].items()};origin='recorded_initial_Y'
 elif 'initial_mass_fractions' in d:
  y={k:float(v)/nuc(k).A for k,v in d['initial_mass_fractions'].items()};origin='recorded_initial_mass_fractions'
 elif trace_path is not None:
  trace_path=Path(trace_path);t=np.load(trace_path);nn=len(t['species']);nz=int(t['zones']);nv=int(t['variables_per_zone'])
  z=t['state'][:,0].reshape((nz,nv));weights=t['cell_mass_kg']/sum(t['cell_mass_kg'])
  y=dict(zip(t['species'].tolist(),(weights@z[:,:nn]).tolist()));origin='mass_weighted_first_trace_state'
  provenance['initial_trace']={'path':str(trace_path),'sha256':hashfile(trace_path),
   'max_cell_initial_Y_deviation_from_mean':float(np.max(np.abs(z[:,:nn]-np.array(list(y.values())))))}
 else:raise ValueError('Actual initial composition required; stage label is insufficient')
 if (initial_Y is not None or initial_counts is not None) and 'initial_Y' in d:
  recorded=d['initial_Y']
  if max(abs(y.get(k,0.)-recorded.get(k,0.)) for k in set(y)|set(recorded))>1e-10:raise ValueError('Explicit initial composition disagrees with recorded actual composition')
 if min(y.values()) < -1e-14:raise ValueError('Negative initial composition')
 if abs(baryons(y)-1)>1e-8:raise ValueError('Initial_Y must satisfy sum A*Y=1')
 cat=sum(y.get(k,0.) for k in states)
 if cat<=0:raise ValueError('No designated circulating catalyst in input')
 inp={k:v/cat for k,v in y.items() if v>0};b=baryons(inp)
 if not math.isclose(sum(inp.get(k,0.) for k in states),1.,rel_tol=1e-12):raise ValueError('Catalyst normalization')
 x,neg,norm=clean_positive(d['output_mass_fractions'])
 hot={k:v*b/nuc(k).A for k,v in x.items() if v>0}
 out,qcold,nucold=cold_decay(hot,decay_map,cold_neutrino_MeV)
 qhot=mass_MeV(inp)-mass_MeV(hot)
 def per_attempt(field,fallback=None):
  if field in d:return d[field]/d['mass_kg']*b*MU/MEV
  if fallback in d:return d[fallback]*b*MU/MEV
  raise ValueError('Missing pulse energy '+field)
 u=per_attempt('initial_energy_J','initial_EOS_energy_J_kg')
 nuhot=per_attempt('weak_energy_escaped_upper_J','weak_energy_escaped_upper_J_kg')
 rad=per_attempt('radiation_escaped_J','radiated_J_kg')
 finalK=d.get('final_kinetic_J',0.);finalU=d.get('final_internal_J',0.)
 if finalK+finalU<=0:raise ValueError('Final kinetic/internal energies required')
 frac=finalK/(finalK+finalU)
 end=d.get('effective_termination_reason',d.get('termination_reason'))
 requested_duration=d.get('mixture_inputs',{}).get('duration_s')
 if requested_duration is not None:
  endpoint_ok=abs(d.get('elapsed_s',-1)-requested_duration)<=max(1e-12,1e-7*requested_duration)
 else:endpoint_ok=d.get('elapsed_sound_times',0)>=4-1e-8
 good=bool(d.get('success')) and end=='final_time' and endpoint_ok
 if require_completed and not good:raise ValueError('Pulse must genuinely complete its requested physical endpoint')
 return {'label':d.get('stage',path.stem),'source':str(path),'provenance':provenance,'initial_composition_source':origin,
  'initial_Y':y,'input_counts':inp,'cold_output':out,'hot_output':hot,'baryons_per_attempt':b,
  'protons_per_catalyst':inp.get('p',0.),'attempts_per_pulse':d['mass_kg']/(b*MU),'mass_kg':d['mass_kg'],
  'Q_hot_MeV':qhot,'Q_cold_MeV':qcold,'Q_total_MeV':qhot+qcold,'nu_hot_MeV':nuhot,'nu_cold_MeV':nucold,'nu_MeV':nuhot+nucold,
  'input_MeV':u,'radiated_MeV':rad,'final_kinetic_fraction':frac,'success':good,
  'hot_baryon_residual':baryons(inp)-baryons(hot),'clipped_negative_mass_fraction':neg,'output_mass_fraction_sum_before_normalization':norm,
  'pulse_energy_residual_fraction_initial':d.get('energy_residual_fraction_initial'),'requested_physical_duration_s':requested_duration,
  'cold_decay_map':DEFAULT_DECAYS if decay_map is None else decay_map,
  'cold_neutrino_MeV':DEFAULT_COLD_NU if cold_neutrino_MeV is None else cold_neutrino_MeV}

def solve_cycle(stages,*,carbon_return=1.,hydrogen_return=None,fresh_carbon_C13=.0107,
                eta_driver=.6,eta_heat=.6,capture=.99,states=DEFAULT_STATES,other_returns=None):
 """Stationary balance using actual catalyst inputs F, not identity.

 (F-rC*M)*attempts = fresh_carbon_vector. All stage input species are counted.
 Other species are supplied fresh unless an explicit other_returns fraction is
 given; their external input and ash are both reported. Net He=ashHe-freshHe.
 Return fractions mean EVENTUAL return after all retained queues/reprocessing.
 Incomplete first-pass extraction alone does not lower them. Energy capture is
 independent. Nuclear nonquartet products remain an explicitly held side stream
 unless a physical removal route is separately specified.
 """
 if hydrogen_return is None:hydrogen_return=carbon_return
 other_returns=other_returns or {}
 for x in [carbon_return,hydrogen_return,fresh_carbon_C13,eta_heat,capture,*other_returns.values()]:
  if not 0<=x<=1:raise ValueError('Fraction outside [0,1]')
 if not 0<eta_driver<=1:raise ValueError('Invalid driver efficiency')
 keys=list(stages);ss=[stages[k] for k in keys]
 if len(ss)!=len(states):raise ValueError('Square stage/state system required; overcomplete recipes need additional dispatch constraints')
 F=np.array([[s['input_counts'].get(k,0.) for s in ss] for k in states]);M=np.array([[s['cold_output'].get(k,0.) for s in ss] for k in states])
 f=np.array([(1-fresh_carbon_C13 if k=='c12' else fresh_carbon_C13 if k=='c13' else 0.) for k in states])
 A=F-carbon_return*M;attempts=np.linalg.solve(A,f)
 if np.min(attempts)<-1e-8:raise ValueError('No nonnegative stationary stage schedule for the specified recipe/feed')
 inputs={};outputs={}
 for s,a in zip(ss,attempts):add(inputs,s['input_counts'],a);add(outputs,s['cold_output'],a)
 returned={k:carbon_return*outputs.get(k,0.) for k in states}
 returned['p']=hydrogen_return*outputs.get('p',0.)
 for k,r in other_returns.items():
  if k in states or k=='p':raise ValueError('Specify CNO/H returns in their dedicated arguments')
  returned[k]=min(inputs.get(k,0.),r*outputs.get(k,0.))
 ext={k:inputs.get(k,0.)-returned.get(k,0.) for k in set(inputs)|set(returned)}
 ash={k:outputs.get(k,0.)-returned.get(k,0.) for k in set(outputs)|set(returned)}
 if min(ext.values()) < -1e-7 or min(ash.values()) < -1e-7:raise ValueError('Negative external/ash flow')
 helium=ash.get('he4',0.)-ext.get('he4',0.)
 if helium<=0 or ext.get('p',0.)<=0:raise ValueError('No positive net He output or fresh H consumption')
 def weighted(key):return float(sum(a*s[key] for a,s in zip(attempts,ss))/helium)
 q=weighted('Q_total_MeV');nu=weighted('nu_MeV');u=weighted('input_MeV');rad=weighted('radiated_MeV');h=q-nu
 K=I=0.
 for s,a in zip(ss,attempts):
  hot=s['input_MeV']+s['Q_hot_MeV']-s['nu_hot_MeV']-s['radiated_MeV']
  K+=a*hot*s['final_kinetic_fraction']/helium;I+=a*hot*(1-s['final_kinetic_fraction'])/helium
 delayed=weighted('Q_cold_MeV')-weighted('nu_cold_MeV')
 qcheck=(mass_MeV(ext)-mass_MeV(ash))/helium;bar=(baryons(ext)-baryons(ash))/helium
 closure=K+I+rad+delayed-h-u
 if abs(bar)>1e-8:raise ValueError(f'Baryon closure {bar}')
 if abs(qcheck-q)>1e-7:raise ValueError(f'Atomic-Q closure {qcheck-q}')
 if abs(closure)>1e-8:raise ValueError('Energy partition closure')
 extn={k:max(v,0.)/helium for k,v in ext.items()};ashn={k:max(v,0.)/helium for k,v in ash.items()}
 mC=sum(v*atomic_kg(k) for k,v in extn.items() if nuc(k).Z==6);mH=extn['p']*atomic_kg('p')
 preaux=eta_heat*capture*(h+u)-u/eta_driver
 pulse={k:{'stage':s['label'],'mass_kg':s['mass_kg'],'attempts_per_pulse':s['attempts_per_pulse'],
  'pulses_per_He':float(a/helium/s['attempts_per_pulse']),'input_counts':s['input_counts'],'protons_per_catalyst':s['protons_per_catalyst'],
  'gross_fuel_kg_per_He':float(a/helium*s['baryons_per_attempt']*MU)} for k,s,a in zip(keys,ss,attempts)}
 return {'feasible':True,'stage_attempts_per_He':dict(zip(keys,(attempts/helium).tolist())),
 'stage_pulse_data':pulse,'input_catalyst_matrix_F':F.tolist(),'output_catalyst_matrix_M':M.tolist(),'states':list(states),'stage_order':keys,
 'carbon_return':carbon_return,'hydrogen_return':hydrogen_return,'other_returns':other_returns,'fresh_carbon_C13_number_fraction':fresh_carbon_C13,
 'external_input_atoms_per_He':extn,'ash_atoms_per_He':ashn,'returned_atoms_per_He':{k:v/helium for k,v in returned.items()},
 'fresh_carbon_kg_per_He':mC,'fresh_hydrogen_kg_per_He':mH,'fresh_C_H_mass_ratio':mC/mH,
 'gross_H_input_atoms_per_He':inputs.get('p',0.)/helium,'gross_H_output_atoms_per_He':outputs.get('p',0.)/helium,
 'nuclear_H_consumed_atoms_per_He':(inputs.get('p',0.)-outputs.get('p',0.))/helium,
 'unreturned_H_atoms_per_He':(1-hydrogen_return)*outputs.get('p',0.)/helium,
 'gross_fuel_kg_per_He':sum(v['gross_fuel_kg_per_He'] for v in pulse.values()),
 'returned_fuel_kg_per_He':sum(atomic_kg(k)*v for k,v in returned.items())/helium,
 'nuclear_Q_MeV_per_He':q,'neutrino_MeV_per_He':nu,'H_nonnu_MeV_per_He':h,'U_prepared_MeV_per_He':u,
 'escaped_prompt_MeV_per_He':rad,'bulk_kinetic_MeV_per_He':K,'hot_internal_MeV_per_He':I,'delayed_decay_heat_MeV_per_He':delayed,
 'eta_driver':eta_driver,'eta_heat':eta_heat,'energy_capture_fraction':capture,'preaux_net_MeV_per_He':preaux,
 'preaux_electric_J_per_He':preaux*MEV,'baryon_residual_per_He':bar,'mass_defect_Q_residual_MeV_per_He':qcheck-q,
 'energy_partition_residual_MeV_per_He':closure,'stationary_matrix_condition_number':float(np.linalg.cond(A))}

def solve_supply_target(stages,*,target=.0056,total_budget=.007,hydrogen_return=None,**kwargs):
 """Largest coupled/independent CNO loss meeting the selected supply ratio.

 A nonmonotone or infeasible custom mixture must be scanned explicitly. The
 bracketed result reports carbon-only structure reserve; other materials have
 independent supply accounts. Fresh hydrogen includes its genuine losses.
 """
 if not 0<target<=total_budget:raise ValueError('Invalid carbon budget split')
 def run(loss):return solve_cycle(stages,carbon_return=1-loss,hydrogen_return=(1-loss if hydrogen_return is None else hydrogen_return),**kwargs)
 low=run(0.)
 if low['fresh_C_H_mass_ratio']>target:return {'feasible':False,'reason':'nuclear/recipe carbon demand exceeds budget at perfect CNO return','perfect_return_case':low}
 # Avoid singular/nonpositive schedules at the meaningless no-recycle endpoint.
 hi=1e-6;high=None
 while hi<=1.:
  try:high=run(hi)
  except (ValueError,np.linalg.LinAlgError):
   if hi>=.1:raise ValueError('Supply target not bracketed before schedule loses feasibility')
   hi*=2;continue
  if high['fresh_C_H_mass_ratio']>=target:break
  if hi==1:break
  hi=min(1.,hi*2)
 if high is None or high['fresh_C_H_mass_ratio']<target:raise ValueError('Supply target not bracketed')
 loss=brentq(lambda x:run(x)['fresh_C_H_mass_ratio']-target,0,hi,xtol=1e-14,rtol=1e-13)
 result=run(loss);result.update({'target_fuel_C_H_mass_ratio':target,'total_C_H_mass_budget':total_budget,
  'net_structural_C_H_allowance':total_budget-target,'net_structural_carbon_kg_per_He':(total_budget-target)*result['fresh_hydrogen_kg_per_He'],
  'minimum_CNO_return_fraction':1-loss,'interpretation':'Maximum allowed eventual unrecoverable-loss boundary; available carbon budget need not be consumed; not predicted or required apparatus loss.'})
 return result
