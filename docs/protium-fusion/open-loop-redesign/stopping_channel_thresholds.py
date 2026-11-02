"""Mass-table nuclear reaction thresholds for illustrative graphite interception.

Q and kinematic thresholds do not determine cross sections or slowing survival.
"""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cycle_ledger import nuc

PROJECTILES=['p','he4','c12','c13','n14','n15']
TARGETS=['c12','c13']
# Atomic masses can be used directly because each listed strong reaction
# conserves Z. Gamma mass is zero. Electron binding uncertainties are tiny
# against the MeV-scale thresholds here, but no precision claim is intended.
ELEMENT={4:'be',5:'b',6:'c',7:'n',8:'o',9:'f',10:'ne',11:'na',12:'mg',13:'al'}
rows=[]
for projectile in PROJECTILES:
 for target in TARGETS:
  pn,tn=nuc(projectile),nuc(target)
  A,Z=pn.A+tn.A,pn.Z+tn.Z
  for emitted in ['gamma','n','p','he4']:
   en=None if emitted=='gamma' else nuc(emitted)
   final_A=A-(en.A if en else 0)
   final_Z=Z-(en.Z if en else 0)
   residue=ELEMENT[final_Z]+str(final_A)
   try: rn=nuc(residue);mr=rn.mass
   except Exception:continue
   final_mass=mr+(en.mass if en else 0)
   initial_mass=pn.mass+tn.mass
   Q=initial_mass-final_mass
   # Exact two-body lab kinetic threshold from invariant s, target at rest.
   threshold=max(0.,(final_mass**2-initial_mass**2)/(2*tn.mass))
   rows.append({'projectile':projectile,'target':target,'emitted':emitted,
      'residue':residue,'Q_MeV':Q,'lab_threshold_MeV':threshold,
      'elastic_channel':residue==target and emitted==projectile,
      'exothermic':Q>0,'interpretation':'Kinematically open at zero projectile energy; probability needs cross section and stopping integral.' if Q>0 else 'Threshold closes this specified channel below the listed energy; other channels may remain open.'})
# Inelastic excitation is nondestructive of carbon, included as energy-loss note.
result={'scope':'Mass-table channel inventory; illustrative strong reactions only, not exhaustive breakup channels or a stopping calculation.',
 'mass_source':'pynucastro installed nuclear mass table via cycle_ledger.nuc; same table as pulse energy audit',
 'formula':'E_lab,threshold=max(0,[(sum final masses)^2-(projectile mass+target mass)^2]/[2 target mass]), masses in MeV/c^2',
 'rows':rows,
 'notes':['No energy cutoff closes all carbon transmutation channels: radiative captures and several charged-particle emission channels are exothermic.',
          'Electromagnetic slowing can suppress charged-ion nuclear reactions through the Coulomb barrier, but a numerical loss claim requires cross sections integrated over the actual slowing path and isotope fractions.',
          'This calculation says nothing about sputtering, recombination, neutral bypass, field efficiency, or capture hardware mass.']}
Path(__file__).with_name('stopping-channel-thresholds.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
for r in rows:
 if r['projectile'] in ['p','he4'] or r['target']=='c12':
  print(f"{r['projectile']} + {r['target']} -> {r['residue']} + {r['emitted']}: Q={r['Q_MeV']:.6f} MeV, E_lab={r['lab_threshold_MeV']:.6f} MeV")
