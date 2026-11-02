"""Physical radiance scope of grey-history Planck shape proxies.
At the saved emission area, a thermal diluted Planck source needs W<=1.
Quantiles below are energy-weighted source-bin distributions, not credible intervals.
"""
from pathlib import Path
import json,hashlib,math
import numpy as np
P=Path(__file__).parent;sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();source=P/'source-family-carbon-median128/integration.json';d=json.loads(source.read_text())['families']['source128'];rows={};hashes={str(source):sha(source)}
def quantile(x,w,ps):
 j=np.argsort(x);xx=x[j];ww=w[j];cw=np.cumsum(ww)/sum(ww);return {str(p):float(np.interp(p,cw,xx))for p in ps}
for cat,item in d['stage_inputs'].items():
 s=item['selected'];r=s['record'];tp=Path(s['receiver']['source_paths']['trace']);tr=np.load(tp);hashes[str(tp)]=sha(tp);nn=len(tr['species']);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);z=tr['state'].T.reshape(-1,nz,nv);m=tr['cell_mass_kg'];t=tr['time_s'];cum=(z[:,:,nn+4]@m)*float(tr['energy_unit_erg_g'])*1e-4;E=np.maximum(np.diff(cum),0);L=E/np.diff(t);Tout=z[:,-1,nn]*r['T0_K'];rr=z[:,-1,nn+1]*r['R0_m'];Tc=.5*(Tout[:-1]+Tout[1:]);radius=.5*(rr[:-1]+rr[1:]);F=L/(4*math.pi*radius**2);Teff=(F/5.670374419e-8)**.25;ok=E>0;E=E[ok];F=F[ok];Tc=Tc[ok];Teff=Teff[ok]
 cases={}
 for name,Tcol in [('outer_T',Tc),('outer_T_half',Tc*.5),('outer_T_quarter',Tc*.25),('surface_flux_equivalent_Teff',Teff)]:
  W=(Teff/Tcol)**4;cases[name]={'emission_energy_weighted_mean_dilution_W':float(E@W/sum(E)),'energy_weighted_W_quantiles':quantile(W,E,[.1,.5,.9,.99]),'fraction_prompt_energy_in_bins_W_gt1':float(sum(E[W>1+1e-12])/sum(E)),'energy_weighted_colour_over_Teff_quantiles':quantile(Tcol/Teff,E,[.1,.5,.9]),'energy_weighted_required_radius_over_saved_radius_quantiles':quantile(np.sqrt(W),E,[.1,.5,.9]),'scope':'W=(Teff/Tcolour)^4. If W>1, this colour cannot emit saved flux as a diluted equilibrium Planck spectrum from saved area. sqrt(W) is the emitting-radius multiplier needed at that colour if luminosity is unchanged.'}
 rows[r['stage']]={'total_prompt_J':float(sum(E)),'cases':cases}
out={'scope':__doc__,'code_sha256':sha(__file__),'input_sha256':hashes,'stages':rows,'interpretation':['Outer_T with W<=1 can be treated as a thermally diluted hard spectrum at the given area; this does not prove actual spectral thermalization.','Flux-equivalent Teff has W=1 by construction. It is the coolest single-temperature equilibrium Planck emitter that can carry the specified bolometric flux from that area.','A softer assigned shape with W>1 is an unconstrained non-LTE/spectral-shape sensitivity or needs a changed emitting area, reprocessing history or energy distribution. It is not automatically an attainable thermal photosphere or a credible worst case.','Non-LTE line/continuum emission and mixed-temperature areas can differ from these single-temperature proxies; actual frequency-dependent transfer remains unresolved.']}
(P/'source128-prompt-radiance.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:{n:v for n,v in z['cases'].items()if n in['outer_T','outer_T_quarter','surface_flux_equivalent_Teff']}for k,z in rows.items()},indent=2))
