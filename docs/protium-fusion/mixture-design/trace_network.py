"""Saved-trace flow and screening diagnostics; no hydrodynamic rerun."""
from pathlib import Path
import argparse,hashlib,json,math,sys
import numpy as np
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
EPS=8.8541878128e-12;QE=1.602176634e-19;KB=1.380649e-23

def weighted_q(x,w):
    ix=np.argsort(x.ravel());z=x.ravel()[ix];s=np.cumsum(w.ravel()[ix])
    return [float(z[min(np.searchsorted(s,p*s[-1]),len(z)-1)]) for p in [.05,.5,.95]] if s[-1]>0 else [None]*3

def reduce(record,pm):
    d=json.loads(record.read_text());z=np.load(d['trace_path']);m=z['cell_mass_kg'];t=z['time_s'];nz=len(m);state=z['state'].reshape(nz,int(z['variables_per_zone']),len(t));names=d['state_variables']
    T=np.maximum(state[:,names.index('T/T0'),:]*d['T0_K'],1e6)
    r=np.vstack([np.zeros(len(t)),state[:,names.index('r/R0'),:]*d['R0_m']]);rho=m[:,None]/(4*np.pi/3*np.diff(r**3,axis=0));Y=np.maximum(state[:,:pm.NN,:],0.)
    dt=np.r_[.5*np.diff(t)[0],.5*(np.diff(t)[:-1]+np.diff(t)[1:]),.5*np.diff(t)[-1]]
    baseweights=m[:,None]/pm.MU*dt[None,:]
    ncat=d['mass_kg']/pm.MU*sum(y for n,y in d['initial_Y'].items() if n!='p')
    te=T/1e9;tf=np.stack([np.ones_like(te),1/te,te**(-1/3),te**(1/3),te,te**(5/3),np.log(te)],axis=0)
    ni=rho/pm.MU;ye=np.einsum('i,zit->zt',pm.ZZ,Y);z2=np.einsum('i,zit->zt',pm.ZZ**2,Y)
    ne=ni*ye;invLD=np.sqrt(QE**2*ni*(ye+z2)/(EPS*KB*T));ae=(3/(4*np.pi*np.maximum(ne,1e-200)))**(1/3);Gammae=QE**2/(4*np.pi*EPS*ae*KB*T)
    total=[];rows=[];stored=[]
    for j,rr in enumerate(pm.RATES):
        rate=sum(np.exp(np.einsum('i,izt->zt',s.a,tf)) for s in rr.sets)
        f=rate*rr.prefactor*(rho/1000)**rr.dens_exp
        for n in rr.reactants:f=f*Y[:,pm.NUCS.index(n),:]
        w=baseweights*f;count=float(w.sum()/ncat);total.append(count)
        if count>1e-10:
            charged=[n.Z for n in rr.reactants if n.Z>0]
            h=np.zeros_like(T);hs=h.copy()
            if len(charged)==2:
                z1,z2c=charged;h=z1*z2c*QE**2*invLD/(4*np.pi*EPS*KB*T)
                hs=.9*Gammae*((z1+z2c)**(5/3)-z1**(5/3)-z2c**(5/3))
            rows.append({'index':j,'reaction':str(rr),'source_label':rr.src,'gross_events_per_initial_catalyst':count,'Q_MeV':float(pm.Q[j]),
                         'flow_weighted_T_GK_5_50_95':weighted_q(T/1e9,w),'flow_weighted_density_kg_m3_5_50_95':weighted_q(rho,w),
                         'flow_weighted_classical_weak_screen_h_5_50_95':weighted_q(h,w),'flow_weighted_ion_sphere_h_5_50_95':weighted_q(hs,w),
                         'flow_weighted_carbon_Gamma_5_50_95':weighted_q(6**(5/3)*Gammae,w),
                         'flow_fraction_T_above_1GK':float(w[T>1e9].sum()/w.sum()),
                         'diagnostic_flow_multiplier_weak':float((w*np.exp(np.clip(h,0,100))).sum()/w.sum()) if len(charged)==2 else 1.})
    total=np.array(total);delta=pm.STOICH@total
    actual=np.sum(m[:,None]*(state[:,:pm.NN,-1]-state[:,:pm.NN,0]),axis=0)/pm.MU/ncat
    neq=float((pm.Q-pm.QNU)@total)
    measured=d['initial_energy_J']*d['deposited_nuclear_over_initial_energy']/pm.MEV/ncat
    nrow=pm.SPECIES.index('n');producers=[rows0 for rows0 in rows if pm.STOICH[nrow,rows0['index']]>0]
    rows.sort(key=lambda x:x['gross_events_per_initial_catalyst'],reverse=True)
    return {'case':record.stem,'source_record':str(record),'record_sha256':sha(record),'trace_sha256':sha(d['trace_path']),
            'Tmax_GK':float(T.max()/1e9),'density_max_kg_m3':float(rho.max()),'initial_catalysts':ncat,
            'neutron_produced_per_initial_catalyst':float(np.maximum(pm.STOICH[nrow],0)@total),
            'neutron_producers':producers,'all_rates':rows,'net_species_change_quadrature_max_absolute_per_initial_catalyst':float(max(abs(delta-actual))),
            'nuclear_deposition_MeV_per_initial_catalyst_from_sampled_flows':neq,'record_nuclear_deposition_MeV_per_initial_catalyst':measured,
            'nuclear_deposition_quadrature_relative':neq/measured-1 if measured else None,
            'screening_scope':'Classical Debye-Huckel and ion-sphere exponents are diagnostics, not a consistent intermediate-coupling correction. Degenerate electron response and Coulomb EOS are omitted here and in the source model.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--record',type=Path,action='append',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();sys.path.insert(0,str(a.root/'docs/protium-fusion'));import pulse_model as pm
    cases=[reduce(p.resolve(),pm) for p in a.record]
    out={'scope':'Trapezoid integration of saved trajectory reaction flows; check quadrature error before quantitative attribution. No modified rates or screening fed back into hydro.','code_sha256':sha(__file__),'cases':cases}
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{k:x[k] for k in ['case','Tmax_GK','density_max_kg_m3','neutron_produced_per_initial_catalyst','net_species_change_quadrature_max_absolute_per_initial_catalyst','nuclear_deposition_quadrature_relative']} for x in cases],indent=2))
if __name__=='__main__':main()
