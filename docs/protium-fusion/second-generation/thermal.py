#!/usr/bin/env python3
"""Finite graphite storage and isolated surface-response screens for new pulses.
The all-phase bound adds independent extrema. Prompt extrapolation is a proxy.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import erfc
P=Path(__file__).resolve().parent;F=P.parent
sys.path.insert(0,str(F));from fast_eos import FastEOS
R=1.1e9;C=299792458.;AREA=4*np.pi*R**2;RHO=2000.;CP=1800.;SIGMA=350.;L=SIGMA/RHO
M=json.loads((P/'rates.json').read_text());names=M['species'];mass=dict(zip(names,M['mass_MeV']))
MU=1.66053906660e-27;MEV=1.602176634e-13
parent=F/'component-design/receiver_thickness.py'
tree=ast.parse(parent.read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['pulse_bins','front_step']],type_ignores=[]),str(parent),'exec'))
rates={'n13':.001159113256783153,'o14':.009817978825996571,'o15':.005681556796002065}
daughters={'n13':'c13','o14':'n14','o15':'n15'};nu={'n13':.706,'o14':1.05,'o15':1.0}
eos=FastEOS()

def source(path):
    d=json.loads(path.read_text());trace=Path(d['trace_path']) if 'trace_path' in d else path.with_suffix('.npz')
    tr=np.load(trace);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);st=tr['state'].reshape(nz,nv,-1)
    nn=d['state_variables'].index('T/T0');m=tr['cell_mass_kg'];A=np.array(M['A'][:27]);Z=np.array(M['Z'][:27])
    ro=np.r_[0.,st[:,nn+1,-1]*d['R0_m']];v=np.r_[0.,st[:,nn+2,-1]*float(tr['velocity_unit_cm_s'])/100]
    y=np.maximum(st[:,:27,-1],0.);bar=y@A;rho=m/(4*np.pi/3*np.diff(ro**3))
    u=eos.state(rho/1000,st[:,nn,-1]*d['T0_K'],(y@Z)/bar,y.sum(axis=1)/bar).e*1e-4
    if np.any(np.diff(v)<=0):raise ValueError('Endpoint is not yet ballistic: '+str(path))
    slope=m/np.diff(v);drdv=np.diff(ro)/np.diff(v);b=R-ro[:-1]+drdv*v[:-1]
    off=float(tr['time_s'][-1])-drdv-R/299792458.
    dec={n:float(m@y[:,names.index(n)])/MU*(mass[n]-mass[daughters[n]]-nu[n])*MEV/AREA for n in rates}
    K=float(np.sum(slope*np.diff(v**3)/6))/AREA;U=float(m@u)/AREA
    return dict(m=m,v=v,b=b,off=off,slope=slope,u=u,decay_fluence=dec,
        fluence=K+U+d['radiation_escaped_J']/AREA+sum(dec.values()),
        prompt_fluence=d['radiation_escaped_J']/AREA,K=K,U=U,record=d,
        record_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest())

def main():
    inp=json.loads((F/'component-design/final-family-dispatch/sigma350/input.json').read_text())
    st={n:source(Path(v['source'])) for n,v in inp['loaded_stages'].items()}
    new=source(P/'pulse-n15-m100-v1-z128-t1.json')
    reference=st['n15'];st['n15']=new
    dt=.5;edges=np.arange(0,20000+dt,dt);time=edges[1:];rows=[]
    period={n:inp['plant']['pulse_interval_s'][n]*(100 if n=='n15' else 1)/16 for n in st}
    q={}; saved={}
    for n,s in st.items():
        K,U,dec,momentum=pulse_bins(s,edges,True)
        q[n]=K+U+dec;saved[n+'_material_W_m2']=K+U;saved[n+'_decay_W_m2']=dec
    qmean=sum(s['fluence']/period[n] for n,s in st.items())
    storage=sum(s['fluence'] for s in st.values())/(SIGMA*CP)
    delay=8951.35*qmean/(SIGMA*CP)
    for k in [50.,200.,500.]:
        ker=np.diff(front_step(edges,L,k,CP,'N'));individual={}
        for n,s in st.items():
            curve=fftconvolve(q[n],ker)[:len(time)]
            individual[n]=dict(nonuniform_material_decay_peak_K=float(max(curve)),
                peak_heat_W_m2=float(max(q[n])),total_fluence_J_m2=s['fluence'],
                material_energy_in_first20000_s_J_m2=float(np.sum(q[n])*dt))
        lower=1600+storage+delay+qmean*L/(2*k)
        rows.append(dict(k_W_m_K=k,bulk_storage_swing_bound_K=storage,dispatch_lateness_storage_K=delay,
            mean_back_pickup_minimum_offset_K=qmean*L/(2*k),
            material_phase_all_extrema_envelope_K=lower+sum(v['nonuniform_material_decay_peak_K'] for v in individual.values()),
            isolated=individual,through_depth_diffusion_s=L*L*RHO*CP/k,
            fundamental_nonuniform_decay_s=L*L*RHO*CP/k/math.pi**2))
    # Same photon spectrum/deposition profile, pulse duration proportional to source radius.
    ratio=new['prompt_fluence']/reference['prompt_fluence'];width=100**(1/3)
    prompt=dict(old_fluence_J_m2=reference['prompt_fluence'],new_fluence_J_m2=new['prompt_fluence'],
        fluence_ratio=ratio,assumed_width_ratio=width,
        old_finite_depth_flash_rise_K=54.81,
        same_spectrum_diffusive_rise_proxy_K=54.81*ratio/math.sqrt(width),
        qualification='Scaling proxy only: deltaT proportional to fluence/sqrt(width) at the same spectral depth and material diffusivity. New multigroup photon transport and source angular/time histories must replace this proxy.')
    np.savez_compressed(P/'thermal-waveforms.npz',time_s=time,**saved)
    out=dict(scope=__doc__,rows=rows,prompt=prompt,all_source_mean_W_m2=qmean,
        source_hashes={n:dict(record=s['record_sha256'],trace=s['trace_sha256']) for n,s in st.items()},
        limitations=['Constant properties and smooth ballistic endpoints. The all-source bound assigns missed-hot photons to hot storage, conservatively.',
            '1600 K is the minimum useful pickup, not the maximum exposed surface temperature.',
            'Skin/post resistance, damaged conductivity, joints, creep, carbon emission and material lifetime require additional design.',
            'Isolated surface maxima are summed for an all-phase screen; these are not a simulated observed maximum in the new calendar.',
            'Existing NRT damage laws and carbon emission fits cannot establish material survival at these temperatures.'])
    (P/'thermal-results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(prompt=prompt,rows=[{k:r[k] for k in r if k!='isolated'} for r in rows]),indent=2))

if __name__=='__main__':main()
