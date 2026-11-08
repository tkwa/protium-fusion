"""Extract growth histories and conservation from axisymmetric Athena++ output."""
import argparse,json,re,sys
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import legval

BASE=Path('/var/tmp/tkwa-chatgpt/protium-stability')
sys.path.insert(0,str(BASE/'athena/vis/python'))
import athena_read

def snapshot(files):
    blocks=[]
    for f in files:
        r,t,p,d=athena_read.vtk(f)
        blocks.append((t[0],r,t,{k:v[0] for k,v in d.items()}))
    blocks.sort(key=lambda b:b[0])
    rr=blocks[0][1]
    tt=np.concatenate([b[2][:-1] for b in blocks]+[blocks[-1][2][-1:]])
    data={k:np.concatenate([b[3][k] for b in blocks],axis=0) for k in blocks[0][3]}
    with open(files[0],'rb') as f:
        head=f.read(300).decode('ascii',errors='replace')
    match=re.search(r'time=\s*([\deE+.-]+)',head)
    return float(match.group(1)),rr,tt,data

def metrics(files,ell):
    t,rf,tf,d=snapshot(files)
    rho=d['dens'];fuel=d['s1'];prod=d['s0'];mom=d['mom']
    dv=np.diff(rf**3)[None,:]/3 * (np.cos(tf[:-1])-np.cos(tf[1:]))[:,None]*2*np.pi
    weights=(np.cos(tf[:-1])-np.cos(tf[1:]))/2
    rcent=.75*np.diff(rf**4)/np.diff(rf**3)
    mu=.5*(np.cos(tf[:-1])+np.cos(tf[1:]))
    coeff=np.zeros(ell+1);coeff[-1]=1;pol=legval(mu,coeff)
    pol-=weights@pol;pol/=np.sqrt(weights@(pol*pol))
    mass=(fuel*dv).sum()
    rhobar=weights@rho
    fbar=weights@fuel
    variance=weights@((rho/rhobar-1)**2)
    radial_weight=fbar*np.diff(rf**3)
    density_rms=np.sqrt(radial_weight@variance/radial_weight.sum())
    velocity=mom/rho[...,None]
    vmean=(weights[:,None]*fuel*velocity[:,:,0]).sum(axis=0)/np.maximum(fbar,1e-99)
    vr_scatter=np.sqrt(np.sum(fuel*dv*(velocity[:,:,0]-vmean)**2)/mass)
    theta_ke=np.sum(.5*rho*dv*velocity[:,:,1]**2)
    radii={}
    for q in [.25,.5,.75,.9]:
        rr=[]
        for row in fuel:
            shellmass=row*np.diff(rf**3);cdf=np.r_[0,np.cumsum(shellmass)]
            # Uniform density within each finite-volume cell; interpolate r^3.
            rr.append(np.interp(q*cdf[-1],cdf,rf**3)**(1/3))
        rr=np.array(rr);mean=weights@rr
        radii[str(q)]=dict(mean_m=float(mean*30),rms_fraction=float(np.sqrt(weights@((rr/mean-1)**2))),
                           mode_fraction=float(weights@((rr/mean-1)*pol)))
    return dict(time_us=t*30,fuel_mass=float(mass),product_mass=float((prod*dv).sum()),
                density_rms=float(density_rms),vr_scatter_m_s=float(vr_scatter*1e6),
                theta_ke=float(theta_ke),fuel_mean_density_kg_m3=float(np.sum(fuel*dv*rho)/mass*1e8),
                radii=radii,max_density_kg_m3=float(rho.max()*1e8))

def analyze(case):
    cfg=json.loads((case/'config.json').read_text())
    h=athena_read.hst(case/'pulse.hst')
    qheat=float((Path(cfg['inputs'])/'rates.txt').read_text().split()[0])
    rows=[]
    groups={}
    for f in case.glob('pulse.block*.out2.*.vtk'):
        number=f.name.split('.')[-2];groups.setdefault(number,[]).append(f)
    if cfg['nt']>1 and cfg.get('geometry','spherical')=='spherical':
        for number,files in sorted(groups.items()):rows.append(metrics(files,cfg['ell']))
    initial=h['E_minus_nuclear'][0]
    clock_offset=cfg.get('reference_time_us',11.) if cfg.get('perturbation')==5 else 0.
    late_start=np.interp(max(h['time'][0],h['time'][-1]-2/30),h['time'],h['product_mass'])
    report=dict(name=case.name,config=cfg,execution=json.loads((case/'execution.json').read_text()) if (case/'execution.json').exists() else None,
        final_time_us=float(h['time'][-1]*30),conversion=float(h['product_mass'][-1]/h['fuel_mass'][0]),
        reference_clock_offset_us=clock_offset,reference_clock_final_time_us=float(h['time'][-1]*30+clock_offset),
        initial_product_fraction=float(h['product_mass'][0]/h['fuel_mass'][0]),
        output_fraction_produced_in_last_2us=float((h['product_mass'][-1]-late_start)/max(h['product_mass'][-1],1e-99)),
        gain=float(qheat*h['product_mass'][-1]/initial),initial_energy_J=float(initial*2.7e24),
        nuclear_output_J=float(qheat*h['product_mass'][-1]*2.7e24),
        initial_mass_kg=float(h['fuel_mass'][0]*2.7e12),
        final_fuel_mass_change=float(h['fuel_mass'][-1]/h['fuel_mass'][0]-1),
        energy_residual_fraction_input=float(h['E_minus_nuclear'][-1]/initial-1),
        maximum_energy_residual_fraction_input=float(np.max(abs(h['E_minus_nuclear']/initial-1))),
        history={k:h[k].tolist() for k in ['time','fuel_mass','product_mass','fuel_rho','fuel_r2','theta_ke','E_minus_nuclear']},
        snapshots=rows)
    if 'EOS_outside_fuel_mass' in h:
        report['maximum_EOS_outside_fuel_fraction']=float(np.max(h['EOS_outside_fuel_mass'])/h['fuel_mass'][0])
        report['peak_temperature_K']=float(np.max(h['fuel_Tmax']))
    if 'EOS_face_limits' in h:
        report['EOS_face_limits']=float(h['EOS_face_limits'][-1])
    if 'EOS_dense_face_limits' in h:
        report['EOS_dense_face_limits']=float(h['EOS_dense_face_limits'][-1])
    if cfg.get('geometry')=='cartesian':
        momentum=np.column_stack([h[f'{axis}-mom'] for axis in [1,2,3]])
        report['initial_com_speed_m_s']=float(np.linalg.norm(momentum[0])/h['mass'][0]*1e6)
        report['maximum_vector_momentum_change_per_initial_mass_m_s']=float(np.max(np.linalg.norm(momentum-momentum[0],axis=1))/h['mass'][0]*1e6)
    for row in rows:
        row['theta_ke_fraction_input']=row['theta_ke']/initial
        row['reference_clock_time_us']=row['time_us']+clock_offset
    (case/'diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['name','final_time_us','conversion','gain','initial_energy_J','energy_residual_fraction_input']}))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('names',nargs='+');a=p.parse_args()
    for name in a.names:analyze(BASE/'runs'/name)
