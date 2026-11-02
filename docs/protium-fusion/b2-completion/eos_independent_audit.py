"""Independent B2 legacy-tail EOS, energy-scale and viscosity-sign audit.

Run with PYTHONPATH=/tmp/protium-python and the repository's recorded Python
numerical dependencies. Uses frozen inputs only; writes eos-independent-audit.json
beside this script. It neither evolves a pulse nor modifies any source EOS.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'phase-two/pulse-transport')]
import pynucastro as pn
from fast_eos import FastEOS, KB_CGS, NA, ARAD_CGS
from microphysics import coulomb, E, MU, KB, HB, ME, C
import radial_pulse as rp


def main():
    old=ROOT/'reassessment-2049/small-stages'
    record_path=old/'scaled-hot-B2-64-R30-tight.json'
    trace_path=old/'scaled-hot-B2-64-R30-tight-compact.npz'
    record=json.loads(record_path.read_text())
    trace=np.load(trace_path)
    species=trace['species'].tolist(); nn=len(species); nv=int(trace['variables_per_zone'])
    nz=int(trace['zones']); aa=np.array([pn.Nucleus(n).A for n in species]); zz=np.array([pn.Nucleus(n).Z for n in species])
    mass_si=record['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3)
    eos=FastEOS(); maximum=None; visited=0
    for j,time in enumerate(trace['time_s']):
        state=trace['state'][:,j].reshape(nz,nv)
        radii_si=np.r_[0,state[:,nn+1]*record['R0_m']]
        rho_si=mass_si/(4*np.pi/3*np.diff(radii_si**3))
        temperature=state[:,nn]*record['T0_K']
        y=np.maximum(state[:,:nn],0);y/=np.sum(y*aa,axis=1)[:,None]
        ye=y@zz;yi=y.sum(axis=1);rho=rho_si/1000
        therm=eos.state(rho,temperature,ye,yi)
        residual=(therm.p-rho*rho*therm.de_drho-temperature*therm.dp_dT)/(temperature*therm.dp_dT)
        for i in range(nz):
            row=dict(time_s=float(time),zone=i,T_K=float(temperature[i]),rho_kg_m3=float(rho_si[i]),relative_thermal_identity_error=float(residual[i]))
            if maximum is None or abs(residual[i])>abs(maximum['relative_thermal_identity_error']):maximum=row
            visited+=1
    # All variables now describe the final archived accepted endpoint.
    comp=pn.Composition(species)
    for i,n in enumerate(species):comp[n]=float(y[-1,i]*aa[i])
    comp.normalize();direct=pn.StellarEOS();spot=[]
    for temp in [1e7,3e6,1e6,float(temperature[-1]),1.00001e5,7.97065182e4,1e4]:
        s=direct.pe_state(float(rho[-1]),temp,comp)
        thermal=s.p/rho[-1]**2-s.de_drho
        row=dict(rho_cgs=float(rho[-1]),T_K=temp,source_eta=float(s.eta),source_p_cgs=float(s.p),source_e_erg_g=float(s.e),source_cv=float(s.c_v),source_adiabatic_coefficient=float(thermal),source_identity_relative_to_thermal=float((s.p-rho[-1]**2*s.de_drho-temp*s.dp_dT)/(temp*s.dp_dT)))
        if temp>=1e5:
            f=eos.state(rho[-1],temp,comp.ye,1/comp.abar)
            coefficient=f.p/rho[-1]**2-f.de_drho;maxwell=temp*f.dp_dT/rho[-1]**2
            row.update(fast_adiabatic_vs_source=float(coefficient/thermal-1),fast_identity_relative_to_thermal=float((coefficient-maxwell)/maxwell),fast_cv_vs_source=float(f.c_v/s.c_v-1),fast_p_vs_source=float(f.p/s.p-1),fast_e_vs_source=float(f.e/s.e-1))
        spot.append(row)
    uc,pc,gamma,ge=coulomb(temperature,rho_si,y,zz)
    # Independent isolated-atom hydrogenic upper scale, n=1 then n=2 shell
    # occupancy. Omits repulsion and all dense-plasma shifts; not actual binding.
    atomic_upper_eV=13.605693123*zz**2*(np.minimum(zz,2)+np.maximum(zz-2,0)/4)
    binding_upper=E/MU*(y@atomic_upper_eV)
    cold=[]
    for i in np.where(temperature<1e8)[0]:
        ne=rho_si[i]/MU*ye[i];pf=HB*(3*np.pi**2*ne)**(1/3)
        ef=ME*C*C*(np.sqrt(1+(pf/(ME*C))**2)-1)
        cold.append(dict(zone=int(i),mass_fraction=float(mass_si[i]/record['mass_kg']),T_K=float(temperature[i]),rho_kg_m3=float(rho_si[i]),fermi_energy_eV=float(ef/E),fermi_temperature_K=float(ef/KB),dominant_ion_Gamma=float(gamma[i]),ion_ion_pressure_fraction=float(pc[i]/(therm.p[i]/10)),ideal_internal_J=float(mass_si[i]*therm.e[i]*1e-4),ion_ion_energy_correction_J=float(mass_si[i]*uc[i]),isolated_binding_upper_J=float(mass_si[i]*binding_upper[i]),fraction_initial_ideal_internal=float(mass_si[i]*therm.e[i]*1e-4/record['initial_energy_J']),fraction_initial_abs_ion_ion=float(abs(mass_si[i]*uc[i])/record['initial_energy_J']),fraction_initial_binding_upper=float(mass_si[i]*binding_upper[i]/record['initial_energy_J']),fraction_nuclear_binding_upper=float(mass_si[i]*binding_upper[i]/record['nuclear_energy_J'])))
    # Reconstruct the legacy scalar viscosity and discrete spherical work.
    initial=eos.state(record['rho0_kg_m3']/1000,record['T0_K'],8/14,3/14)
    cs0=float(np.sqrt(initial.gamma1*initial.p/(record['rho0_kg_m3']/1000)))
    radii=radii_si*100;velocity=np.r_[0,state[:,nn+2]*cs0];mass=mass_si*1000
    compress=np.maximum(-np.diff(velocity),0)
    sound=np.sqrt(np.maximum(therm.gamma1*therm.p/rho,0))
    q=rho*(record['artificial_viscosity_quadratic']*compress**2+record['artificial_viscosity_linear']*sound*compress)
    volume_dot=4*np.pi*np.diff(radii*radii*velocity);rho_dot=-rho*rho*volume_dot/mass
    adiabatic=(therm.p/rho**2-therm.de_drho)*rho_dot/therm.c_v
    viscous=q*rho_dot/(rho*rho*therm.c_v)
    corrected_q=np.where(volume_dot<0,q,0.)
    assert np.all(-corrected_q*volume_dot/mass>=0)
    opacity=record['opacity_m2_kg']*10;er=ARAD_CGS*temperature**4
    rc=.5*(radii[1:]+radii[:-1]);rf=.5*(rho[1:]+rho[:-1]);ef=.5*(er[1:]+er[:-1]);grad=np.diff(er)/np.diff(rc)
    ratio=np.abs(grad)/(opacity*rf*np.maximum(ef,1e-300));lim=(2+ratio)/(6+3*ratio+ratio**2)
    flux=-C*100*lim*grad/(opacity*rf)
    flux_outer=C*100*er[-1]/(4+3*opacity*rho[-1]*(radii[-1]-rc[-1]))
    luminosity=np.r_[0,4*np.pi*radii[1:-1]**2*flux,4*np.pi*radii[-1]**2*flux_outer]
    qrad=np.diff(luminosity)/mass
    # Floor here is a diagnostic to suppress known invalid cold fit extrapolations;
    # no such modification is proposed as a measured rate or hydro result.
    dy,qtot,qnu=rp.network(np.maximum(temperature,1e7),rho,state[:,:nn])
    e_yi=1.5*KB_CGS*NA*temperature
    e_ye=(therm.e+rho*therm.de_drho-e_yi*yi)/ye
    rhs=dict(zone=nz-1,T_K=float(temperature[-1]),rho_cgs=float(rho[-1]),rho_dot_cgs_s=float(rho_dot[-1]),volume_dot_cm3_s=float(volume_dot[-1]),outer_minus_inner_velocity_cm_s=float(np.diff(velocity)[-1]),legacy_q_cgs=float(q[-1]),q_over_physical_pressure=float(q[-1]/therm.p[-1]),adiabatic_Tdot_K_s=float(adiabatic[-1]),legacy_viscous_Tdot_K_s=float(viscous[-1]),radiative_Tdot_K_s=float(-qrad[-1]/therm.c_v[-1]),floor_diagnostic_nuclear_Tdot_K_s=float((qtot[-1]-qnu[-1])/therm.c_v[-1]),floor_diagnostic_composition_Tdot_K_s=float((-e_yi[-1]*dy[-1].sum()-e_ye[-1]*(dy[-1]@zz))/therm.c_v[-1]),compression_only_q_cgs=float(corrected_q[-1]),compression_only_work_nonnegative_all_cells=True)
    inputs=[record_path,trace_path,ROOT/'radial_pulse.py',ROOT/'pulse_model.py',ROOT/'fast_eos.py',ROOT/'electron-eos.npz',ROOT/'phase-two/pulse-transport/microphysics.py',Path(__file__)]
    out=dict(scope=__doc__,visited_saved_cell_states=visited,max_thermal_identity_residual=maximum,source_EOS_spot_checks=spot,cold_endpoint_scales=cold,outer_endpoint_rhs=rhs,input_sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},limitations=['The old endpoint is affected by negative artificial viscous dissipation; its cold-state numbers are diagnostics, not a corrected trajectory.','Ideal source EOS assumes complete ionization; source agreement is not validation of dense cold matter.','Coulomb estimate retains only the ion-ion rigid-background term. Atomic upper is an isolated-atom scale. Neither bounds trajectory changes or opacity.'])
    (HERE/'eos-independent-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('visited_saved_cell_states','max_thermal_identity_residual','cold_endpoint_scales','outer_endpoint_rhs')},indent=2))

if __name__=='__main__':main()
