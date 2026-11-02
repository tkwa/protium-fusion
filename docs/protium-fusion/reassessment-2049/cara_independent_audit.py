"""Independent arithmetic audit of saved CARA and selected-plant results.

Python 3 standard library only. Run beside the audited JSON files. Reconstructs
selected equations without importing either generating script; writes
cara-independent-audit.json. Does not validate a physical accelerator design.
"""
from pathlib import Path
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MEV = 1.602176634e-13
E = 1.602176634e-19
M = 9.1093837139e-31
C = 299792458.0
EPS = 8.8541878188e-12
checks = []

def check(name, x, y, tolerance=2e-12):
    error = abs(x-y)/max(abs(x),abs(y),1e-100)
    assert error <= tolerance, (name,x,y,error)
    checks.append(dict(name=name,relative_error=error))

def main():
    d = json.loads((HERE/'ion-cara-screen.json').read_text())
    p = json.loads((HERE/'selected-ion-point.json').read_text())
    k = d['kinematics'][1]
    z = d['selected_parametric_case']
    ref = next(x for x in json.loads((ROOT/'network-plant.json').read_text())['cases']
               if x['assumptions']['B']==140 and x['assumptions']['heat']==.6
               and x['assumptions']['end_in']==.999)
    mc2=M*C*C/MEV; omega=E*k['B_T']/M; a=E*k['E_V_m']/(M*C)
    u=math.sqrt(2*k['target_kinetic_MeV']/mc2)
    retarded=u/a; lab=retarded*(1+u*u/6); distance=C*(lab-retarded)
    pulse=EPS*k['E_V_m']**2*math.pi*k['radius_m']**2*C*retarded
    cool=3*math.pi*EPS*M**3*C**3/(E**4*k['B_T']**2)
    # Simpson integral of P(tau)*(dt/dtau), instead of the saved closed form.
    n=10000; step=retarded/n
    integrand=lambda t:mc2/(2*cool)*(a*t)**2*(1+(a*t)**2/2)
    radiation=(integrand(0)+integrand(retarded)+sum(
        (4 if i%2 else 2)*integrand(i*step) for i in range(1,n)))*step/3
    for name,independent,saved in [
        ('retarded pulse',retarded,k['retarded_pulse_s']),
        ('particle transit',lab,k['particle_transit_s']),
        ('acceleration distance',distance,k['acceleration_length_m']),
        ('empty pulse energy',pulse,k['empty_pulse_J']),
        ('synchrotron integral',radiation,k['first_order_synchrotron_MeV'])]:
        check(name,independent,saved)
    areas=[x['channels']*math.pi*x['radius_m']**2 for x in ref['channel_allocations']]
    ti=next(x for x in ref['channel_allocations'] if x['ZA']==[22,49])
    factor=sum(areas)/(ti['channels']*math.pi*ti['radius_m']**2)
    check('all-processor guide factor',factor,z['network_guide_factor'])
    rate=d['primary_rate_per_active_Ti_m_s']; b=z['pair_supply_per_primary']; cap=z['capture_fraction_per_complete_attempt']
    A=factor*2*pulse*(.01+1/1000)/.5/(rate*MEV)
    B=b*100*(1/.95-.6)*(1/cap-.5)
    f=math.sqrt(B/A); work=b*120/.5; total=2*math.sqrt(A*B)+work
    check('RF coefficient A',A,z['A_MeV_s']); check('population coefficient B',B,z['B_MeV_per_s'])
    check('optimized attempt rate',f,z['optimized_attempts_s']); check('total RF envelope',total,z['total_added_MeV_per_primary'])
    # Explicit deterministic birth-and-removal recurrence after 1000 periods.
    births=b*rate; period=1/f; population=0.
    for _ in range(1000):
        integrated_population=population*period+.5*births*period**2
        population=(population+births*period)*(1-cap)
    residence=integrated_population/period/births
    check('discrete population residence',residence,z['mean_susceptible_positron_residence_s'])
    width=2*a/(omega*u); bands=[]
    for row in d['monochromatic_bandwidth_discriminator']:
        T=row['illustrative_axial_temperature_eV']
        sigma=math.sqrt(T/(mc2*1e6)); fraction=min(1.,2*math.erf(width/(math.sqrt(2)*sigma)))
        Bw=b*100*(1/.95-.6)*(1/fraction-.5)
        cost=2*math.sqrt(A*Bw)+work
        check('band fraction '+str(T),fraction,row['maximum_fraction_two_optimally_placed_bands'])
        check('band budget '+str(T),cost,row['corrected_budget']['total_added_MeV_per_primary'])
        bands.append(dict(kT_eV=T,two_band_fraction_upper=fraction,cost_MeV_primary=cost))
    # Independently reconstruct closures in all 20 selected-plant comparison rows.
    sa=p['retained_plant_assumptions']; year=365.25*86400
    for i,row in enumerate(p['cases']):
        mass=sum(row['mass_breakdown_kg'].values())
        check('case mass '+str(i),mass,row['total_mass_kg'])
        maintenance=sa['annual_replacement_fraction']*sa['embodied_J_kg']*mass/year
        check('case manufacturing '+str(i),maintenance,row['replacement_manufacturing_W'])
        check('case electrical '+str(i),row['net_W']+maintenance,row['bus_before_replacement_W'])
        check('case internal heat '+str(i),row['internal_waste_W'],row['nuclear_after_capture_W']-row['net_W'])
        check('case external heat '+str(i),row['radiator_area_m2']*row['radiator_emission_flux_W_m2'],row['internal_waste_W']+row['external_absorbed_W'])
        check('case specific power '+str(i),1000*row['net_W']/mass,row['net_mW_kg'])
        check('case processing power '+str(i),p['source_rate_s']*MEV*row['extra_processing_MeV_per_primary'],row['extra_processing_W'])
    hashes=[]
    for label,data in [('CARA',d),('selected plant',p)]:
        for filename,expected in data['input_sha256'].items():
            actual=hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()
            assert actual==expected,(label,filename,'stale input')
            hashes.append(dict(packet=label,path=filename,sha256=actual))
    selected=p['selected_conservative_case']
    out=dict(scope=__doc__,check_count=len(checks),all_checks_passed=True,
        max_relative_error=max(x['relative_error'] for x in checks),checks=checks,
        checked_input_hashes=hashes,
        audited_packet_sha256={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in
            [HERE/'ion-cara-screen.py',HERE/'ion-cara-screen.json',HERE/'selected_ion_point.py',HERE/'selected-ion-point.json']},
        independent_values=dict(retarded_s=retarded,lab_s=lab,distance_m=distance,pulse_J=pulse,
            radiation_MeV=radiation,guide_factor=factor,A_MeV_s=A,B_MeV_s_inv=B,
            attempt_rate_s_inv=f,mean_residence_s=residence,cost_MeV_primary=total,
            fractional_detuning_halfwidth=width,illustrative_bands=bands),
        selected_plant={key:selected[key] for key in ('net_W','total_mass_kg','net_mW_kg','extra_processing_W','construction_energy_payback_years')},
        limitations='Arithmetic and provenance audit only. No kinetic capture, radiation-reaction phase, RF source, segmentation, isotope-processing, or manufactured hardware validation.')
    (HERE/'cara-independent-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('check_count','all_checks_passed','max_relative_error','independent_values','selected_plant')},indent=2))

if __name__=='__main__':main()
