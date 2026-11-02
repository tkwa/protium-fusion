"""Charge-flow, electrostatic and gated-loss screens for the finite end system.

The source is the explicit 140 T / 70% / 99.9% network plant. Outgoing reaction
products are fully stripped. Prescribed equal energy per charge sets their
relativistic transport speeds. Uniform cylindrical space charge is a diagnostic
for a long unscreened beam; it is not a finite electrode or kinetic solution.
"""
from pathlib import Path
from collections import defaultdict
from math import pi,sqrt,log,exp
import json,argparse
from catalyst_network import load_data,stage,HOT,label
from catalyst_catalog import nucmass
from nonneutral_screen import MEV,C,EPS0,E_CHARGE
ROOT=Path(__file__).parent
ME=9.1093837139e-31


def speed(za,voltage):
    z,a=za;m=nucmass(z,a)*MEV/C**2
    u=z*E_CHARGE*voltage/(m*C*C)
    return C*sqrt(u*(2+u))/(1+u)


def streams(plant,energy,data):
    rows=[];aggregate=defaultdict(float)
    for za in HOT:
        key=label(za);flow=energy['network']['processor_pn_rates'][key]
        if not flow:continue
        src=plant['selected_source_profiles'][key]
        geom=next(r for r in plant['channel_allocations'] if r['ZA']==list(za))
        nuclides=defaultdict(float);s=stage(za,data[za],energy['T_MeV'])
        for e in s['events']:
            if e['mt']==103:continue
            nuclides[tuple(e['residual'])]+=e['rate']
            for p,n in e['light'].items():
                if p[0]>0:nuclides[p]+=e['rate']*n
        charge=sum(p[0]*v for p,v in nuclides.items())
        expected=(za[0]+1)*next(r for r in energy['stages'] if tuple(r['ZA'])==za)['extraction_events_per_pn']
        products=[]
        for product,f in sorted(nuclides.items(),key=lambda x:-x[1]):
            rate=plant['primary_pn_rate_s']*flow*f;aggregate[product]+=rate
            products.append({'ZA':list(product),'per_processor_pn':f,'rate_s':rate,'current_A':rate*product[0]*E_CHARGE})
        current=sum(p['current_A'] for p in products)
        rows.append({'processor':key,'products':products,'charge_per_pn':charge,'charge_ledger_error_per_pn':charge-expected,
          'current_A':current,'channels':geom['channels'],'current_per_channel_A':current/geom['channels'],
          'aperture_current_density_A_m2':current/(geom['channels']*pi*src['radius_m']**2)})
    return rows,[{'ZA':list(p),'rate_s':r,'current_A':p[0]*E_CHARGE*r} for p,r in aggregate.items()]


def transport_case(rows,aggregate,primary_rate,voltage=40e6,beam_to_wall_ratio=1.,potential_limit=4e4):
    # phi(axis)-phi(wall) = lambda/(4*pi*eps0)*(1+2 ln(b/a)).
    factor=1+2*log(1/beam_to_wall_ratio);line=0.;electron_K=0.;curr=0.;details=[]
    for p in aggregate:
        za=tuple(p['ZA']);v=speed(za,voltage);kinetic_ratio=za[0]*voltage/(1e6*nucmass(*za))
        line+=p['current_A']/v;curr+=p['current_A']
        electron_K+=p['rate_s']*za[0]*kinetic_ratio*ME*C*C/MEV/primary_rate
    for row in rows:
        lam=sum(p['current_A']/speed(tuple(p['ZA']),voltage) for p in row['products'])/row['channels']
        phi=lam/(4*pi*EPS0)*factor
        # All product species in this group are at least as fast as the slowest.
        slowest=min(row['products'],key=lambda p:speed(tuple(p['ZA']),voltage))
        vmin=speed(tuple(slowest['ZA']),voltage)
        # Sufficient planar bound: require delta_phi <= beta*U; use the
        # corresponding minimum speed to bound rho everywhere in the gap.
        beta=.1;j=row['aperture_current_density_A_m2']
        safe_gap=sqrt(8*EPS0*speed(tuple(slowest['ZA']),voltage*(1-beta))*beta*voltage/j)
        details.append({'processor':row['processor'],'channel_potential_V':phi,
          'maximum_identical_channels_per_long_beam_at_limit':potential_limit/phi,
          'slowest_product_ZA':slowest['ZA'],'slowest_product_speed_m_s':vmin,
          'sufficient_planar_gap_m_for_10percent_spacecharge':safe_gap})
    return {'energy_per_charge_V':voltage,'beam_radius_over_wall_radius':beam_to_wall_ratio,
      'summed_outgoing_current_A':curr,'summed_line_charge_C_m':line,
      'long_merged_beam_axis_potential_V':line/(4*pi*EPS0)*factor,
      'co_moving_electron_kinetic_MeV_per_primary':electron_K,
      'neutralization_electrons_per_primary':curr/E_CHARGE/primary_rate,
      'channel_details':details,'allowed_long_beam_potential_V':potential_limit}


def gates(plant,energy):
    L=plant['active_length_m'];T=energy['T_MeV'];mp=nucmass(1,1)*MEV/C**2;out=[]
    for key,src in plant['selected_source_profiles'].items():
        rate=src['reaction_rate_per_m_s']*L;flux=2*src['proton_inventory_per_m']*sqrt(T*MEV/(2*pi*mp))
        open_loss=flux/rate
        out.append({'processor':key,'open_H_end_flux_s_per_channel':flux,
          'processor_pn_s_per_channel':rate,'open_H_losses_per_processor_pn':open_loss,
          'closed_H_losses_at_25kT_per_processor_pn':open_loss*exp(-25),
          'gate_duty_for_one_H_loss_per_processor_pn':1/open_loss,
          'catalyst_consumption_lifetime_s':src['catalyst_inventory_per_m']/src['reaction_rate_per_m_s'],
          'proton_consumption_lifetime_s':src['proton_inventory_per_m']/src['reaction_rate_per_m_s'],
          'daughter_to_catalyst_number_for_1ms_residence':.001*src['reaction_rate_per_m_s']/src['catalyst_inventory_per_m']})
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--outdir',type=Path,default=ROOT);args=ap.parse_args()
    root=args.root;outdir=args.outdir
    energy=json.loads((root/'network-ion-budget.json').read_text())
    plant=next(r for r in json.loads((root/'network-plant.json').read_text())['cases'] if r['assumptions']['B']==140 and r['assumptions']['heat']==.7 and r['assumptions']['end_in']==.999)
    rows,agg=streams(plant,energy,load_data())
    cases=[transport_case(rows,agg,plant['primary_pn_rate_s'],voltage=U,potential_limit=min(4e4,.1*U)) for U in [1e4,4e4,1e5,1e6,40e6,55e6]]
    out={'scope':__doc__,'reference_primary_rate_s':plant['primary_pn_rate_s'],'processor_streams':rows,'aggregate_products':agg,'transport_cases':cases,'axial_gate_screen':gates(plant,energy),
      'limits':['No selective excitation, collector, electron suppression or beam neutralization design is established.',
      'Current counts isotope-changing reaction products only. Nonreacting leakage and wall-ion circulation are separate.',
      'Space-charge potential assumes a uniform long cylinder with the beam filling the conducting wall.',
      'Comoving electron kinetic energy excludes creation, phase-space matching, recombination, instability and collection costs.',
      'Planar bound assumes prescribed electrode potentials and positive-energy transmitted particles; transverse optics are separate.']}
    margin=json.loads((root/'network-plant-sensitivity.json').read_text())['additional_electrical_MeV_per_primary_at_1mW_kg']
    nch=sum(r['channels'] for r in rows);scale=plant['primary_pn_rate_s']*MEV
    out['power_allowances']={'total_channels':nch,'extra_electrical_MeV_per_primary_at_1mW_kg':margin,
      'extra_electrical_W_at_1mW_kg':margin*scale,'extra_electrical_W_per_channel_at_1mW_kg':margin*scale/nch,
      'already_reserved_extraction_W':plant['source_loss_MeV_per_primary']['extraction_loss_MeV_per_pn']*scale,
      'already_reserved_atomic_preparation_W':plant['source_loss_MeV_per_primary']['atomic_preparation_reserve_MeV_per_pn']*scale,
      'already_reserved_general_auxiliary_W':plant['other_auxiliary_W']}
    out['input_files']=['network-ion-budget.json','network-plant.json','network-plant-sensitivity.json','beam-data/tendl2025-proton-recycling.zip']
    (outdir/'end-transport-screen.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for d in cases:print('U',d['energy_per_charge_V'],'current MA',d['summed_outgoing_current_A']/1e6,'merged phi GV',d['long_merged_beam_axis_potential_V']/1e9,'electron MeV',d['co_moving_electron_kinetic_MeV_per_primary'],'Ti',d['channel_details'][0])
    print('Ti gate',out['axial_gate_screen'][0])
