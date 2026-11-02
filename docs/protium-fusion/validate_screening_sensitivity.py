"""Verify consistent logarithmic removal and independent distant-limit response."""
from pathlib import Path
from math import log
import argparse,json,time
import numpy as np
from scipy.special import exp1
from guiding_center_collisions import coefficients,species_table
from electron_wall_orbits import ME,E_CHARGE as EV


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).parent);ap.add_argument('--data',type=Path,default=Path(__file__).parent);args=ap.parse_args();start=time.monotonic()
    ref=json.loads((args.root/'matched-streaming-warm-cut2-512-192-32.json').read_text());data=json.loads((args.data/'screening-sensitivity.json').read_text());checks=[]
    def check(name,value,tol):checks.append({'name':name,'passed':bool(value<tol),'metric':float(value),'upper_limit':tol})
    v=ref['species'][1]['vth_m_s'];S=np.diag([1/v,1/v,1/v**2])
    def exact_log(factor):
        a=(factor*ref['kappa_m-1']/ref['kcut_m-1'])**2
        return .5*((1+a)*np.exp(a)*exp1(a)-1)
    for cold in [True,False]:
        group=[r for r in data['cases'] if r['cold']==cold];base=next(r for r in group if r['kappa_multiplier']==1)
        vz=0 if cold else ref['vz_m_s'];gc=coefficients(vz,ref['local_state']['relative_azimuthal_m_s'],species_table(ref['local_state']['species'],1.))
        Agc=np.array([gc[0],gc[1],0.]);Dgc=np.array([[gc[3],gc[5],0.],[gc[5],gc[4],0.],[0.,0.,0.]])
        for row in group:
            label=f'cold={cold}, factor={row["kappa_multiplier"]}, nk={row["nk"]}'
            eig=np.linalg.eigvalsh(S@np.array(row['diffusion'])@S)
            check('Positive covariance '+label,-eig.min()/eig.max(),1e-13)
            predicted_change=-log(row['kappa_multiplier'])-exact_log(row['kappa_multiplier'])+exact_log(1.)
            actual_change=row['species'][0]['complementary_log']-base['species'][0]['complementary_log']
            check('Close-range logarithm follows exact finite-filter correction '+label,abs(actual_change-predicted_change),1e-11)
            check('Filtered logarithm matches independent analytic integral '+label,abs(row['species'][0]['filtered_log']-exact_log(row['kappa_multiplier'])),1e-11)
            if row['kappa_multiplier']>1:
                for key,expected in [('drift',Agc),('diffusion',Dgc)]:
                    actual=np.array(base[key])-np.array(row[key]);expected=expected*log(row['kappa_multiplier'])
                    if key=='drift':actual=S@actual;expected=S@expected
                    else:actual=S@actual@S;expected=S@expected@S
                    error=np.linalg.norm(actual-expected)/np.linalg.norm(expected)
                    check('Removed distant interval matches analytic GC '+key+' '+label,error,.005)
        a=next(r for r in group if r['kappa_multiplier']==10 and r['nk']==512);b=next(r for r in group if r['kappa_multiplier']==10 and r['nk']==768)
        for key in ['perpendicular_heating_eV_s','parallel_heating_eV_s','total_heating_eV_s','electron_radial_drift_m_s']:
            check(f'Radial refinement cold={cold}, '+key,abs(a[key]-b[key])/max(abs(b[key]),1e-100),1e-10)
        if not cold:
            for key in ['perpendicular_heating_eV_s','parallel_heating_eV_s','total_heating_eV_s','electron_radial_drift_m_s']:
                check('Reproduce prior warm factor-2 source '+key,abs(base[key]-ref['matched'][key])/abs(ref['matched'][key]),1e-12)
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,'seconds':time.monotonic()-start}
    (args.data/'screening-sensitivity-validation.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:out[k] for k in ['passed','total','seconds']}))
    for c in checks:
        if not c['passed']:print('FAILED',c)
    if out['passed']!=out['total']:raise SystemExit(1)


if __name__=='__main__':main()
