"""Summarize independent encounters and unresolved source quadrature/matching."""
from pathlib import Path
import argparse,json
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,default=Path(__file__).parent);args=ap.parse_args();root=args.outdir
    rows=[];checks=[]
    def check(name,ok,metric=None):checks.append({'name':name,'passed':bool(ok),'metric':metric})
    for path in sorted(root.glob('matched-helical-*.json')):
        data=json.loads(path.read_text());rows.append({'file':path.name,**data})
        A=sum(np.array(s['matched']['drift']) for s in data['species']);D=sum(np.array(s['matched']['diffusion']) for s in data['species'])
        check(f'Species drift/covariance sums {path.name}',np.array_equal(A,data['matched']['drift']) and np.array_equal(D,data['matched']['diffusion']))
        check(f'Positive matched covariance {path.name}',min(data['matched']['scaled_covariance_eigenvalues_s-1'])>=0)
        check(f'Positive complementary short interval {path.name}',all(s['short_logarithm']>0 for s in data['species']))
        if data['vperp_m_s']==0 and data['vz_m_s']==0:
            err=abs(data['matched']['parallel_heating_eV_s']/data['unmagnetized']['parallel_heating_eV_s']-1)
            check(f'Cold parallel impulse unaffected by B {path.name}',err<1e-9,err)
    coarse=next(r for r in rows if r['file']=='matched-helical-warm-cut1-96-64-64.json')
    fine=next(r for r in rows if r['file']=='matched-helical-warm-cut1-160-128-96.json')
    metrics=['perpendicular_heating_eV_s','parallel_heating_eV_s','total_heating_eV_s','electron_radial_drift_m_s']
    refinement={k:{'coarse':coarse['matched'][k],'fine':fine['matched'][k],'relative_to_fine':abs(coarse['matched'][k]/fine['matched'][k]-1)} for k in metrics}
    out={'scope':__doc__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,
         'source_cases':rows,'warm_cut1_quadrature_refinement':refinement,
         'source_operator_qualified':False,'reason':'Perpendicular heating changes by about 2.7% under refinement; filter convergence and source-wide resolution remain to be demonstrated.'}
    (root/'finite-gyro-pilot-analysis.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(out['passed'],'/',out['total'],'data and cold-limit checks pass; source operator remains preliminary.')
    print(json.dumps(refinement,indent=2))
    if out['passed']!=out['total']:raise SystemExit(1)

if __name__=='__main__':main()
