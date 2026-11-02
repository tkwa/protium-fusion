"""Cross-version rate comparison. Differences are not statistical errors."""
from pathlib import Path
import json,gzip
from catalyst_network import load_data
from catalyst_catalog import thermal
from catalyst_thermal import CURVES
ROOT=Path(__file__).parent

if __name__=='__main__':
    now=load_data();old=json.loads(gzip.decompress((ROOT/'recycling-comparison-input.json.gz').read_bytes()))
    rows=[]
    for el,z,a in [('Ti',22,49),('V',23,51),('Cr',24,52),('Mn',25,53)]:
        o=next(d for d in old if d['source_file']==f'p-{el}{a:03d}-mt.tendl')
        for version,d in [('TENDL-2023 alternate MT format',o),('TENDL-2025 standard',now[(z,a)])]:
            cs={int(k):v for k,v in d['MF3'].items()};sv=thermal(cs[4],1.6,z,a)
            rows.append({'target':f'{a}{el}','version':version,'T_MeV':1.6,'pn_sv_m3_s':sv,
              'ratios_to_total_pn':{str(mt):thermal(cs[mt],1.6,z,a)/sv for mt in [102,107,16,28,104,111] if mt in cs},
              'source_file':d['source_file'],'source_sha256':d['source_sha256']})
        padf=CURVES.get(f'p{z:02d}{a:03d}.ndf')
        if padf:
            chans=padf['channels'];ground=thermal(chans[f'{(z+1)*1000+a}_0'],1.6,z,a)
            total=sum(thermal(t,1.6,z,a) for key,t in chans.items() if key.startswith(f'{(z+1)*1000+a}_'))
            rows.append({'target':f'{a}{el}','version':'PADF-2 2026','T_MeV':1.6,
              'pn_ground_sv_m3_s':ground,'pn_all_listed_states_sv_m3_s':total,
              'source_file':padf['source_file'],'source_sha256':padf['source_sha256']})
    (ROOT/'recycling-data-comparison.json').write_text(json.dumps({'scope':__doc__,'cases':rows},indent=2)+'\n')
    print('comparison rows',len(rows))
