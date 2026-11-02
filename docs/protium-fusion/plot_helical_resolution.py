"""Plot separate integration and matching-scale changes at one source state."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=Path(__file__).parent)
    ap.add_argument('--original-root',type=Path,default=Path(__file__).parent);args=ap.parse_args()
    def read(name,old=False):return json.loads(((args.original_root if old else args.data)/name).read_text())
    old1=[read(f'matched-helical-warm-cut1-{g}.json',True) for g in ['96-64-64','160-128-96']]
    new1=[read(f'matched-moment-warm-cut1-{g}.json') for g in ['256-64-32','384-192-32']]
    aligned=read('matched-aligned-warm-cut1-128-64-64.json')
    cuts=[read(n) for n in ['matched-moment-warm-cut1-384-192-32.json','matched-moment-warm-cut2-768-192-32.json','matched-streaming-warm-cut4-1536-256-32.json','matched-streaming-warm-cut8-3072-512-48.json']]
    fig,ax=plt.subplots(2,2,figsize=(11.8,8.2),layout='constrained')
    a=ax[0,0]
    for rows,label,color in [(old1,'Original pilot (angular grids also change)','#b54a34'),(new1,'Azimuth moments','#226d9b')]:
        a.plot([r['quadrature']['nk'] for r in rows],[r['matched']['perpendicular_heating_eV_s']/1e6 for r in rows],'o-',label=label,color=color)
    a.axhline(aligned['matched']['perpendicular_heating_eV_s']/1e6,color='#40834b',ls='--',label='Independent resonance-aligned integral')
    a.set(title='A. Radial resolution at filter factor 1',xlabel='Radial wave-number nodes',ylabel='Perpendicular heating (MeV/s)')
    a.legend(fontsize=8,loc='lower right')
    a=ax[0,1];c=[r['cut_factor'] for r in cuts]
    a.plot(c,[r['matched']['perpendicular_heating_eV_s']/1e6 for r in cuts],'o-',color='#226d9b')
    a.set(title='B. Resolved short/long matching',xlabel=r'Filter factor $k_{cut}v_{th,Ti}/\Omega_e$',ylabel='Perpendicular heating (MeV/s)',xscale='log')
    a.set_xticks(c,[str(int(x)) for x in c]);a.xaxis.set_minor_formatter(NullFormatter())
    a=ax[1,0]
    keys=['perpendicular_heating_eV_s','parallel_heating_eV_s','total_heating_eV_s']
    labels=['Perpendicular','Parallel','Total']
    for key,label in zip(keys,labels):
        vals=[r['matched'][key] for r in cuts]
        a.plot([2,4,8],[100*abs(vals[i]-vals[i-1])/abs(vals[i]) for i in range(1,4)],'o-',label=label)
    a.set(title='C. Change on doubling the filter scale',xlabel='Larger filter factor',ylabel='Absolute relative change (%)',yscale='log',xscale='log')
    a.set_xticks([2,4,8],['2','4','8']);a.xaxis.set_minor_formatter(NullFormatter());a.legend(fontsize=8)
    a=ax[1,1];pos=np.arange(2)
    for j,(key,label) in enumerate(zip(keys,labels)):
        errors=[]
        for cf,new in zip([2,4],cuts[1:3]):
            old=read(f'matched-helical-warm-cut{cf}-96-64-64.json',True)
            errors.append(100*abs(old['matched'][key]-new['matched'][key])/abs(new['matched'][key]))
        a.bar(pos+(j-1)*.24,errors,width=.22,label=label)
    a.set(title='D. Errors hidden by cancellation in old pilots',xlabel='Filter factor',ylabel='Old pilot error relative to resolved value (%)',yscale='log')
    a.set_xticks(pos,['2','4']);a.legend(fontsize=8)
    for a in ax.flat:a.grid(alpha=.2)
    fig.suptitle('Finite-gyroradius numerical resolution at one 3.468 keV electron state\nMedian Ti source radius; total logarithm 15; static-screened, straight-ion model',fontsize=13)
    for ext in ['png','svg']:
        path=args.data/f'helical-resolution-summary.{ext}'
        fig.savefig(path,dpi=160)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__=='__main__':main()
