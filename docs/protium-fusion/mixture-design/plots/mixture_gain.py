"""Export the declared 64-zone mixture sweep and median-carbon128 comparison.

Default plotting uses the compact CSV and provenance JSON only. --rebuild-data
re-reads the exact run directories, checks constraints, and rewrites compact data.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, sys
P=Path(__file__).resolve().parent
ROOT=P.parent.parent
STAGES=['A','B1','B2','C13']
CAT={'A':'n14','B1':'n15','B2':'c12','C13':'c13'}
BASE_H={'A':1.,'B1':1.,'B2':2.,'C13':1.}
SERIES=['frozen64_equal_energy','frozen64_fixed_flow','median_carbon128_equal_energy']

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def near(a,b):return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-15)

def rebuild(runroot, rateroot):
    baseline_files={'A':'regression-A-compression_cap.json','B1':'regression-B1-compression_cap.json','B2':'B2-64-compression-cap.json','C13':'regression-C13-compression_cap.json'}
    baselines={s:json.loads((ROOT/'b2-completion'/f).read_text()) for s,f in baseline_files.items()}
    records=[];ignored=[];by_key={};provenance={}
    directories=[(runroot/x,64) for x in ['control','grid','refine']]+[(rateroot/'nominal128',128)]
    for directory,zones in directories:
        for path in sorted(directory.glob('*.json')):
            if any(path.stem.endswith('.'+suffix) for suffix in ['initial','inputs','receiver','pretrace']):continue
            d=json.loads(path.read_text())
            if not isinstance(d,dict) or 'mixture_inputs' not in d:continue
            i=d['mixture_inputs'];s=d['stage'];duration=i.get('duration_s',i.get('physical_duration_s'))
            conditions={'stage':s in STAGES,'zones':d['zones']==zones,'opacity':near(d['opacity_m2_kg'],.001),'temperature_factor':i.get('temperature_factor',1.)==1.,'prepared_energy_factor':i.get('prepared_energy_factor',1.)==1.,'completed':d['success'] and d['termination_reason']=='final_time' and near(d['elapsed_s'],duration),'reference_duration':near(duration,baselines[s]['elapsed_s'])}
            ru=d.get('rate_update',i.get('rate_update'))
            if zones==64:
                conditions['rate_family']=ru is None
                conditions['constraint']=i['constraint'] in ['fixed_energy','fixed_flow']
                series=SERIES[0] if i['constraint']=='fixed_energy' else SERIES[1]
            else:
                conditions['rate_family']=bool(ru and ru.get('mode')=='updated' and ru.get('quantile')=='median' and ru.get('revision')=='luna2023-public-repository-2026-10-31')
                conditions['constraint']=i['constraint']=='fixed_energy'
                series=SERIES[2]
            if not all(conditions.values()):
                ignored.append({'path':str(path),'failed':{k:v for k,v in conditions.items() if not v}});continue
            assert d['viscosity_mode']=='compression_cap' and d['radiation_diffusion']
            assert not d['reduced_three_species_A_network']
            assert near(d['mass_kg'],baselines[s]['mass_kg'])
            h=i['hydrogen_per_catalyst'];y=d['initial_Y']
            assert near(y['p']/y[CAT[s]],h)
            assert near(d['deposited_nuclear_over_initial_energy'],(d['nuclear_energy_J']-d['weak_energy_escaped_upper_J'])/d['initial_energy_J'])
            if i['constraint']=='fixed_energy':assert near(d['initial_energy_J'],baselines[s]['initial_energy_J'])
            else:assert near(i['velocity_scale'],1.)
            key=(s,series,h)
            if key in by_key:raise ValueError('Duplicate case '+str(key))
            row={'stage':s,'series':series,'H_per_catalyst':h,'gain':d['deposited_nuclear_over_initial_energy'],'zones':zones,'recorded_constraint':i['constraint'],'shared_reference':False,'mass_kg':d['mass_kg'],'initial_energy_J':d['initial_energy_J'],'nuclear_energy_J':d['nuclear_energy_J'],'modeled_weak_escape_J':d['weak_energy_escaped_upper_J'],'deposited_nuclear_energy_J':d['nuclear_energy_J']-d['weak_energy_escaped_upper_J'],'duration_s':duration,'elapsed_s':d['elapsed_s'],'opacity_m2_kg':d['opacity_m2_kg'],'T0_K':d['T0_K'],'rho0_kg_m3':d['rho0_kg_m3'],'velocity_scale':i['velocity_scale'],'temperature_factor':i.get('temperature_factor',1.),'prepared_energy_factor':i.get('prepared_energy_factor',1.),'source_record':str(path),'source_sha256':sha(path)}
            by_key[key]=row;records.append(row)
            provenance[str(path)]={'sha256':sha(path),'profile_source':i['profile_source'],'profile_source_sha256':i['profile_source_sha256'],'baseline_source':i['baseline_source'],'baseline_source_sha256':i['baseline_source_sha256'],'solver_source_sha256':i['solver_source_sha256'],'network_source_sha256':i['network_source_sha256'],'runner_source_sha256':i['runner_source_sha256'],'rate_update':ru,'rtol':d['rtol'],'atol':d['atol']}
    for s in STAGES:
        key=(s,SERIES[1],BASE_H[s])
        if key not in by_key:
            original=by_key[(s,SERIES[0],BASE_H[s])]
            assert near(original['velocity_scale'],1.)
            row=dict(original,series=SERIES[1],shared_reference=True)
            records.append(row);by_key[key]=row
    records.sort(key=lambda x:(STAGES.index(x['stage']),SERIES.index(x['series']),x['H_per_catalyst']))
    counts={s:{f:sum(r['stage']==s and r['series']==f for r in records) for f in SERIES} for s in STAGES}
    expected={'A':[11,4,2],'B1':[9,4,3],'B2':[9,4,2],'C13':[10,4,3]}
    assert {s:[counts[s][f] for f in SERIES] for s in STAGES}==expected,counts
    with (P/'mixture-gain.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    out={'scope':__doc__,'input_directories':[str(x[0]) for x in directories],'plotted_rows':len(records),'unique_source_records':len(provenance),'counts':counts,'shared_reference_rows':4,'excluded_records':ignored,'baseline_durations_s':{s:b['elapsed_s'] for s,b in baselines.items()},'baseline_record_sha256':{s:sha(ROOT/'b2-completion'/f) for s,f in baseline_files.items()},'comparison':{'all_series':'Prescribed spherical hot inflow, fixed pulse mass and enclosed-mass density/temperature profiles within each stage; fixed physical duration; compression-cap viscosity; full 27-isotope/194-rate network; grey opacity 0.001 m2/kg.','frozen64_equal_energy':'Frozen nuclear rates, 64 zones; total prepared energy fixed by adjusting radial velocity.','frozen64_fixed_flow':'Frozen nuclear rates, 64 zones; radial velocity fixed, so total prepared energy varies with composition. Shared baseline rows reuse the identical equal-energy input at H/catalyst 1,1,2,1.','median_carbon128_equal_energy':'128 zones and the public median 12C(p,gamma)13N / 13C(p,gamma)14N rate updates together; total prepared energy fixed. These points change both numerical resolution and the two rate pairs relative to the faint curve.','gain':'Recorded nuclear release minus modeled weak-energy escape, divided by initial internal-plus-bulk-kinetic prepared energy; (nuclear_energy_J - weak_energy_escaped_upper_J) / initial_energy_J. Subsequent radiative escape is not subtracted from this deposited-at-production gain. It is a hot-pulse diagnostic, not a whole-plant electrical gain.','line_scope':'Straight segments join sampled 64-zone cases as a visual guide. No smoothed optimum or unrun case is inferred.','uncertainty':'No error bars: these datasets are different modeling choices, not draws defining a confidence or credible interval. Mesh convergence and source formation remain separate studies.'},'source_records':provenance,'csv_sha256':sha(P/'mixture-gain.csv'),'builder_sha256':sha(__file__)}
    (P/'mixture-gain-provenance.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    return out

def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FixedFormatter, MultipleLocator
    provenance=json.loads((P/'mixture-gain-provenance.json').read_text())
    assert sha(P/'mixture-gain.csv')==provenance['csv_sha256']
    with (P/'mixture-gain.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
    for r in rows:
        for key in ['H_per_catalyst','gain']:r[key]=float(r[key])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.labelcolor':'#253247','text.color':'#253247','axes.edgecolor':'#bbc3cc','xtick.color':'#46556a','ytick.color':'#46556a','axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'protium-mixture-gain-2026-10-31','pdf.fonttype':42,'savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(11.8,9.6),sharex=True,sharey=True)
    fig.subplots_adjust(left=.105,right=.975,bottom=.225,top=.785,wspace=.12,hspace=.25)
    titles={'A':'A · nitrogen-14','B1':'B1 · nitrogen-15','B2':'B2 · carbon-12 return','C13':'C13 · carbon-13'}
    styles={SERIES[0]:dict(color='#9ca8b5',marker='o',markersize=3.6,linewidth=1.5,alpha=.9,zorder=2),SERIES[1]:dict(color='#ba7434',marker='^',markersize=5.2,linewidth=1.3,linestyle=(0,(4,3)),markerfacecolor='white',markeredgewidth=1.25,zorder=3),SERIES[2]:dict(color='#165a95',marker='o',markersize=7.0,linestyle='None',markeredgecolor='white',markeredgewidth=.7,zorder=5)}
    for ax,s in zip(axes.flat,STAGES):
        ax.set_xscale('log',base=2);ax.set_xlim(.45,18);ax.set_ylim(-.5,27.5)
        ax.xaxis.set_major_locator(FixedLocator([.5,1,2,4,8,16]));ax.xaxis.set_major_formatter(FixedFormatter(['½','1','2','4','8','16']))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.grid(axis='y',color='#e4e9ee',linewidth=.7,zorder=0)
        ax.axvline(BASE_H[s],color='#d5dce3',linestyle=(0,(2,3)),linewidth=1,zorder=1)
        for series in SERIES:
            subset=sorted([r for r in rows if r['stage']==s and r['series']==series],key=lambda r:r['H_per_catalyst'])
            ax.plot([r['H_per_catalyst'] for r in subset],[r['gain'] for r in subset],**styles[series])
        ax.set_title(titles[s],loc='left',fontsize=12.2,fontweight='bold',pad=11)
        duration=provenance['baseline_durations_s'][s]*1e6
        ax.text(1.,1.035,f'{duration:.3f} μs',transform=ax.transAxes,ha='right',va='bottom',fontsize=9.2,color='#647287')
        ax.tick_params(axis='both',which='major',length=3,pad=5)
    fig.text(.105,.95,'Pulse gain versus fuel mixture',fontsize=20,fontweight='bold',ha='left')
    fig.text(.105,.915,'Completed spherical-pulse simulations at prescribed hot inflow',fontsize=11.5,color='#536378',ha='left')
    legend=[Line2D([0],[0],**styles[SERIES[0]],label='64 zones · frozen rates · equal energy'),Line2D([0],[0],**styles[SERIES[1]],label='64 zones · frozen rates · fixed flow'),Line2D([0],[0],**styles[SERIES[2]],label='128 zones · median carbon rates · equal energy'),Line2D([0],[0],color='#bcc6d1',linestyle=(0,(2,3)),linewidth=1,label='Reference H/catalyst: A 1 · B1 1 · B2 2 · C13 1')]
    fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(.097,.892),frameon=False,ncol=2,fontsize=10.1,handlelength=2.7,columnspacing=2.4,labelspacing=.75)
    fig.text(.54,.178,'Hydrogen nuclei per catalyst nucleus (log₂ scale)',ha='center',fontsize=11.5)
    fig.text(.029,.53,'Deposited nuclear energy / prepared energy',rotation=90,va='center',ha='center',fontsize=11.5)
    footer=['Mass, density and temperature profiles are fixed within each stage; grey opacity κ = 0.001 m²/kg.',
            'Equal energy adjusts inflow velocity. Fixed flow retains velocity and changes prepared energy. Physical duration is fixed.',
            'The 128-zone points change mesh and two carbon capture rates together. Lines connect tested 64-zone cases.',
            'Fixed-flow reference points reuse identical baseline inputs. Pulse gain precedes whole-plant energy accounting.']
    for j,line in enumerate(footer):fig.text(.105,.132-.025*j,line,ha='left',fontsize=9.1,color='#58677b')
    fig.savefig(P/'mixture-gain.png',dpi=210)
    fig.savefig(P/'mixture-gain.pdf',metadata={'Title':'Pulse gain versus fuel mixture','Author':'Protium fusion research','CreationDate':None,'ModDate':None})
    fig.savefig(P/'mixture-gain.svg',metadata={'Title':'Pulse gain versus fuel mixture','Creator':'Protium fusion research','Date':None})
    plt.close(fig)
    manifest={'script_sha256':sha(__file__),'data_sha256':sha(P/'mixture-gain.csv'),'provenance_sha256':sha(P/'mixture-gain-provenance.json'),'matplotlib_version':matplotlib.__version__,'python_version':sys.version,'exports':{f:sha(P/f) for f in ['mixture-gain.png','mixture-gain.pdf','mixture-gain.svg']},'point_counts':provenance['counts']}
    (P/'mixture-gain-exports.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'unique_source_records':provenance['unique_source_records'],'exports':list(manifest['exports'])},indent=2))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--rebuild-data',action='store_true');parser.add_argument('--run-root',type=Path,default=Path('/home/tkwa/src/protium-fusion/output/mixture-runs'));parser.add_argument('--rate-root',type=Path,default=Path('/home/tkwa/src/protium-fusion/output/mixture-rate-runs'));args=parser.parse_args()
    if args.rebuild_data:rebuild(args.run_root,args.rate_root)
    plot()
if __name__=='__main__':main()
