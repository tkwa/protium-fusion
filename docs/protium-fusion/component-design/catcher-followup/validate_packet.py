#!/usr/bin/env python3
"""Read-only aggregate checks for the catcher follow-up; writes its own report.
These checks establish file consistency, numerical identities and provenance.
They do not turn conditional material properties into measured performance.
"""
from pathlib import Path
import ast,hashlib,json,math,re,subprocess
P=Path(__file__).resolve().parent;R=P.parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=[];failures=[]
def check(name,ok,detail=None):
    checks.append({'name':name,'passed':bool(ok),**({'detail':detail} if detail is not None else {})})
    if not ok:failures.append(name)
def read(p):return json.loads(Path(p).read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Non-JSON constant '+x)))
for f in sorted(P.rglob('*.py')):
    try:ast.parse(f.read_text());check('python_syntax:'+str(f.relative_to(P)),True)
    except Exception as e:check('python_syntax:'+str(f.relative_to(P)),False,str(e))
for f in sorted(P.rglob('*.json')):
    if f.name in ['validation.json','manifest.json'] and f.parent==P:continue
    try:read(f);check('strict_json:'+str(f.relative_to(P)),True)
    except Exception as e:check('strict_json:'+str(f.relative_to(P)),False,str(e))
for f in sorted(P.rglob('manifest.json')):
    if f.parent==P:continue
    d=read(f);k=next((k for k in ['sha256','files','files_sha256']if k in d),None)
    if k is None:check('manifest_schema:'+str(f.relative_to(P)),False);continue
    entries=d[k]
    if isinstance(entries,list):entries={x['name']:x['sha256']for x in entries}
    for name,val in entries.items():
        target=f.parent/name;expected=val['sha256']if isinstance(val,dict) else val
        check('component_hash:'+str(target.relative_to(P)),target.is_file() and sha(target)==expected)
# Resolve simple Markdown local links; external URLs, anchors and code templates
# remain under the source-specific evidence review.
for f in sorted(P.rglob('*.md')):
    for target in re.findall(r'\]\(([^\s]+?)\)',f.read_text()):
        if '://'in target or target.startswith('#') or target.startswith('mailto:'):continue
        target=target.split('#')[0].strip('<>')
        if not target or any(x in target for x in ['{','}','*','`']):continue
        q=Path(target) if target.startswith('/') else f.parent/target
        check('markdown_target:'+str(f.relative_to(P))+':'+target,q.exists())
sp=read(P/'stopping/surface-spectra.json');rf=read(P/'refurbishment/results.json');pt=read(R/'mixture-design/selected-plant/point.json')['plant']
A=4*math.pi*(1.1e9)**2;mu=1.66053906660e-27
check('sphere_area',abs(sp['area_m2']/A-1)<1e-12)
check('backing_areal_mass',sp['slab_sigma_kg_m2']==350.)
check('backing_depth',abs(sp['slab_thickness_m']-.175)<1e-14)
for species,mass,field in [('p',1.00782503223,'neutral_H_kg_s'),('he4',4.00260325413,'neutral_He_kg_s')]:
    rate=sum(v['species'][species]['nuclei_m2']/v['period_s']for v in sp['stages'].values())
    check('independent_spectrum_flow:'+species,abs(rate*A*mu*mass/rf[field]-1)<1e-12)
check('returned_H_matches_plant',abs(rf['neutral_H_kg_s']/pt['recovered_H_kg_s']-1)<1e-9)
check('exact_dispatch_count',sum(rf['dispatch_370day_events'].values())==4228)
pa=read(P/'plant-account/results.json')
for row in pa['rows']:
    if 'plant'not in row:check('explicit_infeasible:'+row['label'],row['feasible']is False and bool(row.get('reasons')));continue
    d=row['plant'];s=row['scenario'];c=d['catcher_named_account'];a=d['area_m2']
    check('plant_mass:'+row['label'],abs(sum(d['mass_breakdown_kg'].values())/d['mass_kg']-1)<1e-10)
    check('plant_electric:'+row['label'],abs((d['gross_electric_W']-d['driver_W']-d['ordinary_auxiliary_W']-d['pump_W']-d['replacement_W']-d['support_W']-d['net_W'])/d['gross_electric_W'])<1e-10)
    check('plant_cold:'+row['label'],abs(sum(d['cold_heat_components_W'].values())/d['cold_rejection_W']-1)<1e-10)
    check('named_aux_fits:'+row['label'],(c['already_named_electrical_W_m2']+c['new_named_electrical_W_m2'])*a<=d['ordinary_auxiliary_W']*(1+1e-10))
    check('held_material_is_not_fresh_debit:'+row['label'],abs(d['fresh_C_kg_s']/pt['fresh_C_kg_s']-1)<1e-10 if 'fresh_C_kg_s'in d else d['helium_per_s']==pt['helium_per_s'])
# This is a preservation check on the inherited archive, independent of whether
# the new packet has already been staged or committed.
repo=P.parents[3]
proc=subprocess.run(['git','diff','--name-only','fb26ac5bc08ee4880daf7f9a815d67b4787aac64','--','docs/protium-fusion'],cwd=repo,text=True,capture_output=True)
if proc.returncode==0:
    outside=[x for x in proc.stdout.splitlines()if not x.startswith('docs/protium-fusion/component-design/catcher-followup/')]
    check('preexisting_research_unchanged',not outside,outside)
else:check('preexisting_research_git_check_available',False,proc.stderr.strip())
report={'scope':__doc__,'checks':checks,'passed':sum(x['passed']for x in checks),'failed':len(failures),'failures':failures,'code_sha256':sha(__file__)}
(P/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':report['passed'],'failed':len(failures),'failures':failures},indent=2))
raise SystemExit(bool(failures))
