"""Extract published LUNA repository tables from archived web-tool responses."""
import hashlib,json,re,pathlib
HERE=pathlib.Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
s=json.loads((HERE/'luna-web-tables.json').read_text())
parts=s.split('13C(p,g)14N table page')
blocks=[parts[0],parts[1]+json.loads((HERE/'luna-web-tail.json').read_text())]
out={'revision':'luna2023-public-repository-2026-10-31','paper_doi':'10.1103/PhysRevLett.131.162701',
 'scope':'Authors/collaboration Bayesian evaluated thermonuclear rates, including 0.007–10 GK. Newer C13 2025 evaluation reports agreement through 10 GK; its numerical supplement was not retrieved. No independent extrapolation above 1 GK is used.',
 'columns':['T9','minus3sigma','minus2sigma','minus1sigma','median','plus1sigma','plus2sigma','plus3sigma'],
 'units':{'T9':'GK','rate':'cm^3 mol^-1 s^-1'},
 'uncertainty_note':'Repository labels are pointwise 1,2,3 sigma bands. Their cross-temperature correlation is not supplied; coherent low/high runs are sensitivities, not independent posterior samples.',
 'retrieval_date':'2026-10-31','retrieval_sha256':{p.name:sha(p) for p in [HERE/'luna-web-tables.json',HERE/'luna-web-tail.json',HERE/'luna2023-source.tar']},'tables':{}}
for name,block,id_ in zip(['C12_p_to_N13_reaclib','C13_p_to_N14_reaclib'],blocks,[178,180]):
 rows={}
 for line in block.splitlines():
  if not re.match(r'L\d+: \d',line):continue
  nums=re.findall(r'\d\.\d{3}e[+-]\d{2}',line)
  if len(nums)==8:
   row=list(map(float,nums));rows[row[0]]=row
 vals=[rows[t] for t in sorted(rows)]
 assert len(vals)==54 and vals[0][0]==.007 and vals[-1][0]==10,(name,len(vals))
 assert all(all(a<b for a,b in zip(r[1:],r[2:])) for r in vals)
 out['tables'][name]={'url':f'https://luna.lngs.infn.it/index.php/scientific-output/reaction-rate-repository%3Flayout%3Dedit%26id%3D{id_}.html','rows':vals}
(HERE.parent/'updated-carbon-rates.json').write_text(json.dumps(out,indent=2)+'\n')
print('Both tables: 54 temperatures, seven quantiles, 0.007–10 GK.')
