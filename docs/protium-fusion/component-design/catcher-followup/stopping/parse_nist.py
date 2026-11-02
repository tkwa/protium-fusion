"""Extract graphite H and He evaluated tables from the archived primary PDF text."""
from pathlib import Path
import json,re,hashlib
P=Path(__file__).resolve().parent/'sources'
def main():
 text=(P/'nist-newstar-2017.txt').read_text();out={}
 for sp,number in [('p',8),('he4',11)]:
  block=text.split(f'Table A.{number}.',1)[1].split(f'Table A.{number+1}.',1)[0];rows=[]
  for line in block.splitlines():
   fields=line.split()
   if len(fields)!=6:continue
   try:v=list(map(float,fields))
   except ValueError:continue
   if v[0]<=0 or v[0]>10000:continue
   rows.append(v)
  assert len(rows)>45 and all(rows[i+1][0]>rows[i][0] for i in range(len(rows)-1))
  assert all(0<v[5]<=1 and abs((v[1]+v[2])/v[3]-1)<.001 for v in rows)
  out[sp]={'table':f'A.{number}','columns':['energy_MeV','electronic_MeV_cm2_g','elastic_MeV_cm2_g','total_MeV_cm2_g','CSDA_g_cm2','projected_to_CSDA'],'rows':rows}
 assert next(v for v in out['p']['rows'] if v[0]==1)[4]==.002803
 assert next(v for v in out['he4']['rows'] if v[0]==.5)[4]==.0003292
 source=json.loads((P/'nist-source.json').read_text());source.update(tables=out,parser_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),visual_check='Printed pages10 and13 rendered and inspected; table headings, units and relevant low-MeV entries verified. The extracted60MeV alpha detour0.9905 is preserved as published text; it lies outside all incident alpha energies in this calculation.')
 (P/'nist-graphite-tables.json').write_text(json.dumps(source,indent=2)+'\n');print({k:len(v['rows']) for k,v in out.items()})
if __name__=='__main__':main()
