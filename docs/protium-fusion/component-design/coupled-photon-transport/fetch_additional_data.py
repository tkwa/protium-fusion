from pathlib import Path
import urllib.request,hashlib,tarfile,json,time
base=Path(__file__).resolve().parent;out=[]
items=[('G4PhotonEvaporation.6.1.2.tar.gz','d80ba9bdefcf9e23487a26adfa273304'),('G4PARTICLEXS.4.2.tar.gz','b63341cc97a472b04ff607b094493174')]
for name,md5 in items:
 target=base/name;url='https://cern.ch/geant4-data/datasets/'+name;t=time.monotonic()
 if not target.exists():
  with urllib.request.urlopen(url,timeout=60) as response,target.open('wb') as f:
   while True:
    b=response.read(1024*1024)
    if not b:break
    f.write(b)
 digest=hashlib.md5(target.read_bytes()).hexdigest();assert digest==md5,(name,digest)
 with tarfile.open(target) as tar:tar.extractall(base/'data',filter='data')
 rec={'archive':name,'source':url,'bytes':target.stat().st_size,'md5':digest,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'seconds':time.monotonic()-t};out.append(rec);print(rec,flush=True)
(base/'additional-data-manifest.json').write_text(json.dumps(out,indent=2)+'\n')
