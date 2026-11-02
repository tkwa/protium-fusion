"""Propagate separately refined B2 traces through the rounded station account."""
import hashlib,json
from pathlib import Path
from station import evaluate
P=Path(__file__).resolve().parent
rows={}
for key,stem in [('64','station-inputs'),('128','station-inputs-b2-128-energy-matched'),('256','station-inputs-b2-256-energy-matched')]:
    path=P/(stem+'.json'); x=json.loads(path.read_text()); c=x['cycle']; r=x['receiver']; R=1.1e9
    f=.8/(r['occupied_s_per_He_at_R1e9']['mass999']*R/1e9)
    d=evaluate(f,R,c,r,absorber_sigma=5000,driver_peak_W_kg=100,service_sigma=1000)
    rows[key]={'input':path.name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'result':d}
out={'scope':'B2-only spatial sensitivity at the rounded operating point; the other three stages remain64zones.', 'station_code_sha256':hashlib.sha256((P/'station.py').read_bytes()).hexdigest(),'cases':rows}
(P/'station-grid-sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
