"""Real Geant4 option4 electromagnetic receiver transport, native energy scorers.

No photonuclear physics. Closed cavity returns are tracked by native geometry.
Escaped calorimetric energy is K for e- and K+2mec2 for e+, preserving positron
annihilation potential. Single-purpose process exits after validation/persistence,
avoiding the observed Python/Geant4 destruction-order crash after completed runs.
Cuts are applied through Geant4 UI because inherited Python physics-list methods
misaddress the base object in this wheel; actual material-couple tables are logged.
"""
from pathlib import Path
import argparse,json,math,os,time,hashlib,statistics,collections,sys
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--normal-exit',action='store_true');ap.add_argument('--sigma',type=float,default=1000);ap.add_argument('--radius-m',type=float,default=1e6);ap.add_argument('--depth-um',type=float,default=25);ap.add_argument('--mode',choices=['isotropic','radial'],default='isotropic');ap.add_argument('--energies',default='.511');ap.add_argument('--histories',type=int,default=100);ap.add_argument('--cut-um',type=float,default=10);ap.add_argument('--seed',type=int,default=86137);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--data',type=Path,default=P/'data');a=ap.parse_args()
for k,v in {'GEANT4_DATA_DIR':'','G4LEDATA':'G4EMLOW8.8','G4ENSDFSTATEDATA':'G4ENSDFSTATE3.0','G4LEVELGAMMADATA':'PhotonEvaporation6.1.2','G4PARTICLEXSDATA':'G4PARTICLEXS4.2'}.items():os.environ[k]=str(a.data/v)
import geant4_pybind as g
start=time.monotonic();R=a.radius_m*g.m;thick=a.sigma/2000*g.m;Ro=R+thick
class Detector(g.G4VUserDetectorConstruction):
 def __init__(self):super().__init__();self.sds=[];self.scorers=[]
 def Construct(self):
  n=g.G4NistManager.Instance();vac=n.FindOrBuildMaterial('G4_Galactic');carbon=n.FindOrBuildElement('C');self.mat=g.G4Material('ReceiverGraphite',2*g.g/g.cm3,1,g.kStateSolid,1000*g.kelvin);self.mat.AddElement(carbon,1)
  self.world=g.G4LogicalVolume(g.G4Orb('World',Ro+10*g.m),vac,'World');phys=g.G4PVPlacement(None,g.G4ThreeVector(),self.world,'World',None,False,0,False)
  self.shell=g.G4LogicalVolume(g.G4Sphere('Shell',R,Ro,0,2*math.pi,0,math.pi),self.mat,'Shell');g.G4PVPlacement(None,g.G4ThreeVector(),self.shell,'Shell',self.world,False,0,False);return phys
 def ConstructSDandField(self):
  for name,lv in [('Hot',self.shell),('Vacuum',self.world)]:
   sd=g.G4MultiFunctionalDetector(name);sc=g.G4PSEnergyDeposit('Edep');sd.RegisterPrimitive(sc);g.G4SDManager.GetSDMpointer().AddNewDetector(sd);self.SetSensitiveDetector(lv,sd);self.sds.append(sd);self.scorers.append(sc)
class Event(g.G4UserEventAction):
 def __init__(self):super().__init__();self.rows=[];self.ids=None;self.particles=collections.Counter();self.creators=collections.Counter()
 def BeginOfEventAction(self,event):self.eventid=event.GetEventID();self.leak=0.;self.escapes={};self.steps=0;self.tracks=0;self.residual=0.
 def EndOfEventAction(self,event):
  if self.ids is None:self.ids=[g.G4SDManager.GetSDMpointer().GetCollectionID(name+'/Edep') for name in ['Hot','Vacuum']]
  hce=event.GetHCofThisEvent();values=[]
  for i in self.ids:
   hc=hce.GetHC(i);values.append(sum(v for _,v in hc) if hc is not None else 0.)
  self.rows.append({'deposited_MeV':values[0]/g.MeV,'vacuum_deposit_MeV':values[1]/g.MeV,'escaped_calorimetric_MeV':self.leak/g.MeV,'residual_inside_kinetic_MeV':self.residual/g.MeV,'steps':self.steps,'tracks':self.tracks,'escape_by_particle_MeV':self.escapes})
class Generator(g.G4VUserPrimaryGeneratorAction):
 def __init__(self):super().__init__();self.gun=g.G4ParticleGun(1);self.gun.SetParticleDefinition(g.G4Gamma.GammaDefinition());self.energy=.511
 def GeneratePrimaries(self,event):
  offset=a.depth_um*g.um if a.depth_um>0 else -1*g.um;self.gun.SetParticlePosition(g.G4ThreeVector(0,0,R+offset));self.gun.SetParticleEnergy(self.energy*g.MeV);mu=2*g.G4UniformRand()-1 if a.mode=='isotropic' else 1.;phi=2*math.pi*g.G4UniformRand();q=math.sqrt(max(0,1-mu*mu));self.gun.SetParticleMomentumDirection(g.G4ThreeVector(q*math.cos(phi),q*math.sin(phi),mu));self.gun.GeneratePrimaryVertex(event)
class Tracking(g.G4UserTrackingAction):
 def __init__(self,event):super().__init__();self.e=event
 def PreUserTrackingAction(self,tr):
  self.e.tracks+=1;self.e.particles[tr.GetParticleDefinition().GetParticleName()]+=1;p=tr.GetCreatorProcess();self.e.creators[p.GetProcessName() if p else 'primary']+=1
 def PostUserTrackingAction(self,tr):
  e=self.e;e.steps+=tr.GetCurrentStepNumber();particle=tr.GetParticleDefinition().GetParticleName();K=tr.GetKineticEnergy()
  if tr.GetPosition().mag()>Ro:
   val=K+(2*g.electron_mass_c2 if particle=='e+' else 0.);e.leak+=val;e.escapes[particle]=e.escapes.get(particle,0)+val/g.MeV
  else:e.residual+=K
class DiagnosticStep(g.G4UserSteppingAction):
 def __init__(self):super().__init__()
 def UserSteppingAction(self,step):
  if event.eventid!=752 or abs(gun.energy-20)>1e-8:return
  tr=step.GetTrack();pre=step.GetPreStepPoint();post=step.GetPostStepPoint();proc=post.GetProcessDefinedStep();sec=step.GetSecondaryInCurrentStep()
  row={'track':tr.GetTrackID(),'parent':tr.GetParentID(),'particle':tr.GetParticleDefinition().GetParticleName(),'step':tr.GetCurrentStepNumber(),'process':proc.GetProcessName() if proc else None,'preK':pre.GetKineticEnergy()/g.MeV,'postK':post.GetKineticEnergy()/g.MeV,'deposit':step.GetTotalEnergyDeposit()/g.MeV,'status':str(tr.GetTrackStatus()),'preRminusR_m':(pre.GetPosition().mag()-R)/g.m,'postRminusR_m':(post.GetPosition().mag()-R)/g.m,'secondaries':[{'particle':x.GetParticleDefinition().GetParticleName(),'K':x.GetKineticEnergy()/g.MeV} for x in sec]}
  print('AUDIT_STEP',json.dumps(row),flush=True)

rm=g.G4RunManagerFactory.CreateRunManager(g.G4RunManagerType.Serial);det=Detector();physics=g.G4VModularPhysicsList();em_constructor=g.G4EmStandardPhysics_option4();physics.RegisterPhysics(em_constructor);physics.SetVerboseLevel(0);rm.SetUserInitialization(det);rm.SetUserInitialization(physics);em=g.G4EmParameters.Instance();em.SetVerbose(0);em.SetFluo(True);em.SetAuger(True);em.SetDeexcitationIgnoreCut(True);em.SetLowestElectronEnergy(100*g.eV);g.G4ProductionCutsTable.GetProductionCutsTable().SetEnergyRange(100*g.eV,100*g.GeV)
event=Event();gun=Generator();tracking=Tracking(event);rm.SetUserAction(gun);rm.SetUserAction(event);rm.SetUserAction(tracking);debugstep=DiagnosticStep();rm.SetUserAction(debugstep);rm.SetVerboseLevel(0);rm.Initialize();ui=g.G4UImanager.GetUIpointer();print('CUT_COMMAND',ui.ApplyCommand(f'/run/setCut {a.cut_um} um'),flush=True);g.G4Random.setTheSeed(a.seed)
records=[]
for E in [float(x) for x in a.energies.split(',')]:
 gun.energy=E;event.rows=[];t=time.monotonic();rm.BeamOn(a.histories);rows=event.rows;leaks=[x['escaped_calorimetric_MeV']/E for x in rows];deps=[x['deposited_MeV']/E for x in rows];balance=[x['deposited_MeV']+x['escaped_calorimetric_MeV']+x['vacuum_deposit_MeV']+x['residual_inside_kinetic_MeV']-E for x in rows];N=len(rows);se=statistics.stdev(leaks)/math.sqrt(N) if N>1 else 0.
 rec={'incident_energy_MeV':E,'histories':N,'deposited_fraction':statistics.mean(deps),'leakage_fraction':statistics.mean(leaks),'sampling_standard_error':se,'sampling_80pct_normal_interval':[max(0,statistics.mean(leaks)-1.281551566*se),min(1,statistics.mean(leaks)+1.281551566*se)],'max_abs_event_energy_residual_MeV':max(abs(x) for x in balance),'mean_energy_residual_MeV':statistics.mean(balance),'max_steps':max(x['steps'] for x in rows),'mean_steps':statistics.mean(x['steps'] for x in rows),'mean_tracks':statistics.mean(x['tracks'] for x in rows),'residual_inside_kinetic_MeV':sum(x['residual_inside_kinetic_MeV'] for x in rows),'vacuum_deposit_MeV':sum(x['vacuum_deposit_MeV'] for x in rows),'escape_by_particle_MeV':{p:sum(x['escape_by_particle_MeV'].get(p,0) for x in rows) for p in ['gamma','e-','e+']},'seconds':time.monotonic()-t,'events':rows};records.append(rec);print('RESULT',json.dumps({k:v for k,v in rec.items() if k!='events'}),flush=True)
 
 if not (N==a.histories and rec['max_abs_event_energy_residual_MeV']<1e-7*E):
  fail=a.out.with_suffix('.energy-failure.json');fail.write_text(json.dumps({'parameters':vars(a)|{'out':str(a.out),'data':str(a.data)},'records':records},allow_nan=False));print('ENERGY_FAILURE_PRESERVED',str(fail),flush=True);os._exit(2)
 assert rec['residual_inside_kinetic_MeV']<1e-7*N*E

print('DUMP_COMMAND',ui.ApplyCommand('/run/dumpCouples'),flush=True);thresholds={p:g.G4ProductionCutsTable.GetProductionCutsTable().ConvertRangeToEnergy(g.G4ParticleTable.GetParticleTable().FindParticle(p),det.mat,a.cut_um*g.um)/g.MeV for p in ['gamma','e-','e+']};assert all(v>0 for v in thresholds.values());tol=g.G4GeometryTolerance.GetInstance();out={'scope':__doc__,'geant4_version':g.G4Version,'geant4_version_number':g.G4VERSION_NUMBER,'binding_version':'0.1.3','parameters':vars(a)|{'out':str(a.out),'data':str(a.data)},'thickness_m':thick/g.m,'graphite_density_kg_m3':2000,'surface_tolerance_m':tol.GetSurfaceTolerance()/g.m,'radial_tolerance_m':tol.GetRadialTolerance()/g.m,'radius_ulp_m':math.ulp(R)/g.m,'requested_cut_m':a.cut_um*1e-6,'cut_configuration':'Geant4 UI /run/setCut after Initialize; actual built material-couple table printed by /run/dumpCouples and verified by campaign script','production_thresholds_MeV':thresholds,'lowest_electron_energy_MeV':em.LowestElectronEnergy()/g.MeV,'source_surface_zero_offset_m':-1e-6,'physics':'G4EmStandardPhysics_option4; fluorescence+Auger+deexcitationIgnoreCut;100eV lowest electron; no photonuclear/hadronic processes','scoring':'Native G4PSEnergyDeposit for hot graphite and vacuum; Python tracking callbacks only tally final outer escapes and residuals. World boundary10m outside receiver.','secondary_creator_counts':dict(event.creators),'particle_track_counts':dict(event.particles),'records':records,'elapsed_seconds':time.monotonic()-start,'process_exit':'All requested histories, energy checks and output persistence completed before os._exit to avoid observed binding destruction-order fault. No physics steps or remaining tracks truncated.','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};a.out.parent.mkdir(parents=True,exist_ok=True)
with a.out.open('w') as f:f.write(json.dumps(out,indent=2,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
print('COMPLETE',str(a.out),flush=True);sys.stdout.flush();sys.stderr.flush();
if not a.normal_exit:os._exit(0)
