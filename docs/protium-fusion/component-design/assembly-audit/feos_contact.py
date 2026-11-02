"""Read-only FEOS helper for the four declared source compositions.
Own material database and runtime directory; no source tables are modified.
"""
from pathlib import Path
import ctypes as c
import os, hashlib
AR=7.5657e-16
RECIPES={'A':dict(A=14,Z=7,H=3,rho=850.,T=150.),'B1':dict(A=15,Z=7,H=2,rho=1062.5,T=250.),'B2':dict(A=12,Z=6,H=2,rho=944.2,T=300.),'C13':dict(A=13,Z=6,H=4,rho=466.4375,T=102.4)}
class FEOS:
 def __init__(self,cache=Path('/home/tkwa/src/protium-fusion/output/assembly-audit/feos'),library_root=Path('/home/tkwa/src/protium-fusion/output/launcher-feos-final/FEOS')):
  cache.mkdir(parents=True,exist_ok=True);db=''
  for j,(stage,z) in enumerate(RECIPES.items(),1):
   n=9800+j;z['entity']=j
   db+=f'\nMaterial-{n}:\n[{n}]_SESAME-Number = 1000\n[{n}]_Treference = {z["T"]/11604.51812}\n[{n}]_Rhoreference = {z["rho"]/1000}\n[{n}]_Bulk-Modulus = 2.0e10\n[{n}]_Number-of-Elements = 2\n[{n}]_A[1] = {z["A"]}\n[{n}]_Z[1] = {z["Z"]}\n[{n}]_X[1] = 1\n[{n}]_A[2] = 1\n[{n}]_Z[2] = 1\n[{n}]_X[2] = {z["H"]}\n'
  self.db=db;(cache/'FEOS_Material-DB.dat').write_text(db)
  table=cache/'FEOS_TF-Table_1197.dat'
  if not table.exists():table.symlink_to(library_root/'EOS-Data/FEOS_TF-Table_1197.dat')
  self.library=library_root/'Code/libfeos.so';self.lib=c.CDLL(str(self.library));self.D=c.c_double;self.I=c.c_int
  DP=c.POINTER(self.D);IP=c.POINTER(self.I);l=self.lib
  l.FEOS_Initialize.argtypes=[self.I,self.I];l.FEOS_Init_Mat.argtypes=[self.I,self.I,self.I,self.I,DP,self.I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[self.I,self.I,self.I,self.D,self.D,IP,DP,DP,DP,DP,DP,DP]
  old=Path.cwd();os.chdir(cache)
  try:
   l.FEOS_Initialize(4,0)
   for j in range(1,5):
    ne=self.I();ni=self.I();tt=self.D();l.FEOS_Init_Mat(j,9800+j,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
  finally:os.chdir(old)
 def state(self,stage,rho,T):
  b=self.I();p=self.D();e=self.D();s=self.D();f=self.D();q=(self.D*2)();qt=self.D()
  self.lib.FEOS_Get_EOS(RECIPES[stage]['entity'],0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt))
  return {'rho_kg_m3':float(rho),'T_K':float(T),'p_Pa':p.value*.1+AR*T**4/3,'u_J_kg':e.value*1e-4+AR*T**4/rho,'s_J_kg_K':s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),'mean_charges':list(q),'status':b.value}
