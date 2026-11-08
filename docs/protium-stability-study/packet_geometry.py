"""Translate angular direction errors into shaped-packet geometry, not fluid stability.

Each packet has one bulk direction. The circular-cap reference has solid angle
4*pi/N and c=cos(alpha)=1-2/N. Its mean sin^2(angle) is 2/N-4/(3*N^2).
A Fibonacci Voronoi tiling gives a realizable, approximately equal-area check.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def directions(n):
    z=1-2*(np.arange(n)+.5)/n
    phi=np.arange(n)*np.pi*(3-np.sqrt(5))
    s=np.sqrt(1-z*z)
    return np.column_stack([s*np.cos(phi),s*np.sin(phi),z])

def antipodal_directions(n):
    assert n%2==0
    z=(np.arange(n//2)+.5)/(n//2)
    phi=np.arange(n//2)*np.pi*(3-np.sqrt(5))
    s=np.sqrt(1-z*z)
    v=np.column_stack([s*np.cos(phi),s*np.sin(phi),z])
    return np.stack([v,-v],axis=1).reshape(-1,3)

def meridional(divisions):
    # High-order Gauss quadrature split at every sector edge.
    x,w=np.polynomial.legendre.leggauss(32)
    second=0.
    for j in range(divisions+1):
        c=j*np.pi/divisions
        lo=max(0,c-np.pi/(2*divisions));hi=min(np.pi,c+np.pi/(2*divisions))
        t=(hi+lo)/2+(hi-lo)/2*x
        second+=np.sum(w*np.sin(t)*np.sin(t-c)**2)*(hi-lo)/4
    n=(1+np.sqrt(1-4*second/3))/second
    return dict(divisions=divisions,direction_sine_rms=float(np.sqrt(second)),
                circular_cap_equivalent_N=float(n))

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    sample=directions(131072)
    rows=[]
    for n in [16,32,64,128,256,512,1024,2048,4096,8192]:
        sites=directions(n)
        distance,index=cKDTree(sites).query(sample,workers=1)
        cosine=1-distance**2/2
        angle2=1-cosine*cosine
        areas=np.bincount(index,minlength=n)/len(sample)
        rows.append(dict(N=n,sine_rms=float(np.sqrt(angle2.mean())),
            cap_sine_rms=float(np.sqrt(2/n-4/(3*n*n))),
            area_rms_fraction=float(np.std(areas)*n),
            max_angle_deg=float(np.degrees(np.arccos(cosine.min())))))
    result=dict(meaning='Geometry only: no fusion-output claim',sample_points=len(sample),
        tilings=rows,meridional=[meridional(k) for k in [4,8,16,32,64]])
    result['antipodal_3d_layouts']=[]
    for n in [32,64,96,128,256]:
        distance,index=cKDTree(antipodal_directions(n)).query(sample,workers=1)
        cosine=1-distance**2/2;areas=np.bincount(index,minlength=n)/len(sample)
        result['antipodal_3d_layouts'].append(dict(N=n,sine_rms=float(np.sqrt(np.mean(1-cosine*cosine))),
            area_rms_fraction=float(np.std(areas)*n)))
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['meridional'],indent=2))
if __name__=='__main__':main()
