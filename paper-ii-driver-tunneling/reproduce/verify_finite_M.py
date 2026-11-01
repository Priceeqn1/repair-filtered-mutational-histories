import numpy as np
from itertools import combinations
from scipy.sparse import lil_matrix,csr_matrix
from scipy.sparse.linalg import spsolve
from math import comb
def xsource(nu):
    # Positive physical root of (2+4nu)x^2-x-4nu=0.
    from math import sqrt
    return (1+sqrt(1+16*nu*(2+4*nu)))/(2*(2+4*nu))
p=dict(one_c=(.3,1),one_l=(.1,2),dc=(.5,1,1),dl=(.5,10,.05))

def states(M,k):
 return [tuple(z) for n in range(1,M) for z in __import__('itertools').product(range(n+1),repeat=k) if sum(z)==n]
# efficient enumeration via stars and bars

def composition(n,k):
 if k==1:
  yield (n,);return
 for a in range(n+1):
  for rest in composition(n-a,k-1):yield (a,)+rest

def make(M,nu):
 one=[s for n in range(1,M) for s in composition(n,2)]
 post=[s for n in range(1,M) for s in composition(n,4)]
 ix={s:i for i,s in enumerate(one)};iy={s:i for i,s in enumerate(post)}
 Q=lil_matrix((len(one),len(one)));K=lil_matrix((len(post),len(post)));B=lil_matrix((len(one),len(post)))
 bv=[p['one_c'][0],p['one_l'][0],p['dc'][0],p['dl'][0]]
 dv=[p['one_c'][1],p['one_l'][1],p['dc'][1],p['dl'][1]]
 tv=[0,0,p['dc'][2],p['dl'][2]]
 def gen(ss,idx,G):
  ts=np.zeros(len(ss))
  for ii,z in enumerate(ss):
   n=sum(z);total=0
   for j,c in enumerate(z):
    if c==0:continue
    for rate,delta in [(c*bv[j],1),(c*dv[j],-1)]:
     total+=rate
     nz=list(z);nz[j]+=delta; nz=tuple(nz)
     if 0<sum(nz)<M:G[ii,idx[nz]]+=rate
    if j in (1,3):
     rate=c*nu;total+=rate
     nz=list(z);nz[j]-=1;nz[j-1]+=1;G[ii,idx[tuple(nz)]]+=rate
    if j>=2:
     rate=c*tv[j];total+=rate;ts[ii]+=rate
   G[ii,ii]-=total
  return ts
 gen(one,ix,Q); t=gen(post,iy,K)
 for ii,z in enumerate(one):
  n=sum(z)
  if n>=M-1:continue
  for j,c in enumerate(z):
   if c==0:continue
   dest=(z[0],z[1],int(j==0),int(j==1))
   B[ii,iy[dest]]+=c*bv[j]
 q=spsolve((-K).tocsr(),t);A=spsolve((-Q).tocsr(),B.tocsr().dot(q))
 return { 'qC':q[iy[(0,0,1,0)]], 'qL':q[iy[(0,0,0,1)]], 'AC':A[ix[(1,0)]], 'AL':A[ix[(0,1)]], 'sz':(len(one),len(post))}

def result(M,nu):
 r=make(M,nu)
 phi=1/(1+xsource(nu))
 r['D']=2*((1-phi)*r['AC']+phi*r['AL'])/((1-phi)*r['qC']+phi*r['qL'])
 return r

if __name__=='__main__':
 import time
 for M in [3,4,5,6,8,10,12]:
  t0=time.time()
  for nu in [0,.5,3,20,100]:
   r=result(M,nu)
   print(f'M={M:2}, nu={nu:5g}, D={r["D"]:.8f}, sizes={r["sz"]}')
  print('dt',time.time()-t0,flush=True)
