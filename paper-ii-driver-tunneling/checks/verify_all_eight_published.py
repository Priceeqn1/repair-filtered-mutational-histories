"""Rational-residual certificate of all eight published M=9 and M=10 witnesses.
No complete root count is claimed. Numerical candidates can have arbitrary
errors: the exact rational residual bounds verify the printed inequalities.
"""
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

ONE_C=(F(306,755),F(6,5)); ONE_L=(F(1,5),F(3,2))
TWO_C=(F(12,5),F(8,5),F(1));TWO_L=(F(2),F(2),F(1,1000000))
bs=[ONE_C[0],ONE_L[0],TWO_C[0],TWO_L[0]]
ds=[ONE_C[1],ONE_L[1],TWO_C[1],TWO_L[1]]
cs=[F(0),F(0),TWO_C[2],TWO_L[2]]
THETA_MIN=F(1,1000000)
a=min(ONE_C[1]-ONE_C[0],ONE_L[1]-ONE_L[0])

def comps(n,k):
 if k==1:
  yield (n,);return
 for j in range(n+1):
  for tail in comps(n-j,k-1):yield (j,)+tail

def prepare(M,nu):
 pre=[z for n in range(1,M) for z in comps(n,2)]
 active=[z for n in range(1,M) for z in comps(n,4) if z[2]+z[3]>0]
 Ip={z:i for i,z in enumerate(pre)}; Ik={z:i for i,z in enumerate(active)}
 def gen(states,ids):
  rows=[{} for _ in states];t=[F(0) for _ in states]
  for i,z in enumerate(states):
   out=F(0)
   for k,cnt in enumerate(z):
    if not cnt:continue
    for rate,delta in ((cnt*bs[k],1),(cnt*ds[k],-1)):
     out+=rate; nz=list(z);nz[k]+=delta;nz=tuple(nz)
     j=ids.get(nz)
     if j is not None:rows[i][j]=rows[i].get(j,F(0))+rate
    if k in (1,3):
     out+=cnt*nu;nz=list(z);nz[k]-=1;nz[k-1]+=1
     j=ids.get(tuple(nz))
     if j is not None:rows[i][j]=rows[i].get(j,F(0))+cnt*nu
    if k>=2:out+=cnt*cs[k];t[i]+=cnt*cs[k]
   rows[i][i]=rows[i].get(i,F(0))-out
  return rows,t
 Q,_=gen(pre,Ip)
 K,t=gen(active,Ik)
 B=[{} for _ in pre]
 for i,z in enumerate(pre):
  if sum(z)>=M-1:continue
  for k,cnt in enumerate(z):
   if not cnt:continue
   zz=list(z)+[0,0];zz[2+k]+=1
   j=Ik[tuple(zz)];B[i][j]=B[i].get(j,F(0))+cnt*bs[k]
 return pre,active,Ip,Ik,Q,K,B,t

def numeric_solve(rows,rhs):
 r=[];c=[];v=[]
 for i,row in enumerate(rows):
  for j,t in row.items():r.append(i);c.append(j);v.append(-float(t))
 M=coo_matrix((v,(r,c)),shape=(len(rows),len(rows))).tocsr()
 return spsolve(M,np.array([float(t) for t in rhs],dtype=float))

def exact_residual(rows,rhs,sol):
 maxabs=F(0)
 for i,row in enumerate(rows):
  v=rhs[i]+sum((coeff*sol[j] for j,coeff in row.items()),F(0))
  maxabs=max(maxabs,abs(v))
 return maxabs

def source_bounds(nu):
 if nu==0:return F(1,2),F(1,2)
 # The positive equilibrium branch solves (2+4nu)x^2-x-4nu=0.
 xf=(1+math.sqrt(1+16*float(nu)*(2+4*float(nu))))/(2*(2+4*float(nu)))
 eps=F(1,10**11)
 xl=F.from_float(xf)-eps;xh=F.from_float(xf)+eps
 def p(x):return (2+4*nu)*x*x-x-4*nu
 assert F(1,2)<xl<xh<1
 assert p(xl)<0<p(xh),('bad source enclosure',nu)
 return xl,xh

def certificate(M,nu):
 pre,post,Ip,Ik,Q,K,B,t=prepare(M,nu)
 qtilde=numeric_solve(K,t)
 qrat=[F.from_float(float(v)) for v in qtilde]
 rq=exact_residual(K,t,qrat)
 Eq=rq/THETA_MIN
 Bq=[sum((coeff*qrat[j] for j,coeff in row.items()),F(0)) for row in B]
 atilde=numeric_solve(Q,Bq)
 arat=[F.from_float(float(v)) for v in atilde]
 ra=exact_residual(Q,Bq,arat)
 # Inequalities: (-K)1 >= theta_min on states with D>0;
 # (-Q)(|z|/a) >= 1, using strictly subcritical one-driver types.
 Bnorm=max((sum(row.values(),F(0)) for row in B),default=F(0))
 Ea=F(M-1,1)/a*(ra+Bnorm*Eq)
 xl,xh=source_bounds(nu)
 ac,al=arat[Ip[(1,0)]],arat[Ip[(0,1)]]
 qc,ql=qrat[Ik[(0,0,1,0)]],qrat[Ik[(0,0,0,1)]]
 nlo=2*(xl*(ac-Ea)+(al-Ea));nhi=2*(xh*(ac+Ea)+(al+Ea))
 dlo=xl*(qc-Eq)+(ql-Eq);dhi=xh*(qc+Eq)+(ql+Eq)
 assert nlo>0 and dlo>0
 low=nlo/dhi; high=nhi/dlo
 sign='>1' if low>1 else '<1' if high<1 else 'INDETERMINATE'
 print('M=',M,'nu=',nu,'active states=',len(post),'pre states=',len(pre))
 print('   exact rational residual norms',float(rq),float(ra))
 print('   conservative certified errors q,A:',float(Eq),float(Ea))
 print('   certified ratio interval: [%.12f, %.12f] %s'%(float(low),float(high),sign),flush=True)
 assert sign!='INDETERMINATE'
 return low,high

if __name__=='__main__':
 published = {
  (9,F(0)):(F('0.9984692'),F('0.9984744')),
  (9,F(1,10)):(F('0.9645279'),F('0.9645337')),
  (9,F(6)):(F('1.0092026'),F('1.0092074')),
  (9,F(1000)):(F('0.9983638'),F('0.9989161')),
  (10,F(0)):(F('1.0079118'),F('1.0079253')),
  (10,F(1,10)):(F('0.9622777'),F('0.9622879')),
  (10,F(10)):(F('1.0166344'),F('1.0166487')),
  (10,F(100)):(F('1.0093700'),F('1.0094907')),
}
if __name__=='__main__':
  for M,rates in [(9,[F(0),F(1,10),F(6),F(1000)]),(10,[F(0),F(1,10),F(10),F(100)])]:
    for nu in rates:
      lo,hi=certificate(M,nu)
      pl,ph=published[(M,nu)]
      assert pl<lo<hi<ph,(M,nu,float(lo),float(hi))
      print('   MANUSCRIPT PUBLISHED ENCLOSURE STRICTLY VERIFIED')
