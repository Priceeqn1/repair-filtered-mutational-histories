import sympy as s
from itertools import product
from time import time
n=s.symbols('n')
from sympy.polys.matrices import DomainMatrix
T=time()
births=[s.Rational(306,755),s.Rational(1,5),s.Rational(12,5),s.Rational(2)]
deaths=[s.Rational(6,5),s.Rational(3,2),s.Rational(8,5),s.Rational(2)]
commits=[0,0,s.Rational(1),s.Rational(1,1000000)]
def compositions(n,k):
 if k==1:
  yield (n,);return
 for j in range(n+1):
  for tail in compositions(n-j,k-1): yield (j,)+tail
M=3
pre=[z for m in range(1,M) for z in compositions(m,2)]
post=[z for m in range(1,M) for z in compositions(m,4)]
def mat(states):
 idx={z:i for i,z in enumerate(states)}
 G=s.zeros(len(states),len(states)); t=s.zeros(len(states),1)
 for row,z in enumerate(states):
  rate=0
  for i,k in enumerate(z):
   if not k: continue
   for intensity,shift in [(k*births[i],1),(k*deaths[i],-1)]:
    rate+=intensity; nz=list(z);nz[i]+=shift;nz=tuple(nz)
    if 0<sum(nz)<M:G[row,idx[nz]]+=intensity
   if i in (1,3):
    rate+=k*n;nz=list(z);nz[i]-=1;nz[i-1]+=1;G[row,idx[tuple(nz)]]+=k*n
   if i>=2: rate+=k*commits[i];t[row]+=k*commits[i]
  G[row,row]-=rate
 return G,t
Q,_=mat(pre); K,t=mat(post)
B=s.zeros(len(pre),len(post)); pi={z:i for i,z in enumerate(post)}
for row,z in enumerate(pre):
 if sum(z)>=M-1:continue
 for i,k in enumerate(z):
  if k:
   nz=z+(0,0)
   nz=list(nz);nz[2+i]+=1
   B[row,pi[tuple(nz)]]+=k*births[i]
print('Exact M3 model sizes: pre=%d, post=%d'%(len(pre),len(post)),flush=True)
Qq=s.Matrix(-K).inv()*t
# q obtained by exact symbolic matrix inversion
Aa=(-Q).inv()*B*Qq
# A obtained by exact symbolic matrix inversion
qa=Qq[pi[(0,0,1,0)]];qb=Qq[pi[(0,0,0,1)]]
aa=Aa[pre.index((1,0))];ab=Aa[pre.index((0,1))]
# The pristine direct coefficient is independent of repair; damaged coefficients are rational.
x=s.symbols('x');nv=x*(2*x-1)/(4*(1-x*x))
D=s.factor(2*(x*aa+ab)/(x*qa+qb));# Full population-weighted ratio D assembled
D_x=s.factor(s.cancel(D.subs(n,nv)));# Uniformize repair via x=(1+t)/(2+t).
num,den=s.fraction(D_x);print('num degree',s.degree(num,x),'den degree',s.degree(den,x));assert s.gcd(num,den)==1
par=s.Poly(num-den,x);stat=s.Poly(s.diff(num,x)*den-num*s.diff(den,x),x)

def vari(p):
 coefs=p.all_coeffs();sg=[s.sign(c) for c in coefs if c!=0]
 return sum(sg[j]!=sg[j+1] for j in range(len(sg)-1))
for name,p in [('parity',par),('stationary',stat)]:
 tr=p.transform(s.Poly(x+1,x),s.Poly(x+2,x));print(name,'Descartes variations',vari(tr),'degree',tr.degree(),'at ends',s.sign(p.eval(s.Rational(1,2))),s.sign(p.eval(1)),flush=True)
for name,p in [('parity',par),('stationary',stat)]:
 tr=p.transform(s.Poly(x+1,x),s.Poly(x+2,x)); print(name,'coefficient signs hi-to-low:', ''.join('+' if c>0 else '-' if c<0 else '0' for c in tr.all_coeffs()))
 print(name,'sample signs:',[(str(u),s.sign(p.eval(u))) for u in [s.Rational(1,2),s.Rational(91,100),s.Rational(92,100),s.Rational(1)]])
 print(name,'degree',p.degree(), 'gcd(n,d)', s.polys.polytools.degree(s.gcd(num,den),x))
assert vari(par.transform(s.Poly(x+1,x),s.Poly(x+2,x)))==0
assert vari(stat.transform(s.Poly(x+1,x),s.Poly(x+2,x)))==1
assert s.sign(stat.eval(s.Rational(1,2)))==1 and s.sign(stat.eval(1))==-1
# Thus exactly one stationary point, necessarily a strict maximum; no parity root.
print('CERTIFIED: exactly one interior maximum and no D3=1 crossing.')
print('Total calculation time %.3f seconds'%(time()-T))
