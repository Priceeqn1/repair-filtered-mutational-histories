"""Independent exact (SymPy) verification of opportunity balance/counterexample."""
import itertools
import sympy as s
nu=s.symbols('nu', nonnegative=True)
STATES=[z for n in (1,2) for z in itertools.product(range(n+1),repeat=2) if sum(z)==n]
idx={z:i for i,z in enumerate(STATES)}
B=(s.Integer(10),s.Integer(4));D=(s.Integer(1),s.Integer(1));M=3
n=len(STATES);Q=s.zeros(n);h=s.zeros(n,1);death=s.zeros(n,1);est=s.zeros(n,1)
for i,z in enumerate(STATES):
    zz=list(z);size=sum(z)
    for k in range(2):
        if not z[k]:continue
        rate=B[k]*z[k]; Q[i,i]-=rate
        if size+1==M:est[i]+=rate
        else:
            w=zz.copy();w[k]+=1;Q[i,idx[tuple(w)]]+=rate
            h[i]+=rate
        rate=D[k]*z[k];Q[i,i]-=rate; death[i]+=rate
        if size-1>0:
            w=zz.copy();w[k]-=1;Q[i,idx[tuple(w)]]+=rate
    rate=nu*z[1];Q[i,i]-=rate
    if rate!=0:
        w=zz.copy();w[1]-=1;w[0]+=1;Q[i,idx[tuple(w)]]+=rate
G=(-Q).inv()
i=idx[(0,1)]
vals={name:s.factor((G*v)[i]) for name,v in [('opp',h),('deaths',death),('est',est)]}
for name,val in vals.items():
    print(name,'value0=',val.subs(nu,0),'derivative0=',s.factor(s.diff(val,nu).subs(nu,0)),'limit=',s.limit(val,nu,s.oo))
print('balance residual:',s.simplify(vals['opp']-vals['deaths']-2*vals['est']+1))
assert s.simplify(vals['opp']-vals['deaths']-2*vals['est']+1)==0
assert s.diff(vals['opp'],nu).subs(nu,0)<0
print('PASS')
