import sympy as s
v=s.symbols('v',nonnegative=True)
x=(1+s.sqrt(64*v*v+32*v+1))/(8*v+4)
P=40000*v**4+2397000*v**3+15730820*v**2+16064982*v+1592283
Q=1000*v**3+8300*v**2+14430*v+6771
DL=s.Rational(2575,11752)*P/((5*v+263)*Q)
DC=s.Rational(2575,1469)
DP=(x*DC+DL)/(1+x)
print('x0',x.subs(v,0),'D0',s.factor(DP.subs(v,0)),float(DP.subs(v,0)))
print('endpoint',DC)
print('root',s.nsolve(DP-1,.2,prec=35))
print('DL0',float(DL.subs(v,0)),'DC',float(DC))
for z in [0,.1,.2,.5,1,2,4]:print(z,float(DP.subs(v,z)))
