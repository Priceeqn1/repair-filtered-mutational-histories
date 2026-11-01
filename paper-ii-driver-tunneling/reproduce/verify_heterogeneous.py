"""Independent verification for Paper II heterogeneous post-hit extension."""
import math
from scipy.optimize import brentq
import sympy as sp
from verify_finite_M import result, p

# Exact parameter choices
p['one_c'] = (8/25, 1.0)
p['one_l'] = (1/10, 2.0)
p['dc'] = (1/2, 1.0, 1.0)
p['dl'] = (1/2, 5.0, 1/10)

bC,dC=p['one_c']; bL,dL=p['one_l']
BC,dc,tc=p['dc']; BL,dl,tl=p['dl']
fc=bC/(dC-bC); fl=bL/(dL-bL)
aL=dL-bL; T=aL*fc-bL
qC=2*tc/((dc+tc-BC)+math.sqrt((dc+tc-BC)**2+4*BC*tc))
def FL(q):return BL*q*q+(dl+tl-BL)*q-tl

def qL(nu):
    a=dl+tl-BL+nu
    return 2*(tl+nu*qC)/(a+math.sqrt(a*a+4*BL*(tl+nu*qC)))
def phi(nu):
    x=(1+math.sqrt(1+16*nu*(2+4*nu)))/(2*(2+4*nu))
    return 1/(1+x)
def D_inf(nu):
    q=qL(nu); al=(bL*q+nu*fc*qC)/(aL+nu);u=phi(nu)
    return 2*((1-u)*fc*qC+u*al)/((1-u)*qC+u*q)

def W(q):return fc*FL(q)-T*q
for nu in [0,.1,.5,2,3,20,100]:
    q=qL(nu);u=phi(nu);Q=(1-u)*qC+u*q
    lhs=D_inf(nu)-2*fc
    rhs=2*u*W(q)/((aL+nu)*Q)
    assert abs(lhs-rhs)<1e-11,(nu,lhs,rhs)
print('Exact ratio identity: verified by floating-point substitution.')

q=sp.symbols('q');F=sp.Rational(1,2)*q*q+sp.Rational(23,5)*q-sp.Rational(1,10)
f=sp.Rational(8,17);t=sp.Rational(27,34)
Qc=(sp.sqrt(17)-3)/2
assert sp.simplify(f*F.subs(q,Qc)-t*Qc - (sp.Rational(36,85)+sp.Rational(113,170)*Qc))==0
print('Symbolic overshoot certificate:', sp.simplify(f*F.subs(q,Qc)-t*Qc))
print('Unbounded: D(0)=%.9f D(2)=%.9f D(infinity)=%.9f' % (D_inf(0),D_inf(2),2*fc))
assert D_inf(0)<1 and D_inf(2)>1 and 2*fc<1
for a,b in [(0,2),(2,100)]:
    print('Unbounded crossing in (%g,%g): %.9f'%(a,b,brentq(lambda nu:D_inf(nu)-1,a,b)))

for M in [5,6,8,10]:
    v0=result(M,0)['D'];v2=result(M,2)['D'];r=result(M,0)
    vinfty=2*r['AC']/r['qC']
    print('M=%d: D(0)=%.9f D(2)=%.9f D(infty)=%.9f'%(M,v0,v2,vinfty))
    if M>=6:assert v0<1 and v2>1 and vinfty<1
print('All checks passed.')
