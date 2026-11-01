"""Audit common-conditioned-Green identity independently of finite differences."""
import numpy as np
from scipy.linalg import solve
from verify_generator import matrices,STATES,IDX,xof
N=len(STATES)
I1=np.array([i for i,z in enumerate(STATES) if z[2]==0])
I2=np.array([i for i,z in enumerate(STATES) if z[2]>0])

def decomp(nu,population=True):
    Qt,Bfull,tall=matrices(nu)
    Q=Qt[np.ix_(I1,I1)];K=Qt[np.ix_(I2,I2)];B=Bfull[np.ix_(I1,I2)];t=tall[I2]
    eps=1e-6
    P=(matrices(nu+eps)[0]-matrices(nu-eps)[0])/(2*eps)
    Qp=P[np.ix_(I1,I1)]; Kp=P[np.ix_(I2,I2)]
    H=np.linalg.inv(-K); G=np.linalg.inv(-Q)
    q=H@t; A=G@B@q
    alpha=np.zeros(len(I1));beta=np.zeros(len(I2))
    if population:
        x=xof(nu)
        xp=(xof(nu+eps)-xof(nu-eps))/(2*eps)
        alpha[list(I1).index(IDX[(1,0,0)])]=x
        alpha[list(I1).index(IDX[(0,1,0)])]=1
        alphap=np.zeros_like(alpha);alphap[list(I1).index(IDX[(1,0,0)])]=xp
    else:
        alpha[list(I1).index(IDX[(0,1,0)])]=1
        alphap=np.zeros_like(alpha)
    beta[list(I2).index(IDX[(0,0,1)])]=alpha.sum()
    betap=np.zeros_like(beta);betap[list(I2).index(IDX[(0,0,1)])]=alphap.sum()
    S=alpha@A;T=beta@q
    qc=H@(Kp@q);Ac=G@(Qp@A+B@qc)
    truth=(alphap@A+alpha@Ac)/S-(betap@q+beta@qc)/T
    vpre=(Qp@A)/A;vpost=(Kp@q)/q
    mseq=(alpha@G@B)*q/S;mdir=beta*q/T
    rho=(alpha@G)*A/S
    HH=np.diag(1/q)@H@np.diag(q)
    components=[alphap@A/S-betap@q/T,rho@vpre,(mseq-mdir)@HH@vpost]
    return truth, components, mseq.sum(),mdir.sum()

for ispop in (False,True):
    for nu in (.1,.8,2.):
        true,components,a,b=decomp(nu,ispop)
        assert abs(true-sum(components))<1e-9,(true,components)
        assert np.isclose(a,1) and np.isclose(b,1)
        print(f'population={ispop} nu={nu} derivative={true:.11f} '+
              'components='+str([round(z,11) for z in components]))
print('ALL GREEN RESPONSE TESTS PASS')
