import sys, os, random
import sympy as sp
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','code-1017c2-coboundary'))
sys.path.insert(0,os.path.dirname(__file__))
import q484_engine as E
import ls as LSM
from fractions import Fraction as Fr
t=sp.Symbol('t')
import census_S as C

for prec in (20,40,80,160):
    LSM.K=prec
    import importlib
    random.seed(1)
    N,m,eps=3,1,1
    w=[]
    for j in range(1,N+1):
        w.append(sp.Integer(eps)+t if j==m else sp.Rational(random.randint(3,40),random.randint(1,7)))
    q0=sp.Rational(3,7)
    P=E.subsets(N); idx={I:k for k,I in enumerate(P)}; n=len(P)
    Phi=[[None]*n for _ in range(n)]
    vals={}
    for M in P:
        for I in P:
            a=E.Atilde(N,w,M,I,X='C').subs({E.q:q0})
            l=C.to_LS(a,prec)
            Phi[idx[M]][idx[I]]=l
            vals[(M,I)]=None if l.is_zero() else l.o
    _,ordet=LSM.solve(Phi,[LSM.LS.zero()]*n,n)
    diagsum=sum(vals[(I,I)] or 0 for I in P)
    minent=min(v for v in vals.values() if v is not None)
    print(f"prec={prec}: ordet={ordet}  sum_diag={diagsum}  min entry val={minent}")
    if prec==20:
        for M in P:
            print("   ",M,[vals[(M,I)] for I in P])
