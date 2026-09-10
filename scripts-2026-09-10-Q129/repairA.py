import sympy as sp
from ribbon import *
t,z = sp.Symbol('t'), sp.Symbol('z')

def build(e,N):
    P=all_parts_upto(N); idx={p:i for i,p in enumerate(P)}; n=len(P)
    R=sp.zeros(n,n)
    for mu in P:
        for lam,ht in add_ribbons(mu,e,N+3):
            if sum(lam)<=N: R[idx[lam],idx[mu]] += t**ht
    return P,idx,R

for e in (2,3,4):
    N = 8 if e<4 else 9
    P,idx,R=build(e,N); n=len(P)
    I=sp.eye(n)
    # Zc = (I - z R)^{-1}, computed as a finite Neumann series (R nilpotent on truncation)
    Zc=sp.eye(n); T=sp.eye(n); k=0
    while True:
        T=(z*R)*T
        T=sp.expand(T)
        if T==sp.zeros(n,n): break
        Zc=Zc+T; k+=1
    Zc=sp.expand(Zc)
    # 1. inverse is two-term?
    chk1 = sp.expand(Zc*(I-z*R)) == I
    # 2. conjugation reciprocity:  Zc(z,1/t) == E D^{-1} Zc(z,t) D E,  E=omega perm, D=t^{(e-1)w_e}
    E=sp.zeros(n,n)
    for p in P: E[idx[conj(p)],idx[p]]=1
    D=sp.diag(*[t**((e-1)*ecore_weight(p,e)[1]) for p in P])
    lhs=sp.expand(Zc.subs(t,1/t).applyfunc(sp.cancel))
    rhs=sp.expand(E*D.inv()*Zc*D*E).applyfunc(sp.cancel)
    chk2 = sp.simplify(lhs-rhs)==sp.zeros(n,n)
    # 3. the ADIN-BAUER INVERSION form: Zc^{-1} =?= eps D^{-1} Zc(1/t) D  for eps in {I,E,-I,-E}
    inv=sp.expand(I-z*R)
    res={}
    for nm,eps in (("I",I),("omega",E),("-I",-I),("-omega",-E)):
        cand=sp.expand(eps*D.inv()*Zc.subs(t,1/t)*D).applyfunc(sp.cancel)
        res[nm]= sp.simplify(cand-inv)==sp.zeros(n,n)
    print(f"e={e} N={N} dim={n} chain-length max={k}")
    print(f"   [A1] Zc * (I - zR) == I                       : {chk1}")
    print(f"   [A2] Zc(z,1/t) == omega D^-1 Zc(z,t) D omega  : {chk2}")
    print(f"   [A3] Zc^-1 == eps D^-1 Zc(1/t) D              : {res}")
