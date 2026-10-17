"""Targeted: lam=mu=() at N=3, type C, q=3/7.  WHICH divisors carry the poles for each chi?"""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1017c2-coboundary')
import sympy as sp
from collections import Counter
import q484_engine as E
qs=E.q; qv=sp.Rational(3,7)
CHI={'id':lambda v:sp.Integer(1),
     'psi':lambda v:(1-qv*v**2)/(1-qv),
     'xpsi':lambda v:v*(1-qv*v**2)/(1-qv),
     'x':lambda v:v}
N=3; w=E.wvars(N); S=sorted(E.subsets(N),key=lambda M:(len(M),M))
A0={(M,I):E.Atilde(N,w,M,I,'C',simplify=False).subs(qs,qv) for M in S for I in S}
for lab,chi in CHI.items():
    A={}
    for M in S:
        P=sp.Integer(1)
        for m in M: P*=chi(w[m-1])
        for I in S: A[(M,I)]=sp.cancel(sp.together(A0[(M,I)]*P))
    lam=mu=()
    c={}
    print(f'--- chi={lab} ---',flush=True)
    for nu in S:
        v=A[(nu,lam)]*A[(nu,mu)]; acc=sp.Integer(0)
        for kap in S:
            if kap==nu: break
            if c.get(kap,0)!=0: acc+=c[kap]*A[(nu,kap)]
        c[nu]=sp.cancel(sp.together((v-acc)/A[(nu,nu)]))
        if c[nu]==0: continue
        num,den=sp.fraction(c[nu])
        fl=[(f,m) for f,m in sp.factor_list(sp.expand(den))[1] if len(sp.Poly(f,*w).monoms())>1]
        mono=[(f,m) for f,m in sp.factor_list(sp.expand(den))[1] if len(sp.Poly(f,*w).monoms())==1]
        print(f'  nu={nu}  nonmono-den: {[(sp.sstr(f),m) for f,m in fl]}   mono-den: {[(sp.sstr(f),m) for f,m in mono]}',flush=True)
