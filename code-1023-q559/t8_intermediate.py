"""Verify the INTERMEDIATE steps the write-up asserts, not just the endpoint."""
import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from itertools import combinations
from psi import B as Bfun, mvec
from sampler import clean_w
import random
exec(open('t4_P4_blocks.py').read().split("rng = random.Random")[0].split('"""',2)[2])

rng=random.Random(8181)
ok1=bad1=0; ok2=bad2=0; ok3=bad3=0
for trial in range(80):
    r=rng.randint(1,3); n=rng.randint(r+1,r+2)
    wq=clean_w(n,rng,lo=2,hi=12,dlo=1,dhi=5)
    wsym=[Q(x) for x in wq[:-1]]+[W]; Mp=wq[:-1]
    # STEP 1: lim_{W->oo} B_M(T' u {n}) = B_{M'}(T')
    for Tp in combinations(range(n-1), r-1):
        T = tuple(sorted(Tp+(n-1,)))
        lim = sp.limit(Bsym(wsym,T), W, sp.oo)
        tgt = Q(Bfun(Mp, Tp))
        if sp.simplify(lim-tgt)==0: ok1+=1
        else:
            bad1+=1
            if bad1<4: print("  S1 FAIL", r,n,T,lim,tgt)
    # STEP 2: for T avoiding n,  lim W^r * B_M(T) = (-1)^r prod_{x in T} x * B_{M'}(T)
    for T in combinations(range(n-1), r):
        lim = sp.limit(W**r * Bsym(wsym,T), W, sp.oo)
        px = sp.Integer(1)
        for i in T: px *= Q(wq[i])
        tgt = (-1)**r * px * Q(Bfun(Mp,T))
        if sp.simplify(lim-tgt)==0: ok2+=1
        else:
            bad2+=1
            if bad2<4: print("  S2 FAIL", r,n,T,lim,tgt)
    # STEP 3: prod_{x in T} x * det(x^m) = det(x^{m+1})
    theta=[rng.randint(-4,4) for _ in range(r)]; m=mvec(theta)
    for T in combinations(range(n-1), r):
        x=[Q(wq[i]) for i in T]
        lhs = sp.Integer(1)
        for xi in x: lhs *= xi
        Mm = sp.Matrix(r,r, lambda i,j: x[i]**m[j])
        Mm1 = sp.Matrix(r,r, lambda i,j: x[i]**(m[j]+1))
        if sp.simplify(lhs*Mm.det() - Mm1.det())==0: ok3+=1
        else: bad3+=1
print("STEP 1  lim B_M(T' u {n}) = B_{M'}(T') :        OK=%d BAD=%d" % (ok1,bad1))
print("STEP 2  lim W^r B_M(T) = (-1)^r (prod x) B_M'(T):OK=%d BAD=%d" % (ok2,bad2))
print("STEP 3  (prod x) det(x^m) = det(x^{m+1}) :       OK=%d BAD=%d" % (ok3,bad3))
