import sympy as sp
from sympy import symbols, expand, Poly, I, im
import math
m = symbols('m', real=True)
def ImP_coeff_b(b):
    s = 1 + I
    def rl(l):
        tot = 0
        for c in range(0, l//2 + 1):
            k = l - 2*c
            ff = 1
            for j in range(0, k + c):
                ff *= (m - j)
            tot += ff / (sp.factorial(k) * sp.factorial(c)) * s**k
        return sp.expand(tot)
    return expand(im(rl(b) - rl(b-1)))
def v2(n):
    n=int(n)
    if n==0: return 10**9
    v=0
    while n%2==0: n//=2; v+=1
    return v
def C(n,k):
    if k<0 or k>n: return 0
    return math.comb(n,k)
def imv2(j):
    if j%4==0: return None
    return j//2
def term_v2(b,mm,j):
    g=imv2(j)
    if g is None: return None
    l=b-j
    beta = l//2 if l%2==0 else (b-1-j)//2
    cm=C(mm,j); bn=C(mm-j,beta)
    if cm==0 or bn==0: return None
    return v2(cm)+g+v2(bn)

# (A) aggregate identity for b=1 mod4, large m
for b in [5,9,13,17,21]:
    R=(b-1)//2; vR=v2(math.factorial(R))
    Ib=Poly(ImP_coeff_b(b),m)
    agg_ok=True; ind_ok=True; aggfail=[]; indfail=[]
    import random
    ms=list(range(b,b+200))+[random.randint(b,20000) for _ in range(400)]
    for mm in ms:
        Iv=int(Ib.eval(mm))
        prod=1
        for r in range(R+1): prod*=(mm-r)
        if v2(Iv)!=v2(prod)-vR: agg_ok=False; aggfail.append(mm)
        t1=term_v2(b,mm,1)
        for j in range(2,b+1):
            tj=term_v2(b,mm,j)
            if tj is not None and tj<=t1: ind_ok=False; indfail.append((mm,j,t1,tj)); break
    print(f"b={b}: aggregate v2 identity {'OK' if agg_ok else 'FAIL '+str(aggfail[:3])}; "
          f"individual-domination {'OK' if ind_ok else 'FAIL '+str(indfail[:3])}")
