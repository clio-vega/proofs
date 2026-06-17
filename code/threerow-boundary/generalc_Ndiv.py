"""
generalc_Ndiv.py -- prove session 2026-06-17

The crux. For deep boundary indices the deficit factors (a-c+2)=2(P+b/2+1) [b even]
and (a-b+1)=2(P+gamma) [c odd, gamma=(c+1)/2] fall OUTSIDE the consecutive run Pi_i,
so they cannot be absorbed by Lemma P.  Hypothesis: N_i then COMPENSATES, i.e. as
INTEGERS, v2(N_i) >= v2 of the un-absorbable deficit factor.

We test, for each (c,i,parity of b), the integer divisibility relations:
   does (P + b/2 + 1) | N_i  (b even)?   does (P + gamma) | N_i (c odd)?
and more refined:  v2(N_i) >= v2(P+b/2+1) ? v2(N_i)>= v2(P+gamma)?
Also: relate content g[c][k] to k via 2*g - (deficit-budget).
"""
import sympy as sp
from math import comb, factorial
from mn import Mj

a, b = sp.symbols('a b')

def collect(c, i, amax=64, bmax=28):
    pts = []
    for bb in range(c, bmax):
        for aa in range(bb, amax):
            if (aa + bb + c) % 2 != 0: continue
            m = (aa + bb + c) // 2
            j = bb + i
            if j > m: continue
            pts.append((aa, bb, Mj((aa, bb, c), j, m)))
    return pts

def fit_poly(pts, maxdeg):
    monos = [(p, q) for p in range(maxdeg+1) for q in range(maxdeg+1) if p+q<=maxdeg]
    if len(pts) < len(monos)+3: return None
    A = sp.Matrix([[aa**p*bb**q for (p,q) in monos] for (aa,bb,_) in pts])
    y = sp.Matrix([M for (_,_,M) in pts])
    sol = A.solve_least_squares(y)
    if any(A*sol-y): return None
    return sp.factor(sp.expand(sum(sol[k]*a**p*b**q for k,(p,q) in enumerate(monos))))

def extract_N(c, i):
    M = fit_poly(collect(c,i), maxdeg=2*c+1)
    if M is None: return None
    bprod = (b-c+1)*sp.prod([b+t for t in range(2,i+1)])
    if i==1: bprod = (b-c+1)*(a-b+1)
    ratio = sp.cancel(M/bprod)
    num, den = sp.fraction(sp.together(ratio))
    num = sp.expand(num)
    poly = sp.Poly(num, a, b)
    g=0
    for cc in poly.coeffs(): g = sp.igcd(g, int(cc))
    return sp.expand(num/g)

def v2int(n):
    n=abs(int(n))
    if n==0: return 999
    k=0
    while n%2==0: n//=2;k+=1
    return k

print("Testing deficit-compensation: does N_i carry the un-absorbable deficit's v2?")
print("ell=ceil((b+c+3)/2); Pi_i lower end = P+k+1, k=c-i.")
print("a-c+2 deficit factor = P+b/2+1 (b even); a-b+1 deficit factor = P+(c+1)/2 (c odd).")
print()

for c in [3,4,5,6]:
    print("="*70)
    print("c=%d"%c)
    for i in range(1,c+1):
        k=c-i
        N = extract_N(c,i)
        if N is None:
            print("  i=%d: no fit"%i); continue
        gamma = (c+1)//2 if c%2==1 else None
        # sweep, classify
        worst_abs_acm2 = 999  # min over (b even) of v2(N) - v2(P+b/2+1)
        worst_acm2_inrange = None
        worst_abp1 = 999
        abp1_inrange_count=0; abp1_outrange_count=0
        acm2_inrange_count=0; acm2_outrange_count=0
        for P in range(0,30):
            for bb in range(c, 60):
                aa=2*P+bb+c
                m=P+bb+c
                if aa<bb: continue
                w=bb+c
                ell = -(-(w+3)//2)  # ceil
                Llo=P+k+1; Lhi=P+ell-1
                Nv = int(N.subs({a:aa,b:bb}))
                vN = v2int(Nv)
                # a-c+2 deficit (b even)
                if bb%2==0:
                    fac = P + bb//2 + 1
                    inrange = (Llo <= fac <= Lhi)
                    if inrange: acm2_inrange_count+=1
                    else:
                        acm2_outrange_count+=1
                        diff = vN - v2int(fac)
                        if diff<worst_abs_acm2: worst_abs_acm2=diff
                # a-b+1 deficit (c odd), present only k<=c-2
                if c%2==1 and k<=c-2:
                    fac = P + gamma
                    inrange = (Llo <= fac <= Lhi)
                    if inrange: abp1_inrange_count+=1
                    else:
                        abp1_outrange_count+=1
                        diff = vN - v2int(fac)
                        if diff<worst_abp1: worst_abp1=diff
        msg="  i=%d (k=%d):"%(i,k)
        if acm2_outrange_count:
            msg+=" [a-c+2 OUT %d times: min(v2N - v2(P+b/2+1))=%d]"%(acm2_outrange_count, worst_abs_acm2)
        else:
            msg+=" [a-c+2 in-range always (%d)]"%acm2_inrange_count
        if c%2==1 and k<=c-2:
            if abp1_outrange_count:
                msg+=" [a-b+1 OUT %d: min(v2N - v2(P+gamma))=%d]"%(abp1_outrange_count, worst_abp1)
            else:
                msg+=" [a-b+1 in-range always (%d)]"%abp1_inrange_count
        print(msg)
