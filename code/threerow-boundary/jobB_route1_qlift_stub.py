"""
jobB_route1_qlift_stub.py -- 2026-06-18 CODE, Job B route-1 (cheap stub).

Question (dream's highest-leverage route): can N_i^{(c)} be written as a ratio of
q-binomials whose cyclotomic factor Phi_2(q)=q+1 has multiplicity at q=-1 equal to
the 2-adic content?  That would be the Gatzweiler-Krattenthaler quotient (2502.06032).

This is a STUB.  It establishes the machinery (q-integer / q-binomial, cyclotomic
(q+1)-multiplicity counter, q=1 sanity) and reports HONESTLY what blocks a direct
q-lift, so the next (post-browse) cycle can drop in the GK determinantal q-form.
"""
import sympy as sp

q = sp.symbols('q')

def qint(n):
    """[n]_q = 1+q+...+q^{n-1} (n>=0)."""
    n = int(n)
    if n <= 0:
        return sp.Integer(0) if n == 0 else None
    return sp.expand(sum(q**i for i in range(n)))

def qfact(n):
    r = sp.Integer(1)
    for i in range(1, int(n) + 1):
        r *= qint(i)
    return sp.expand(r)

def qbinom(n, k):
    """Gaussian binomial [n,k]_q as a polynomial in q."""
    n, k = int(n), int(k)
    if k < 0 or k > n:
        return sp.Integer(0)
    num = sp.Integer(1)
    for i in range(k):
        num *= qint(n - i)
    return sp.expand(sp.cancel(num / qfact(k)))

def phi2_mult(poly):
    """multiplicity of (q+1) in poly  =  order of vanishing at q=-1."""
    poly = sp.expand(poly)
    if poly == 0:
        return None
    m = 0
    while True:
        val = poly.subs(q, -1)
        if sp.simplify(val) != 0:
            return m
        poly = sp.cancel(poly / (q + 1))
        m += 1

def v2(n):
    n = abs(int(n)); k = 0
    if n == 0: return 10**9
    while n % 2 == 0:
        n //= 2; k += 1
    return k

print("="*72)
print("Route-1 q-lift stub: what ord_{q=-1} can and cannot see (corrected)")
print("="*72)
from math import comb

# FACT 1 (corrected -- my first guess that ord_{q=-1}[n,k]_q == v2(C(n,k)) was FALSE).
# The Phi_2-multiplicity of a SINGLE Gaussian binomial is at most 1:
#    ord_{q=-1}[n,k]_q = floor(n/2)-floor(k/2)-floor((n-k)/2)  =  1 iff (n even & k odd) else 0.
print("\n FACT 1: ord_{q=-1}[n,k]_q vs the closed form floor(n/2)-floor(k/2)-floor((n-k)/2):")
bad1 = 0; in01 = True
for n in range(1, 18):
    for k in range(0, n + 1):
        qm = phi2_mult(qbinom(n, k))
        cf = n//2 - k//2 - (n-k)//2
        if qm != cf: bad1 += 1
        if qm not in (0,1): in01 = False
print(f"   closed-form match: mismatches={bad1};  always in {{0,1}}: {in01}")
print(f"   => a single q-binomial's (q+1)-mult is BOUNDED BY 1 -- it is NOT the v2 content")
print(f"      (e.g. v2(C(4,2))=1 but ord_q=-1[4,2]=0; v2(C(8,4))=1 ord=0).")

# FACT 2: q-factorial carries the leading layer only.
#    ord_{q=-1}[m]_q!  =  floor(m/2)      (= #even factors),
#    whereas the true 2-adic valuation is  v2(m!) = m - s2(m) = sum_r floor(m/2^r).
print("\n FACT 2: ord_{q=-1}[m]_q!  vs  floor(m/2)  vs  v2(m!):")
print(f"   {'m':>3} {'ord_q=-1':>9} {'floor(m/2)':>10} {'v2(m!)':>7}")
for m in [4,5,6,8,12,16]:
    om = phi2_mult(qfact(m))
    print(f"   {m:>3} {om:>9} {m//2:>10} {m - bin(m).count('1'):>7}")
print("   => ord_{q=-1} recovers ONLY the leading binary layer floor(m/2^1), NOT v2.")

print("""
STRUCTURAL VERDICT (a clean NEGATIVE with a reason)
---------------------------------------------------
1. ord_{q=-1} (the Phi_2 = q+1 multiplicity) sees ONLY the lowest binary layer:
   for a q-binomial it is in {0,1}; for a q-factorial it is floor(m/2), i.e. the
   single term floor(m/2^1) of v2(m!) = sum_r floor(m/2^r).  So evaluating the
   GK q-quotient at q=-1 alone CANNOT reproduce contents of size 5..8 that the
   master deficits actually have -- those need the higher binary layers.

2. M_{b+i} is moreover a signed S_3 ALTERNANT, and the Job-A Kummer census found
   CANCELLATION-BORN content (c=7,k=5: alternating sum lifts v2(M) by +3 over the
   minimal term).  A termwise q-lift's (q+1)-mult is min-of-parts and misses this.

3. WHAT THE ROUTE ACTUALLY NEEDS (browse read of 2502.06032, Gatzweiler-Krattenthaler):
   - the boundary minor M_{b+i} as a SINGLE q-determinant (LGV q-lift of the 3x3
     Jacobi-Trudi matrix), so there is no alternating cancellation to hide content;
   - and then the content must be read as the FULL 2-adic valuation of the q=1
     value, i.e. via the tower of cyclotomics Phi_{2^r} at 2^r-th roots of unity
     (r=1,2,3,...), NOT just Phi_2 at q=-1.  The floor 2*floor(k/2) is the r=1
     contribution; the surplus the Job-A table shows (content > 2*floor(k/2)) is
     exactly the r>=2 layers.

STATUS: machinery verified; the q=-1-only route is provably insufficient for the
        full content (it captures the floor's leading layer only).  A genuine
        q-proof needs the GK q-determinant AND the cyclotomic tower, not q=-1 alone.
""")
