"""
jobB_route1_qdet.py -- 2026-06-19 CODE, Job B route-1 q-determinant probe.

Question (07-05 crown route 1): does the boundary deficit N_i^{(c)} have a
q-analogue whose FULL cyclotomic-tower multiplicity  Sum_r mult_{Phi_{2^r}}
equals its exact 2-adic content g[c][k] -- including the r>=2 surplus that a
single q-binomial (q=-1 only) misses (06-18 stub)?

KEY THEORY.  For any q-polynomial P(q) in Z[q] with P(1)!=0,
    Sum_{r>=1} mult_{Phi_{2^r}}(P)  =  v2(P(1))  -  leak(P),
where leak(P) >= 0 is the 2-content of P's NON-cyclotomic factors / leading coeff
at q=1.  Because every q-INTEGER [n]_q is a product of cyclotomics and
Phi_{2^r}(1)=2 (Phi_d(1)=1 for non-prime-power d, =p odd for odd-prime powers),
a q-lift that is a genuine RATIO OF q-INTEGERS (q-hook-length formula) has
leak=0: its tower multiplicity is EXACTLY v2 of the q=1 value.  A SIGNED
DETERMINANT q-lift can leak.  So "tower == content" tests whether M_{b+i}'s
2-content is FULLY CYCLOTOMIC.

CONSTRUCTION.  M_j = f^{(a-j,b-j,c-j)} (dimension continuation; M_j vs MN verified
in the interior).  q-analogue = q-fake-degree determinant
    D_q(mu) = [n]_q! * det( 1/[mu_i - i + j]_q! ),   n=|mu|,  [neg]! := inf.
D_q(1) = f^mu = M_j.  We compute tower(D_q) and assemble
    tower(N^q) = tower(D_q) + tower([c]_q![k]_q!) - tower(den_q)
and compare to the exact integer content g[c][k].
"""
import sympy as sp
from math import factorial

q = sp.symbols('q')

# ====================== alternant engine (self-contained) ======================
def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 10**9
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def _f(n):
    return factorial(n)

def Ecoef(j, A, B, C):
    if A < 0 or B < 0 or C < 0 or (A + B + C) != 2 * j:
        return 0
    s = A + B - C
    if s & 1:
        return 0
    al = s // 2; be = B - al; ga = A - al
    if al < 0 or be < 0 or ga < 0 or al + be + ga != j:
        return 0
    return _f(j) // (_f(al) * _f(be) * _f(ga))

def _mult3(n, i, j, k):
    if i < 0 or j < 0 or k < 0 or i + j + k != n:
        return 0
    return _f(n) // (_f(i) * _f(j) * _f(k))

def Tfast(p, qq, r, j, s):
    if p < 0 or qq < 0 or r < 0 or p + qq + r != 2 * j + s:
        return 0
    tot = 0
    for P in range(s + 1):
        for Q in range(s - P + 1):
            R = s - P - Q
            e = Ecoef(j, p - P, qq - Q, r - R)
            if e:
                tot += _mult3(s, P, Q, R) * e
    return tot

_S3 = [((2,1,0),1),((1,2,0),-1),((2,0,1),-1),((0,2,1),1),((1,0,2),1),((0,1,2),-1)]

def Malt(a, b, c, j, m):
    s = 2 * m - 2 * j
    base = (a + 2, b + 1, c)
    tot = 0
    for (perm, sgn) in _S3:
        tot += sgn * Tfast(base[0]-perm[0], base[1]-perm[1], base[2]-perm[2], j, s)
    return tot

def Ni_data(c, i, Pv, Bv, even):
    """Exact integer N_i, plus (M, den, mu=(a-j,b-j,c-j), n) on the slice."""
    k = c - i
    bv = 2 * Bv if even else 2 * Bv + 1
    av = 2 * Pv + bv + c
    m = Pv + bv + c
    j = bv + i
    if bv < c + 1 or j > m:
        return None
    M = Malt(av, bv, c, j, m)
    den = (bv - c + 1)
    for t in range(2, i + 1):
        den *= (bv + t)
    num = M * factorial(c) * factorial(k)
    if den == 0 or num % den != 0:
        return None
    mu = (av - j, bv - j, c - j)
    return num // den, M, den, mu, av, bv, m, j

# ====================== q-machinery ======================
def qint(n):
    n = int(n)
    return sp.Integer(0) if n <= 0 else sp.expand(sum(q**i for i in range(n)))

def qfact(n):
    n = int(n)
    r = sp.Integer(1)
    for i in range(1, n + 1):
        r *= qint(i)
    return sp.expand(r)

_PHI = {r: sp.Poly(sp.cyclotomic_poly(2**r, q), q) for r in range(1, 14)}

def phi2r_mult(polyP, r):
    """multiplicity of Phi_{2^r} in poly P (a sympy Poly in q)."""
    Phi = _PHI[r]
    m = 0
    while True:
        quo, rem = sp.div(polyP, Phi, q)
        if not rem.is_zero:
            return m
        polyP = quo
        m += 1

def tower_mult(poly, RMAX=13):
    poly = sp.expand(poly)
    if poly == 0:
        return None
    P = sp.Poly(poly, q)
    return sum(phi2r_mult(P, r) for r in range(1, RMAX + 1))

def Dq(mu):
    """q-fake-degree det: [n]_q! det(1/[mu_i-i+j]_q!).  Returns sympy poly in q."""
    L = len(mu); n = sum(mu)
    rows = []
    for i in range(L):
        row = []
        for jc in range(L):
            e = mu[i] - (i + 1) + (jc + 1)
            row.append(sp.Integer(0) if e < 0 else 1 / qfact(e))
        rows.append(row)
    det = sp.Matrix(rows).det()
    return sp.expand(sp.cancel(qfact(n) * det))

# ====================== SECTION 0: interior sanity (zero leakage) ======================
print("="*78)
print("SECTION 0: q-hook tower == v2  for genuine dimensions (zero leakage check)")
print("="*78)
bad0 = 0
for mu in [(3,2,1),(4,2,2),(5,3,1),(6,4,2),(7,4,3),(8,5,2),(6,5,4)]:
    D = Dq(mu); f1 = int(D.subs(q, 1)); tw = tower_mult(D); v = v2int(f1)
    ok = (tw == v)
    bad0 += (0 if ok else 1)
    print(f"  mu={mu}: D_q(1)={f1}, tower={tw}, v2={v}  {'OK' if ok else 'LEAK'}")
print(f"  anomalies={bad0}  (tower==v2 confirms the cyclotomic-tower principle)\n")

# ====================== SECTION 1: boundary deficits -- tower vs content ======================
print("="*78)
print("SECTION 1: boundary N_i^{(c)} -- does tower(N^q) == content g[c][k] ?")
print("  tower(N^q) = tower(D_q(mu)) + tower([c]!_q[k]!_q) - tower(den_q)")
print("  den_q = [b-c+1]_q * prod_{t=2}^i [b+t]_q  (ratio of q-integers, leak=0)")
print("="*78)

def content_and_qtower(c, i, even, Bmax=40, Pmax=40):
    """find content minimizer over the slice; at it, compute the q-tower of N^q."""
    k = c - i
    g = 10**9; best = None
    for Bv in range((c + 1) // 2 + 1, Bmax + 1):
        for Pv in range(0, Pmax + 1):
            res = Ni_data(c, i, Pv, Bv, even)
            if res is None:
                continue
            val = res[0]; vv = v2int(val)
            if vv < g:
                g = vv; best = (Pv, Bv, res)
    if best is None:
        return None
    Pv, Bv, res = best
    Nval, M, den, mu, av, bv, m, j = res
    # q-tower of M via D_q(mu).  mu may have negative parts (boundary): Dq handles e<0.
    Dpoly = Dq(mu)
    D1 = int(Dpoly.subs(q, 1))
    den_q = qint(bv - c + 1)
    for t in range(2, i + 1):
        den_q *= qint(bv + t)
    ck_q = qfact(c) * qfact(k)
    tw_D = tower_mult(Dpoly)        # None if Dpoly==0
    tw_den = tower_mult(den_q)
    tw_ck = tower_mult(ck_q)
    tw_N = None if tw_D is None else tw_D + tw_ck - tw_den
    return dict(g=g, Pv=Pv, Bv=Bv, mu=mu, M=M, D1=D1, Nval=Nval,
                tw_D=tw_D, tw_ck=tw_ck, tw_den=tw_den, tw_N=tw_N,
                v2M=v2int(M), leakM=(None if tw_D is None else v2int(M) - tw_D))

print(f"  {'c':>2} {'k':>2} {'i':>2} {'sl':>4} | {'g':>3} {'twN':>4} | "
      f"{'v2M':>3} {'twD':>3} {'leakM':>5} | {'D1==M?':>7} | match")
rows = []
for c in range(4, 10):
    for k in (4, 5, 6):
        i = c - k
        if i < 1:
            continue
        for even in (True, False):
            r = content_and_qtower(c, i, even)
            if r is None:
                continue
            sl = 'even' if even else 'odd'
            d1m = (abs(r['D1']) == abs(r['M']))
            match = (r['tw_N'] is not None and r['tw_N'] == r['g'])
            rows.append((c, k, i, sl, r, d1m, match))
            twN = 'None' if r['tw_N'] is None else r['tw_N']
            twD = 'None' if r['tw_D'] is None else r['tw_D']
            lk = 'None' if r['leakM'] is None else r['leakM']
            print(f"  {c:>2} {k:>2} {i:>2} {sl:>4} | {r['g']:>3} {str(twN):>4} | "
                  f"{r['v2M']:>3} {str(twD):>3} {str(lk):>5} | "
                  f"{str(d1m):>7} | {'YES' if match else 'no'}    "
                  f"D1={r['D1']} M={r['M']}")

nmatch = sum(1 for *_, m in rows if m)
nD1 = sum(1 for *_, d, _ in rows if d)
nNone = sum(1 for *_, r, _, _ in rows if r['tw_D'] is None)
print(f"\n  q=1 check |D_q(1)|==|M|: {nD1}/{len(rows)} pass  "
      f"(D_q(mu) reproduces the boundary alternant?)")
print(f"  D_q(mu) == 0 while M != 0: {nNone}/{len(rows)} cases (naive q-hook continuation FAILS)")
print(f"  tower(N^q)==g matches: {nmatch}/{len(rows)}")
