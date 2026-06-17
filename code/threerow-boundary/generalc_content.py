"""
generalc_content.py -- prove session 2026-06-17

Goal: extract the deficit polynomials N_i^{(c)} for the three-row boundary
indices j=b+i (i=1..c), and compute their 2-adic CONTENT
   g = min over integer (P>=0, b) of v2(N_i(a=2P+b+c, b)),
split by parity of b.  Fit g(k) as a function of depth k = c - i, uniform in c.

Convention (matches c=3,c=4 notes):
   N_i^{(c)} := a-involving part of M_{b+i}, with the HOOK factors
               (a-b+1) and (a-c+2) removed (they cancel in R_i = G_{b+i}/G_0).
   top index i=c has N_c = 1 (depth 0).
"""
import sympy as sp
from math import comb, factorial
from mn import Mj

a, b, P = sp.symbols('a b P')

def collect(c, i, amax=60, bmax=26):
    pts = []
    for bb in range(c, bmax):
        for aa in range(bb, amax):
            if (aa + bb + c) % 2 != 0:
                continue
            m = (aa + bb + c) // 2
            j = bb + i
            if j > m:
                continue
            pts.append((aa, bb, Mj((aa, bb, c), j, m)))
    return pts

def fit_poly(pts, maxdeg):
    monos = [(p, q) for p in range(maxdeg + 1) for q in range(maxdeg + 1)
             if p + q <= maxdeg]
    if len(pts) < len(monos) + 3:
        return None
    A = sp.Matrix([[aa**p * bb**q for (p, q) in monos] for (aa, bb, _) in pts])
    y = sp.Matrix([M for (_, _, M) in pts])
    sol = A.solve_least_squares(y)
    if any(A * sol - y):
        return None
    expr = sum(sol[k] * a**p * b**q for k, (p, q) in enumerate(monos))
    return sp.factor(sp.expand(expr))

def extract_deficit(c, i):
    """fit M_{b+i}, factor, strip hook factors (a-b+1),(a-c+2); return N_i (a-involving residue)."""
    deg = 2 * c + 1
    M = fit_poly(collect(c, i), maxdeg=deg)
    if M is None:
        return None, None
    cst, facs = sp.factor_list(M)
    hookA = sp.expand(a - b + 1)
    hookB = sp.expand(a - c + 2)
    N = sp.Integer(1)
    bprod = sp.Integer(1)
    removedA = removedB = 0
    for fac, mult in facs:
        fac_e = sp.expand(fac)
        deg_a = sp.Poly(fac_e, a).degree() if fac_e.has(a) else 0
        if deg_a == 0:
            bprod *= fac**mult
            continue
        # a-involving factor
        if removedA == 0 and sp.simplify(fac_e - hookA) == 0:
            removedA = 1
            mult -= 1
        if removedB == 0 and sp.simplify(fac_e - hookB) == 0:
            removedB = 1
            mult -= 1
        if mult > 0:
            N *= fac**mult
    N = sp.expand(N)
    return M, N

def v2int(n):
    n = abs(int(n))
    if n == 0:
        return 999
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

def content_v2(N, c, parity):
    """min v2(N(a=2P+b+c, b)) over a grid, b of given parity."""
    g = 999
    witness = None
    for Pv in range(0, 40):
        for bv in range(c, 80):
            if parity is not None and bv % 2 != parity:
                continue
            av = 2 * Pv + bv + c
            val = int(N.subs({a: av, b: bv}))
            vv = v2int(val)
            if vv < g:
                g = vv; witness = (Pv, bv, av, val)
    return g, witness

print("="*80)
print("DEFICIT POLYNOMIALS N_i^{(c)} and their 2-adic contents")
print("k = depth = c - i ;  i=c is top (N=1, k=0)")
print("="*80)
table = {}
for c in range(2, 6):
    print("\n" + "#"*60 + "  c = %d" % c)
    for i in range(1, c + 1):
        k = c - i
        M, N = extract_deficit(c, i)
        if N is None:
            print("  i=%d (k=%d): NO FIT" % (i, k))
            continue
        adeg = sp.Poly(sp.expand(N), a).degree() if N.has(a) else 0
        g_even, w_e = content_v2(N, c, 0)
        g_odd, w_o = content_v2(N, c, 1)
        g_all = min(g_even, g_odd)
        table[(c, k)] = (g_all, g_even, g_odd)
        print("  i=%d (k=%d): a-deg(N)=%d  content g=%d  [b even: %d, b odd: %d]"
              % (i, k, adeg, g_all, g_even, g_odd))
        if c <= 5:
            print("         N = %s" % N)

print("\n" + "="*80)
print("TABLE of content g[c][k]  (min over parities)")
print("="*80)
print("      " + "".join("k=%-7d" % k for k in range(0, 6)))
for c in range(2, 6):
    row = "c=%-3d " % c
    for k in range(0, 6):
        if (c, k) in table:
            row += "%-9d" % table[(c, k)][0]
        else:
            row += "%-9s" % "."
    print(row)

print("\nTABLE  [b even / b odd] split:")
for c in range(2, 6):
    for k in range(0, 6):
        if (c, k) in table:
            ga, ge, go = table[(c, k)]
            print("  c=%d k=%d : g=%d  (even %d / odd %d)" % (c, k, ga, ge, go))
