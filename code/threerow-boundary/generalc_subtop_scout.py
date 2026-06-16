"""
generalc_subtop_scout.py  -- Job B (2026-06-16 code session)

Does the c=3 boundary subtop structure persist uniformly in c?  We derive the
top and next-to-top closed forms M_{b+c}, M_{b+c-1} for c=4,5 (and re-confirm
c=2,3), test the conjectured uniform patterns

    TOP   :  M_{b+c}   = (b-c+1) * prod_{t=2}^{c}(b+t) / c!            (a-FREE hook, Lemma T)
    SUBTOP:  M_{b+c-1} = (b-c+1) * prod_{t=2}^{c-1}(b+t) * N_{c-1}^{(c)} / c!
             with  N_{c-1}^{(c)} = a(b+c) - ( b^2 + (c-1)b + c(c-2) ).

Then we check, numerically, whether the boundary lemma Delta(b+i) > -theta holds
for c=4,5 with the SAME shape of telescope (positive-linear + carries - 2 v2(N))
and whether the two-factor Lemma F2 (Q>=6, remove two designated members) is the
uniform keystone for all c.
"""
import sympy as sp
from math import comb, factorial
from mn import Mj

a, b = sp.symbols('a b')

# --------------------------------------------------------------------------
# 1. fit M_{b+i} as an exact integer polynomial in (a,b), high degree allowed
# --------------------------------------------------------------------------
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

def top_form(c):
    return sp.factor((b - c + 1) * sp.prod([b + t for t in range(2, c + 1)]) / factorial(c))

def subtop_form(c):
    N = a * (b + c) - (b**2 + (c - 1) * b + c * (c - 2))
    return sp.factor((b - c + 1) * sp.prod([b + t for t in range(2, c)]) * N / factorial(c))

print("=" * 72)
print("1. TOP and SUBTOP closed forms vs conjectured uniform pattern")
print("=" * 72)
for c in [2, 3, 4, 5]:
    # top
    fitted_top = fit_poly(collect(c, c), maxdeg=c + 1)
    conj_top = top_form(c)
    ok_top = (fitted_top is not None and sp.simplify(fitted_top - conj_top) == 0)
    # subtop
    deg = 2 * c          # generous
    fitted_sub = fit_poly(collect(c, c - 1), maxdeg=deg)
    conj_sub = subtop_form(c)
    ok_sub = (fitted_sub is not None and sp.simplify(fitted_sub - conj_sub) == 0)
    print("c=%d:" % c)
    print("   TOP    fitted = %s" % fitted_top)
    print("          match Lemma T pattern: %s" % ok_top)
    print("   SUBTOP fitted = %s" % fitted_sub)
    print("          match N_{c-1}^{(c)}=a(b+c)-(b^2+(c-1)b+c(c-2)) pattern: %s" % ok_sub)

# --------------------------------------------------------------------------
# 1b. the SECOND subtop M_{b+c-2}: where uniformity BREAKS
# --------------------------------------------------------------------------
print()
print("=" * 72)
print("1b. SECOND subtop M_{b+c-2} -- a-degree of the deficit (uniformity check)")
print("=" * 72)
for c in [3, 4, 5]:
    f = fit_poly(collect(c, c - 2, amax=60, bmax=28), maxdeg=2 * c + 1)
    if f is None:
        print("c=%d: no fit" % c); continue
    fl = sp.factor_list(f)[1]
    nlin = sum(mult for fac, mult in fl if sp.Poly(fac, a).degree() == 1)
    adeg = sp.Poly(sp.expand(f), a).degree()
    print("c=%d  M_{b+%d} = %s" % (c, c - 2, f))
    print("        a-degree=%d, #a-linear factors=%d  -> %s"
          % (adeg, nlin,
             "TWO a-linear (c=3 factor-in-product pattern)" if nlin >= 2
             else "IRREDUCIBLE a-quadratic deficit (pattern BREAKS for c>=4)"))

# --------------------------------------------------------------------------
# 2. the uniform N_{c-1}^{(c)} and its a-even 2-adic fact
# --------------------------------------------------------------------------
print()
print("=" * 72)
print("2. a-even 2-adic fact for the subtop deficit N_{c-1}^{(c)} (uniform in c?)")
print("=" * 72)
def v2(n):
    n = abs(n)
    if n == 0:
        return 99
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k
for c in [2, 3, 4, 5, 6]:
    dist = {}
    for m in range(c + 4, 80):
        for bb in range(c, m):
            aa = 2 * m - c - bb
            if aa < bb or aa % 2 != 0:
                continue
            if (aa + bb + c) % 2 != 0:
                continue
            N = aa * (bb + c) - (bb**2 + (c - 1) * bb + c * (c - 2))
            dist[v2(N)] = dist.get(v2(N), 0) + 1
    print("c=%d  a-even  v2(N_{c-1}^{(c)}) distribution: %s" % (c, dict(sorted(dist.items()))))

# --------------------------------------------------------------------------
# 3. boundary lemma holds for c=4,5 with margin (direct val sweep) + F2 keystone
# --------------------------------------------------------------------------
print()
print("=" * 72)
print("3. boundary lemma Delta(b+i) > -theta for c=4,5  (direct val, MN engine)")
print("=" * 72)
def val(lam, j, m):
    M = Mj(lam, j, m)
    if M == 0:
        return None
    return j + 2 * v2(comb(m, j) * M)

for c in [4, 5]:
    cnt = 0; bad = 0; minD = 99; minshape = None
    for m in range(c + 4, 42):
        for bb in range(c, m):
            aa = 2 * m - c - bb
            if aa < bb:
                continue
            P = m - bb - c
            if P < 0:
                continue
            # crude box-interior gate: require interior argmin box to fit; use bb>=2c
            if bb < 2 * c:
                continue
            theta = 0 if aa % 2 == 0 else (c * (c - 1)) // 2  # tau-style offset guess
            lam = (aa, bb, c)
            v0 = val(lam, 0, m)
            if v0 is None:
                continue
            for i in range(1, c + 1):
                j = bb + i
                if j > m:
                    continue
                vj = val(lam, j, m)
                if vj is None:
                    continue
                D = vj - v0
                cnt += 1
                if D <= minD:
                    minD = D; minshape = (aa, bb, m, i)
                if not (D > -theta):
                    bad += 1
    print("c=%d:  indices=%d  min Delta=%d at (a,b,m,i)=%s  violations(theta-guess)=%d"
          % (c, cnt, minD, minshape, bad))

# --------------------------------------------------------------------------
# 4. uniform two-factor Lemma F2 (Q>=6 threshold) -- already c-independent
# --------------------------------------------------------------------------
print()
print("=" * 72)
print("4. two-factor Lemma F2 is c-INDEPENDENT (pure statement on consec ints)")
print("=" * 72)
import itertools
fails = 0; checked = 0
for Q in range(6, 30):
    for R in range(0, 400):
        tot = sum(v2(R + t) for t in range(1, Q + 1))
        for t1, t2 in itertools.combinations(range(1, Q + 1), 2):
            checked += 1
            if tot - v2(R + t1) - v2(R + t2) < 1:
                fails += 1
print("F2 (Q>=6, any 2 designated members removed, rest even): checked=%d failures=%d"
      % (checked, fails))
print("=> F2 is a statement about consecutive integers ONLY; it is the SAME")
print("   keystone for every c.  The c-dependence lives entirely in (i) which")
print("   members are designated and (ii) the count Q = b//2 + (offset), which")
print("   grows with b, so Q>=6 holds throughout each box interior.")
