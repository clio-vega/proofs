"""THE DECISIVE TEST (PROVE 2026-09-19 c2).

Theorem to test.  Let Gamma = {(alpha, m*1 - alpha)} subset Z^{2n} and let
S subset Gamma be finite nonempty with x-projection A subset {0..m}^n.  Then

    S is M-convex   <==>   A is an integer BOX  prod_i [l_i, u_i].

Three tests, in this order:
  (0) PURE test, independent of any vertex model: random abstract A.
  (1) MODEL test on the 769 non-monomial Z of structural.py (untuned: those
      supports were computed yesterday for a different purpose).
  (2) The brief's refutation test: does any M-convex passer have constant Sigma?
  (3) Positive control from the source: five-vertex Schur, where lem:rowdeg does
      NOT apply -- support must be off Gamma and M-convex (HMMS Thm 2).
Also instrumented: does the equal-total-degree guard inside is_M_convex ever fire?
"""
import sympy as sp, itertools, random
from fractions import Fraction
from sixvertex import Z, schur_weights
from lorentzian import is_M_convex
import lorentzian

# ---- instrument the vacuous-guard question ---------------------------------
GUARD_FIRED = [0]
_orig = lorentzian.is_M_convex
def is_M_convex_traced(support):
    J = set(support)
    if J and len({sum(a) for a in J}) != 1:
        GUARD_FIRED[0] += 1
    return _orig(J)

def is_box(A):
    """A subset of Z^n is an integer box iff it equals the product of its
    coordinate-wise min..max intervals."""
    A = set(A)
    n = len(next(iter(A)))
    lo = [min(a[i] for a in A) for i in range(n)]
    hi = [max(a[i] for a in A) for i in range(n)]
    size = 1
    for i in range(n):
        size *= hi[i] - lo[i] + 1
    if size != len(A):
        return False
    return all(all(lo[i] <= a[i] <= hi[i] for i in range(n)) for a in A)

def cdict(expr, vars_):
    p = sp.Poly(sp.expand(expr), *vars_)
    return {tuple(mm): Fraction(int(c)) for mm, c in zip(p.monoms(), p.coeffs())}

# ============================================================================
print("=== (0) PURE TEST: random abstract A on Gamma, no vertex model ===")
random.seed(2026)
pure_tot = pure_bad = 0
pure_pos = pure_neg = 0
for trial in range(4000):
    n = random.choice([1, 2, 2, 3])
    m = random.choice([1, 2, 3])
    pts = list(itertools.product(range(m + 1), repeat=n))
    k = random.randint(1, min(len(pts), 6))
    A = set(random.sample(pts, k))
    S = {tuple(a) + tuple(m - ai for ai in a) for a in A}
    mc = is_M_convex_traced(S)
    bx = is_box(A)
    pure_tot += 1
    if mc != bx:
        pure_bad += 1
        if pure_bad <= 5:
            print("   MISMATCH n=%d m=%d A=%s  M-convex=%s box=%s" % (n, m, sorted(A), mc, bx))
    if bx: pure_pos += 1
    else:  pure_neg += 1
print("  random A tested            : %d   (box: %d, non-box: %d)" % (pure_tot, pure_pos, pure_neg))
print("  M-convex(S) != isbox(A)    : %d   <-- must be 0" % pure_bad)

# also: exhaustive for small n,m
print("\n  exhaustive over ALL nonempty A for (n,m) in {(1,1),(1,2),(1,3),(2,1),(2,2)}:")
ex_tot = ex_bad = 0
for (n, m) in [(1,1),(1,2),(1,3),(2,1),(2,2)]:
    pts = list(itertools.product(range(m + 1), repeat=n))
    for r in range(1, len(pts) + 1):
        for A in itertools.combinations(pts, r):
            S = {tuple(a) + tuple(m - ai for ai in a) for a in A}
            ex_tot += 1
            if is_M_convex_traced(S) != is_box(set(A)):
                ex_bad += 1
                if ex_bad <= 5:
                    print("     MISMATCH n=%d m=%d A=%s" % (n, m, A))
print("    subsets tested: %d   mismatches: %d   <-- must be 0" % (ex_tot, ex_bad))

# ============================================================================
print("\n=== (1) MODEL TEST: the 769 non-monomial Z of structural.py ===")
random.seed(11)   # EXACTLY structural.py's stream
tot = mono = nontriv = notMconv = xlevel1 = 0
agree = disagree = 0
box_and_mconv = nonbox_and_notmconv = 0
passers_constant_sigma = []
failers_constant_sigma = 0
ffcount = 0
for trial in range(60):
    n, m = random.choice([(2,2),(2,3),(3,2),(2,4)])
    xs = sp.symbols('x1:%d' % (n+1)); ys = sp.symbols('y1:%d' % (n+1))
    pat = [(random.randint(0,3), random.randint(0,3)) for _ in range(6)]
    pat = [(1,0) if p == (0,0) else p for p in pat]
    def mk(i):
        xi, yi = xs[i], ys[i]
        nm = ['a1','a2','b1','b2','c1','c2']
        return {k: (p[0]*xi + p[1]*yi) for k, p in zip(nm, pat)}
    w = [mk(i) for i in range(n)]
    if sp.expand(w[0]['a1']*w[0]['a2'] + w[0]['b1']*w[0]['b2'] - w[0]['c1']*w[0]['c2']) == 0:
        ffcount += 1
    vars_ = list(xs) + list(ys)
    for top in itertools.product([0,1], repeat=m):
        for bot in itertools.product([0,1], repeat=m):
            z = Z(n, m, list(top), list(bot), [0]*n, [0]*n, w)
            if z == 0: continue
            cd = cdict(z, vars_); tot += 1
            for a in cd:
                for i in range(n):
                    assert a[i] + a[n+i] == m
            if len(cd) == 1: mono += 1; continue
            nontriv += 1
            mc = is_M_convex_traced(set(cd))
            A = {tuple(a[:n]) for a in cd}
            bx = is_box(A)
            sigma_const = (len({sum(t) for t in A}) == 1)
            if sigma_const: xlevel1 += 1
            if not mc: notMconv += 1
            if mc == bx:
                agree += 1
                if mc: box_and_mconv += 1
                else:  nonbox_and_notmconv += 1
            else:
                disagree += 1
                if disagree <= 5:
                    print("   MISMATCH pat=%s n=%d m=%d A=%s M-conv=%s box=%s"
                          % (pat, n, m, sorted(A), mc, bx))
            if sigma_const:
                if mc: passers_constant_sigma.append((pat, n, m, sorted(A)))
                else:  failers_constant_sigma += 1
print("  free-fermion families in stream : %d" % ffcount)
print("  nonzero Z / monomial / non-monomial : %d / %d / %d" % (tot, mono, nontriv))
print("  not M-convex                   : %d" % notMconv)
print("  A is NOT a box                 : %d" % (nontriv - box_and_mconv - 0 if False else nontriv - (agree - nonbox_and_notmconv) - disagree + 0))
print("  PREDICTION  M-convex <=> A box : agree %d / %d   DISAGREE %d   <-- must be 0"
      % (agree, nontriv, disagree))
print("     (M-convex AND box)          : %d" % box_and_mconv)
print("     (not M-convex AND not box)  : %d" % nonbox_and_notmconv)

# ============================================================================
print("\n=== (2) REFUTATION TEST (from the brief): constant-Sigma passers ===")
print("  non-monomial Z with constant Sigma_i deg_{x_i} : %d / %d" % (xlevel1, nontriv))
print("    of which NOT M-convex : %d" % failers_constant_sigma)
print("    of which     M-convex : %d   <-- (T2) is REFUTED if nonzero" % len(passers_constant_sigma))
for p in passers_constant_sigma: print("      ", p)

# ============================================================================
print("\n=== (3) POSITIVE CONTROL from the source: five-vertex Schur ===")
print("    HMMS 1906.09633 Thm 2: N(s_lambda) is Lorentzian, hence M-convex.")
print("    lem:rowdeg does NOT apply (four weights have degree 0), so supp is OFF Gamma.")
xv = sp.symbols('x1:6')
for (n, m) in [(2,4),(2,5),(3,5)]:
    ws = [schur_weights(xv[i]) for i in range(n)]
    vars_ = list(xv[:n])
    for topset in itertools.combinations(range(m), n):
        top = [1 if j in topset else 0 for j in range(m)]
        for botset in itertools.combinations(range(m), n):
            bot = [1 if j in botset else 0 for j in range(m)]
            z = Z(n, m, top, bot, [0]*n, [0]*n, ws)
            if z == 0: continue
            cd = cdict(z, vars_)
            if len(cd) == 1: continue
            mc = is_M_convex(set(cd))
            offGamma = len({sum(a) for a in cd}) == 1  # homogeneous, but in ONE alphabet
            print("    n=%d m=%d top=%s bot=%s  |supp|=%d  M-convex=%s"
                  % (n, m, topset, botset, len(cd), mc))
            assert mc, "INSTRUMENT IS WRONG: Schur control failed"

print("\n=== equal-total-degree guard inside is_M_convex ===")
print("  times the guard `len(degs)!=1` FIRED : %d" % GUARD_FIRED[0])
print("  (on Gamma every point has total degree n*m, so on the model test it is VACUOUS)")
