"""CROSS-CHECKS against the two code-disjoint engines of the 2026-09-08 self-review.
   (A) side: ribbon.py  -- border strips on Young diagrams, no abacus.
   (B) side: hbasis.py  -- h-basis dictionaries, no power sums.
Neither engine was written this session and neither knows anything about Q105."""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/reviews/2026-09-08-selfreview-code')
import ribbon as RB
import hbasis as HB
from abacus import kappa, to_beads, parts as aparts

# ---- (A): every nonzero one-bead element of [R_e(t),R_f(1/t)] has ord_{t=-1} == 1,
#           and it is zero exactly when kappa == 0.   Computed on ribbon.py's engine.
t = RB.t; s = RB.s
nz = zer = bad = 0; orders = {}
for n in range(0, 6):
    for lam in RB.parts(n):
        for e in range(1, 5):
            for f in range(1, 5):
                out = RB.comm(lam, e, f)
                L = len(lam) + 2*(e+f) + 4
                M0 = to_beads(lam, L)
                for mu, val in out.items():
                    Mp = to_beads(mu, L)
                    if len(M0 ^ Mp) != 2:
                        continue
                    (aa,) = M0 - Mp
                    v = sp.cancel(sp.together(val.subs(s, 1/t)))
                    k, N = kappa(M0, aa, e, f)
                    if v == 0:
                        zer += 1
                        if k != 0: bad += 1; print("XA-FAIL zero but kappa!=0", lam, e, f, aa, k)
                        continue
                    if k == 0: bad += 1; print("XA-FAIL nonzero but kappa==0", lam, e, f, aa, v)
                    num = sp.Poly(sp.numer(v), t)
                    o = 0
                    while num.eval(-1) == 0 and num.degree() > 0:
                        o += 1; num = num.diff(t)
                    orders[o] = orders.get(o, 0)+1
                    c1 = sp.limit(sp.diff(v, t), t, -1)
                    if sp.simplify(c1 - (-1)**N * k) != 0:
                        bad += 1; print("XA-FAIL lead", lam, e, f, aa, v, c1, (-1)**N*k)
                    nz += 1
print("[A] ribbon.py cross-check: nonzero=%d order-histogram=%s  zero=%d  failures=%d"
      % (nz, orders, zer, bad))

# ---- (B): Q99 Thm B on hbasis.py, then its order of vanishing at t=-1.
th = HB.sp.symbols('t')
def hb_defect(m, n, f):
    A = HB.mode(HB.mode(f, 1/th, n), th, m)
    B = HB.mode(HB.mode(f, th, m), 1/th, n)
    return HB.sf_add(A, HB.sf_scale(B, -1/th))

def ordv(p):
    p = sp.Poly(sp.numer(sp.cancel(sp.together(p))), th)
    if p.is_zero: return None
    o = 0
    while p.eval(-1) == 0:
        o += 1; p = p.diff(th)
        if p.is_zero: return None
    return o

hist = {}
for m in range(0, 5):
    for n in range(0, 5):
        D = hb_defect(m, n, HB.sf_one())
        N = m+n
        os_ = [ordv(v) for v in D.values() if sp.cancel(v) != 0]
        if not os_: continue
        hist.setdefault(N, set()).update(os_)
print("[B] hbasis.py cross-check, ord_{t=-1} of the defect's h-basis coefficients:")
for N in sorted(hist):
    print("     m+n=%d -> orders %s   (N mod 2 = %d)" % (N, sorted(hist[N]), N % 2))
