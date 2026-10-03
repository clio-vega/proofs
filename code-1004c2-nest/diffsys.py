"""THE DIFFERENCE-SYSTEM MODEL for condition (A).

Derivation (all from the definitions, no new input):
  bead coords a_i = lam_{i-1}+1, g_i = lam_i - lam_{i-1} - 1 >= 0, sum g_i = G = n-m.
  y_i := nu_i - a_i      r_i := kappa_i - a_i
  mu -< nu  h-strip  <=>  mu_i <= nu_i <= mu_{i+1}-1            (a BOX in y)
  nu -< kappa h-strip <=> nu_i <= kappa_i < nu_{i+1}            <=> y_i <= r_i <= g_i + y_{i+1}
  kappa -< lam h-strip <=> lam_{i-1} < kappa_i <= lam_i         <=> 0 <= r_i <= g_i
  b = S_nu - |mu| ;  sigma := sum y_i = |mu| + b - sum a_i
  a = |kappa/nu| = sum (r_i - y_i) = sum r_i - sigma

  ==>  k(a,b) = # { (y,r) in Z^m x Z^m :
                      P_i <= y_i <= Q_i,            sum y_i = sigma,
                      0   <= r_i <= g_i,            sum r_i = a + sigma,
                      y_i <= r_i,   r_i - y_{i+1} <= g_i   (cyclic) }

Every constraint is a BOUND ON A SINGLE COORDINATE or a BOUND ON A DIFFERENCE OF
TWO COORDINATES.  So the ambient set is the lattice-point set of an ALCOVED
polytope and is closed under componentwise max and min.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from gen import lam_prev, nu_next
from itertools import product
from collections import Counter

def data(mu, lam, n, m):
    a = [lam_prev(lam, i, n, m) + 1 for i in range(m)]
    g = [lam[i] - lam_prev(lam, i, n, m) - 1 for i in range(m)]
    P = [mu[i] - a[i] for i in range(m)]
    Q = [nu_next(mu, i, n, m) - 1 - a[i] for i in range(m)]
    return a, g, P, Q

def k_diffsys(mu, lam, n, m, b):
    """{a : k(a,b)} by direct enumeration of the difference system."""
    a_, g, P, Q = data(mu, lam, n, m)
    sig = sum(mu) + b - sum(a_)
    out = Counter()
    for y in product(*[range(P[i], Q[i]+1) for i in range(m)]):
        if sum(y) != sig: continue
        lo = [max(0, y[i]) for i in range(m)]
        hi = [min(g[i], g[i] + y[(i+1) % m]) for i in range(m)]
        if any(lo[i] > hi[i] for i in range(m)): continue
        for r in product(*[range(lo[i], hi[i]+1) for i in range(m)]):
            out[sum(r) - sig] += 1
    return dict(out)

if __name__ == '__main__':
    st = Counter(); wit = []
    for m in range(2, 6):
        for n in range(m, m+5):
            for (mu, lam) in gen.pairs(n, m, 7):
                for b in range(0, 2*n+1):
                    off, co = gen.slice_sum(mu, lam, n, m, b)
                    ref = {off+j: c for j, c in enumerate(co) if c}
                    mine = {k: v for k, v in k_diffsys(mu, lam, n, m, b).items() if v}
                    st['slices'] += 1
                    if ref == mine: st['agree'] += 1
                    else:
                        st['DISAGREE'] += 1
                        if len(wit) < 4: wit.append((m, n, mu, lam, b, ref, mine))
    print("vs gen.slice_sum (box-slice formula):", dict(st))
    for w in wit: print("   ", w)

    # independent instrument: direct chain enumeration, no L_i/R_i, no y/r
    st2 = Counter()
    for m in range(2, 5):
        for n in range(m, m+4):
            for (mu, lam) in gen.pairs(n, m, 6):
                tab = gen.k_table_chains(mu, lam, n, m)
                bs = sorted({b for (a, b) in tab})
                for b in bs:
                    ref = {a: v for (a, bb), v in tab.items() if bb == b and v}
                    mine = {k: v for k, v in k_diffsys(mu, lam, n, m, b).items() if v}
                    st2['agree' if ref == mine else 'DISAGREE'] += 1
    print("vs gen.k_table_chains (independent chain enumeration):", dict(st2))

    # REFUSAL CONTROL: perturb one constraint, the model must stop agreeing
    def k_pert(mu, lam, n, m, b, which):
        a_, g, P, Q = data(mu, lam, n, m)
        sig = sum(mu) + b - sum(a_); out = Counter()
        for y in product(*[range(P[i], Q[i]+1) for i in range(m)]):
            if sum(y) != sig: continue
            if which == 'nocycle':   # r_i <= g_i + y_i  instead of y_{i+1}
                hi = [min(g[i], g[i] + y[i]) for i in range(m)]
            elif which == 'shift':   # r_i <= g_i + y_{i+1} + 1
                hi = [min(g[i], g[i] + y[(i+1) % m] + 1) for i in range(m)]
            lo = [max(0, y[i]) for i in range(m)]
            if any(lo[i] > hi[i] for i in range(m)): continue
            for r in product(*[range(lo[i], hi[i]+1) for i in range(m)]):
                out[sum(r) - sig] += 1
        return dict(out)
    for which in ('nocycle', 'shift'):
        st3 = Counter()
        for m in (2, 3, 4):
            for n in range(m, m+4):
                for (mu, lam) in gen.pairs(n, m, 6):
                    for b in range(0, 2*n+1):
                        off, co = gen.slice_sum(mu, lam, n, m, b)
                        ref = {off+j: c for j, c in enumerate(co) if c}
                        mine = {k: v for k, v in k_pert(mu, lam, n, m, b, which).items() if v}
                        if not ref and not mine: continue
                        st3['agree' if ref == mine else 'REFUSED'] += 1
        print(f"refusal control '{which}':", dict(st3))
