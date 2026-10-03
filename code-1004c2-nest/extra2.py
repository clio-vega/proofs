import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen, diffsys
from itertools import product
from collections import Counter
st = Counter()
for m in range(2, 6):
    for n in range(m, m+5):
        for (mu, lam) in gen.pairs(n, m, 7):
            a_, g, P, Q = diffsys.data(mu, lam, n, m)
            for b in range(0, 2*n+1):
                sig = sum(mu) + b - sum(a_)
                Y = [y for y in product(*[range(P[i], Q[i]+1) for i in range(m)]) if sum(y) == sig]
                Ye = [y for y in Y if all(g[i]+1-max(y[i],0)-max(-y[(i+1)%m],0) >= 1 for i in range(m))]
                if not Ye: continue
                st['slices'] += 1
                nt = len(Ye) >= 2
                st['nontrivial'] += nt
                up = all(min(y) >= 0 for y in Ye)
                dn = all(max(y) <= 0 for y in Ye)
                if up: st['y_ge_0'] += 1; st['y_ge_0_nt'] += nt
                if dn: st['y_le_0'] += 1; st['y_le_0_nt'] += nt
                if up or dn: st['covered'] += 1; st['covered_nt'] += nt
                # is the width set a BOX SLICE?  (necessary: single half-width)
                W = [tuple(g[i]+1-max(y[i],0)-max(-y[(i+1)%m],0) for i in range(m)) for y in Ye]
                hw = {sum(w) for w in W}
                if len(hw) == 1:
                    st['single_hw'] += 1; st['single_hw_nt'] += nt
                    # box slice test: W == box(min,max) cap {sum = T}
                    lo = [min(w[i] for w in W) for i in range(m)]
                    hi = [max(w[i] for w in W) for i in range(m)]
                    T = sum(W[0])
                    full = [w for w in product(*[range(lo[i], hi[i]+1) for i in range(m)]) if sum(w) == T]
                    if set(full) == set(W):
                        st['W_is_boxslice'] += 1; st['W_is_boxslice_nt'] += nt
                        if not (up or dn): st['boxslice_NOT_covered'] += 1; st['boxslice_NOT_covered_nt'] += nt
print(dict(st))
