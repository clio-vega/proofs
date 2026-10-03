"""THE ONE-SIDED CASE -- exactly the regime where the BRIEFED hypothesis (H1) is
satisfiable, so (H2) is the whole question there.

If P_i >= 0 for all i (all y_i >= 0 on the slice) then (y_i)_+ = y_i and
(-y_{i+1})_+ = 0, so eq:w becomes LINEAR:   w_i = g_i + 1 - y_i.
Hence the width-vector set is
      W = { w : L_i <= w_i <= U_i,  sum w_i = T }     (L_i = g_i+1-Q_i, U_i = g_i+1-P_i,
                                                       T = G+m-sigma)
a BOX SLICE, i.e. an M-CONVEX set -- and every half-width equals (T-m)/2, a single
value, exactly as Theorem R of [1004] demands.  So:

   (Q)  Is  F(z) = sum_{w in Box cap {sum w = T}}  prod_i [w_i]_z   log-concave?

Equivalently, with l_i=L_i-1, h_i=U_i-1, D=T-m:
   k(a) = # { (x,u) in Z_{>=0}^m x Z_{>=0}^m : sum x = a, sum u = D-a,
                                               l_i <= x_i+u_i <= h_i }
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ycoord import centred_sum
from gen import is_pf2, conv_intervals
from itertools import product
from collections import Counter

def boxslice(L, U, T):
    return [w for w in product(*[range(L[i], U[i]+1) for i in range(len(L))]) if sum(w) == T]

def F(L, U, T):
    W = boxslice(L, U, T)
    if not W: return None, 0
    acc = Counter()
    for w in W:
        co = conv_intervals(list(w))
        for j, c in enumerate(co): acc[j] += c
    hi = max(acc)
    return [acc.get(i, 0) for i in range(hi+1)], len(W)

def twosided_count(l, h, D):
    """k(a) via the (x,u) model -- an INDEPENDENT instrument for the same number."""
    m = len(l); out = Counter()
    for v in product(*[range(l[i], h[i]+1) for i in range(m)]):
        if sum(v) != D: continue
        co = conv_intervals([vi+1 for vi in v])
        for j, c in enumerate(co): out[j] += c
    return out

if __name__ == '__main__':
    print("=== (Q): F(z) = sum over a BOX SLICE of prod_i [w_i]_z ===")
    st = Counter(); wit = []
    for m in range(2, 6):
        for Umax in range(1, 8):
            for L in product(range(1, 4), repeat=m):
                for U in product(*[range(L[i], Umax+1) for i in range(m)]):
                    for T in range(sum(L), sum(U)+1):
                        f, nW = F(list(L), list(U), T)
                        if f is None or nW < 1: continue
                        st['cases'] += 1
                        st['maxW'] = max(st.get('maxW', 0), nW)
                        st['maxdeg'] = max(st.get('maxdeg', 0), len(f)-1)
                        if is_pf2(f): st['PF2'] += 1
                        else:
                            st['FAIL'] += 1
                            if len(wit) < 5: wit.append((m, L, U, T, nW, f))
    print(dict(st))
    for w in wit: print("   FAIL", w)
