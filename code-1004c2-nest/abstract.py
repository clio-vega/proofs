"""SUFFICIENCY BEFORE SURPLUS.

Does (H3) suffice?  (H3)'s hypotheses are satisfied by EVERY triple (g, B, sigma):
Theorem M says every sublevel set {y in B : sum y = sigma, sum (-y_i)_+ <= j} is
M-convex, for every box B, and the width map is eq:w, which needs only g >= 0.
So the abstract class (H3) carves out is

    C = { (g, B, sigma) : g in Z_{>=0}^m,  B = prod [P_i,Q_i],
                          Y = B cap {sum y = sigma} nonempty,
                          w_i(y) >= 1 for every y in Y and every i }

(the last clause because on a real slice every summand is nonzero -- the effective
slice IS the box, A-M-profile-box-slice).  Every real slice lies in C.

QUESTION: is sum_{y in Y} Trap_{w(y)} always PF_2 on C?
If NO, (H3) is INSUFFICIENT and the route is dead.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ycoord import wmap, centred_sum, defect
from gen import is_pf2
from itertools import product
from collections import Counter

def slice_pts(P, Q, sig):
    m = len(P)
    rngs = [range(P[i], Q[i]+1) for i in range(m)]
    return [y for y in product(*rngs) if sum(y) == sig]

def test(g, P, Q, sig):
    """-> (status, ws, sum) ; status in {'empty','zerowidth','ok','FAIL'}"""
    Y = slice_pts(P, Q, sig)
    if not Y: return ('empty', None, None)
    ws = [wmap(y, g) for y in Y]
    if any(min(w) <= 0 for w in ws): return ('zerowidth', None, None)
    cs = centred_sum(ws)
    if cs is None: return ('parity', None, None)
    return ('ok' if is_pf2(cs) else 'FAIL', ws, cs)

def sweep(m, grange, brange, label, maxfail=6):
    st = Counter(); fails = []
    gs = [gg for gg in product(range(grange+1), repeat=m)]
    for g in gs:
        for P in product(range(-brange, brange+1), repeat=m):
            for Q in product(*[range(P[i], brange+1) for i in range(m)]):
                for sig in range(sum(P), sum(Q)+1):
                    s, ws, cs = test(g, P, Q, sig)
                    st[s] += 1
                    if s == 'FAIL':
                        st['sizes'] = max(st.get('sizes',0), len(ws))
                        if len(fails) < maxfail:
                            fails.append((g, P, Q, sig, ws, cs))
    print(f"[{label}] m={m} g<={grange} box in [-{brange},{brange}]: {dict(st)}")
    for f in fails:
        print(f"   FAIL g={f[0]} P={f[1]} Q={f[2]} sigma={f[3]}")
        print(f"        |Y|={len(f[4])} widths={sorted(f[4])}")
        print(f"        half-widths={sorted((sum(w)-m)//2 for w in f[4])}")
        print(f"        sum={f[5]}")
    return st, fails

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'm2'
    if which == 'm2':
        sweep(2, 8, 8, 'abstract m=2')
    elif which == 'm3':
        sweep(3, 4, 4, 'abstract m=3')
    elif which == 'm4':
        sweep(4, 3, 3, 'abstract m=4')
