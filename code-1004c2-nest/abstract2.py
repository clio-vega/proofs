"""SUFFICIENCY of (H3) in the abstract class, efficient parametrisation.

Reparametrise: s_i(y) := (y_i)_+ + (-y_{i+1})_+ >= 0.  Then eq:w reads
      w_i(y) = (g_i + 1) - s_i(y),
so the width multiset of a slice is determined by the s-vectors of Y together
with g, and w_i >= 1 on Y iff g_i >= max_{y in Y} s_i(y).  Hence the abstract
class is enumerated by (box, sigma) -> Y, then g_i >= maxs_i, free above that.

Theorem M holds for every box, so EVERY such triple satisfies (H3)'s hypotheses.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ycoord import centred_sum
from gen import is_pf2
from itertools import product
from collections import Counter

pos = lambda t: t if t > 0 else 0

def svec(y):
    m = len(y)
    return tuple(pos(y[i]) + pos(-y[(i+1) % m]) for i in range(m))

def check(Y, g):
    m = len(g)
    ws = [tuple(g[i] + 1 - s[i] for i in range(m)) for s in (svec(y) for y in Y)]
    if any(min(w) <= 0 for w in ws): return None
    cs = centred_sum(ws)
    return (is_pf2(cs), ws, cs)

def random_sweep(m, N, R, slack, seed=0):
    rng = random.Random(seed)
    st = Counter(); fails = []
    for _ in range(N):
        P = [rng.randint(-R, R) for _ in range(m)]
        Q = [p + rng.randint(0, R) for p in P]
        lo, hi = sum(P), sum(Q)
        sig = rng.randint(lo, hi)
        Y = [y for y in product(*[range(P[i], Q[i]+1) for i in range(m)]) if sum(y) == sig]
        if len(Y) < 2: st['tiny'] += 1; continue
        S = [svec(y) for y in Y]
        maxs = [max(s[i] for s in S) for i in range(m)]
        g = tuple(maxs[i] + rng.randint(0, slack) for i in range(m))
        r = check(Y, g)
        if r is None: st['zerowidth'] += 1; continue
        ok, ws, cs = r
        st['ok' if ok else 'FAIL'] += 1
        st['G_max'] = max(st.get('G_max', 0), sum(g) - m if False else sum(g))
        st['Y_max'] = max(st.get('Y_max', 0), len(Y))
        if not ok and len(fails) < 8:
            fails.append((g, tuple(P), tuple(Q), sig, ws, cs))
    print(f"[random m={m} N={N} R={R} slack={slack}] {dict(st)}")
    for f in fails:
        print(f"   FAIL g={f[0]} P={f[1]} Q={f[2]} sigma={f[3]} |Y|={len(f[4])}")
        print(f"        widths={sorted(f[4])}  halfw={sorted((sum(w)-m)//2 for w in f[4])}")
        print(f"        sum={f[5]}")
    return fails

if __name__ == '__main__':
    allf = []
    for m in (3, 4, 5):
        for R, slack in ((3, 2), (5, 3), (7, 4)):
            allf += random_sweep(m, 20000, R, slack, seed=m*100+R)
    print("total fails:", len(allf))
