"""Type-a flip lemma check on every recorded step in align.pkl (<=5): for type a, removed tiles g and added tiles n:
is n = { t + v_t : t in g } with each t moved by a lattice vector, and is v_t the zonotope generator NOT an edge of t?
For type b: how many tiles of g have a translate in n?"""
from replay import *
from collections import Counter
from itertools import combinations
UNIT = {(1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1)}
def evecs(t):
    out = set()
    for P, Q in combinations(t[1], 2):
        d = tuple(b-a for a, b in zip(P, Q))
        if (d[2:] == (0, 0) and d[:2] in UNIT) or (d[:2] == (0, 0) and d[2:] in UNIT):
            out.add(d if d > tuple(-x for x in d) else tuple(-x for x in d))
    return out
def T(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
st = Counter(); seen = set()
for r in R:
    for J in (r['Jx'], r['Jy'], r['Jxy'], r['Jyx']):
        for s in J:
            if s in seen: continue
            seen.add(s); typ, g, n = s
            gens = set().union(*[evecs(t) for t in g])
            moved = 0; lacks = 0
            for t in g:
                ok = [v for v in (tuple(b-a for a, b in zip(min(t[1]), min(u[1]))) for u in n if u[0] == t[0]) if T(t, v) in n]
                if ok:
                    moved += 1
                    v = ok[0]; vv = v if v > tuple(-x for x in v) else tuple(-x for x in v)
                    lacks += (vv in gens and vv not in evecs(t))
            st[(typ, '#generators', len(gens), 'tiles translated', moved, 'by lacking generator', lacks)] += 1
for k in sorted(st, key=str): print(st[k], k)
