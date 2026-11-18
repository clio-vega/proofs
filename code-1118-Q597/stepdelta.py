"""Per-step pure-cell bookkeeping. For a step (typ, gone, new): pure cells of gone tri tiles vs new tri tiles.
Q1: does every step preserve the SET of pure-cell keys (side, P_side-triangle)?  Q2: per-step delta delta(s): key -> value change.
Q3: is a journey's net pure delta the sum of its step deltas?  (checked on all journeys of all records)"""
import pickle
from collections import Counter
exec(open('analyse_align.py').read().split('T = Counter()')[0])
def pure(S):
    out = {}
    for k, vs in S:
        if k != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: out[(0, p0)] = next(iter(p1))
        else: out[(1, p1)] = next(iter(p0))
    return out
def sdelta(step):
    _, g, n = step; a, b = pure(g), pure(n)
    if set(a) != set(b): return None
    return {k: tuple(y-x for x, y in zip(a[k], b[k])) for k in a if a[k] != b[k]}
def addd(D, E):
    out = dict(D)
    for k, v in E.items():
        w = tuple(a+b for a, b in zip(out.get(k, (0, 0)), v)); out[k] = w
        if w == (0, 0): del out[k]
    return out
st = Counter(); kinds = Counter()
for r in R:
    for J, net in ((r['Jy'], r['Y']), (r['Jx'], r['X'])):
        tot = {}; ok = True
        for s in J:
            d = sdelta(s)
            if d is None: ok = False; break
            kinds[(s[0], len(pure(s[1])), tuple(sorted(k[0] for k in pure(s[1]))), len(d))] += 1
            tot = addd(tot, d)
        st[('Q1 every step keyset-preserving', ok)] += 1
        if ok: st[('Q3 net == sum of step deltas', tot == net)] += 1
for k in sorted(st, key=str): print(st[k], k)
print('step kinds (typ, #tri, sides of tris, #cells changed):')
for k in sorted(kinds, key=str): print(kinds[k], k)
