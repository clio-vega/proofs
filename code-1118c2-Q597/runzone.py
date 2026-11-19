"""Run-zone lemma check: in every journey (align.pkl, <=5), split into maximal runs of type-a steps.
Within a run the rhombus shape P(e0,e1) is fixed (check).  For each displaced non-rh tile t (t in g, t+v in n):
classify by which rhombus edge it contains (e1-zone: has an edge parallel to e1; e0-zone: parallel to e0).
Lemma: within a run, all e1-zone displacements are ONE vector (= +-e0) and all e0-zone displacements ONE vector (= +-e1)."""
from replay import *
from collections import Counter
from itertools import combinations
UNIT = {(1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1)}
def ev(t):
    out = set()
    for P, Q in combinations(t[1], 2):
        d = tuple(b-a for a, b in zip(P, Q))
        if (d[2:] == (0, 0) and d[:2] in UNIT) or (d[:2] == (0, 0) and d[2:] in UNIT):
            out.add(max(d, tuple(-x for x in d)))
    return out
def Tt(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
st = Counter(); seen = set()
for r in R:
    for J in (r['Jx'], r['Jy'], r['Jxy'], r['Jyx']):
        key = tuple(J)
        if key in seen: continue
        seen.add(key)
        runs = []; cur = []
        for s in J:
            if s[0] == 'a': cur.append(s)
            else:
                if cur: runs.append(cur)
                cur = []
        if cur: runs.append(cur)
        for run in runs:
            shapes = {frozenset(ev(rh(s[1]) and ('rh', rh(s[1])))) for s in run}
            st[('rhombus edge-set fixed in run', len(shapes) == 1)] += 1
            e = ev(('rh', rh(run[0][1]))); e0 = next(v for v in e if v[2:] == (0, 0)); e1 = next(v for v in e if v[:2] == (0, 0))
            V0, V1 = set(), set()
            for s in run:
                for t in s[1]:
                    if t[0] == 'rh': continue
                    vs = [v for v in {tuple(b-a for a, b in zip(min(t[1]), min(u[1]))) for u in s[2] if u[0] == t[0]} if Tt(t, v) in s[2]]
                    if not vs: st['untranslated tile in type a'] += 1; continue
                    v = vs[0]; E = ev(t)
                    if e1 in E and e0 not in E: V1.add(v)
                    elif e0 in E and e1 not in E: V0.add(v)
                    else: st['tile has both/neither rhombus edge'] += 1
            st[('run len', min(len(run), 6), 'e1-zone vectors', len(V1), 'e0-zone vectors', len(V0))] += 1
            if len(V1) == 1: st[('e1-zone vector = +-e0', next(iter(V1)) in (e0, tuple(-x for x in e0)))] += 1
            if len(V0) == 1: st[('e0-zone vector = +-e1', next(iter(V0)) in (e1, tuple(-x for x in e1)))] += 1
for k in sorted(st, key=str): print(st[k], k)
