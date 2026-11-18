"""Wake translation test. For each solo journey J (x in M; y in M) let Rm = tiles removed net (M minus Mx), Ad = tiles added net
(excluding rhombi).  Is there a matching Rm -> Ad by translations, and how many distinct translation vectors per journey?
Also: is each translation vector an edge vector of the migrating rhombus (at some time)?"""
import copy
from collections import Counter
from q597lib import *
def net(J, M):
    T = set((k, frozenset(vs)) for k, vs in M.values())
    T0 = set(T)
    for _, g, n in J: T = (T - g) | n
    rm = {t for t in T0 - T if t[0] != 'rh'}; ad = {t for t in T - T0 if t[0] != 'rh'}
    return rm, ad
def mv(t, v): return (t[0], frozenset(tuple(a+b for a, b in zip(P, v)) for P in t[1]))
st = Counter(); vecs = Counter()
for r in R:
    if r['kind'] != 'BC': continue
    for J in (r['Jx'], r['Jy']):
        rm, ad = net(J, r['M'])
        # candidate vectors: differences of vertices of rm tile and ad tile of same kind
        V = Counter()
        assign = {}
        for t in rm:
            p0 = min(t[1]); cands = set()
            for u in ad:
                if u[0] != t[0]: continue
                for q in u[1]:
                    v = tuple(b-a for a, b in zip(p0, q))
                    if mv(t, v) == u: cands.add(v)
            assign[t] = cands
        ok = all(len(c) >= 1 for c in assign.values())
        uniq = all(len(c) == 1 for c in assign.values())
        st[('every removed tile is a translate of an added tile', ok, 'unique', uniq, 'balanced', len(rm) == len(ad))] += 1
        if ok and uniq:
            vs = Counter(next(iter(c)) for c in assign.values())
            st[('#distinct vectors', len(vs))] += 1
            for v in vs: vecs[v] += 1
for k in sorted(st, key=str): print(st[k], k)
print(vecs.most_common(12))
