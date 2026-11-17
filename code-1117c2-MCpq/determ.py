"""Is the interaction set I a function of SOLO data (X, Y on M)?  For every star(u) present in the cell set, the solo signature
sig(u) = tuple over the 6 star positions of (X, Y) (vectors relative), and the label 'u in U' (star(u) <= I).
Test: label is a function of (side, sig).  Run on all BC records (overlapping + disjoint) and on BB (control)."""
import pickle
from collections import Counter, defaultdict
exec(open('stars.py').read().split('st = Counter()')[0].replace("BCo, BBn = pickle.load(open('sets.pkl', 'rb'))", ""))
recs = pickle.load(open('recs.pkl', 'rb'))
def pos(u):
    m, n = u
    ups = [(m, n), (m-1, n), (m, n-1)]; dns = [(m-1, n), (m, n-1), (m-1, n-1)]
    return [frozenset({(a, b), (a+1, b), (a, b+1)}) for a, b in ups] + [frozenset({(a+1, b), (a, b+1), (a+1, b+1)}) for a, b in dns]
for kind in ('BC', 'BB'):
    lab = defaultdict(Counter); nI = Counter(); exact = Counter()
    for r in recs:
        if r['kind'] != kind: continue
        X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
        K = set(X) | set(Y) | set(Yx) | set(Xy)
        I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
        for side in (0, 1):
            tris = {k[1] for k in r['P'] if k[0] == side}
            verts = {u for T in tris for u in T}
            for u in verts:
                P6 = pos(u)
                if not all(T in tris for T in P6): continue
                sig = tuple((X.get((side, T), Z), Y.get((side, T), Z)) for T in P6)
                inU = all((side, T) in I for T in P6)
                lab[(side, sig)][inU] += 1
        # law check: on I, Yx = s - X, Xy = s - Y
        if I:
            s = {0: (0, -1), 1: (0, 1)}
            ok = all(tuple(a+b for a, b in zip(X.get(k, Z), Yx.get(k, Z))) == s[k[0]] and
                     tuple(a+b for a, b in zip(Y.get(k, Z), Xy.get(k, Z))) == s[k[0]] for k in I)
            exact[('law on I holds', ok)] += 1
    amb = {k: v for k, v in lab.items() if len(v) > 1}
    pos_sigs = {k: v for k, v in lab.items() if v.get(True)}
    print(kind, 'distinct star signatures', len(lab), 'ambiguous', len(amb), 'signatures with u in U', len(pos_sigs), dict(exact))
    for k, v in pos_sigs.items(): print('   U-signature side', k[0], dict(v), [tuple(map(tuple, t)) for t in k[1]])
    for k, v in list(amb.items())[:5]: print('   AMBIG', k, dict(v))
