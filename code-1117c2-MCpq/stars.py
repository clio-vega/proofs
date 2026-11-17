"""Test: interaction set I (cells whose composite differs from no-interaction) is the union of stars of a set U of interior
vertices, all cells of I carrying one value v; U consecutive along a line.  Also record where X, Y sit on I."""
import pickle, itertools
from collections import Counter
BCo, BBn = pickle.load(open('sets.pkl', 'rb'))
def vec(D): return {k: (n[0]-o[0], n[1]-o[1]) for k, (o, n) in D.items()}
Z = (0, 0)
def star(u):
    m, n = u
    ups = [(m, n), (m-1, n), (m, n-1)]; dns = [(m-1, n), (m, n-1), (m-1, n-1)]
    return {frozenset({(a, b), (a+1, b), (a, b+1)}) for a, b in ups} | {frozenset({(a+1, b), (a, b+1), (a+1, b+1)}) for a, b in dns}
# check star definition: each triangle contains u and the 6 are distinct
assert all(len(star((5, 3))) == 6 and all((5, 3) in T for T in star((5, 3))) for _ in [0])
st = Counter()
for r in BCo:
    X, Y, Yx, Xy = map(vec, (r['X'], r['Y'], r['Yx'], r['Xy']))
    K = set(X) | set(Y) | set(Yx) | set(Xy)
    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
    if not I: continue
    sides = {k[0] for k in I}; vals = {r['P'][k] for k in I}
    st[('one side', len(sides) == 1, 'one value', len(vals) == 1)] += 1
    side = next(iter(sides)); tris = {k[1] for k in I}
    # vertices whose full star lies in I
    cand = {u for T in tris for u in T if star(u) <= tris}
    cover = set().union(*[star(u) for u in cand]) if cand else set()
    st[('I == union of full stars', cover == tris)] += 1
    U = sorted(cand)
    if len(U) >= 2:
        d = {(b[0]-a[0], b[1]-a[1]) for a, b in zip(U, U[1:])}
        st[('|U|', len(U), 'steps', tuple(sorted(d)), 'side', side)] += 1
    else: st[('|U|', len(U), 'side', side)] += 1
    # is every cell with value v around U ... the whole star region assigned v also BEFORE; and are there other cells with value v?
    st[('values', len(vals), 'U', len(U))] += 1
    v = next(iter(vals))
    samev = {k[1] for k in r['P'] if k[0] == side and r['P'][k] == v}
    st[('cells with value v == I', samev == tris, 'I subset', tris <= samev)] += 1
for k, c in sorted(st.items(), key=str): print(c, k)
