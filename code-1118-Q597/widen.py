"""PROVE 11-18 step 0. SCEL at a wider sampler cap, signatures FROZEN from 11-17 c2 (rolesigs.pkl) -> genuinely out of sample.
usage: python3 widen.py CAP.  Reports: B/C commute count, SCEL predicted-I == observed-I, on-I formula, |U| histogram, B/B control."""
import sys, copy, pickle
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *
from collections import Counter
CAP = int(sys.argv[1])
S = pickle.load(open('rolesigs.pkl', 'rb'))
Z = (0, 0); D = {0: (1, -1), 1: (1, 0)}; SV = {0: (0, -1), 1: (0, 1)}
def pq(b):
    q, p = qp(b); return {**{(0, t): v for t, v in q.items()}, **{(1, s): v for s, v in p.items()}}
def dl(A, B): return {k: tuple(x-y for x, y in zip(B[k], A[k])) for k in A if A[k] != B[k]}
def mv(b, x, T, DV):
    bb = copy.deepcopy(b); mig1(bb, 0, x, DV[T]); return bb
def pos(u):
    m, n = u
    ups = [(m, n), (m-1, n), (m, n-1)]; dns = [(m-1, n), (m, n-1), (m-1, n-1)]
    return [frozenset({(a, b), (a+1, b), (a, b+1)}) for a, b in ups] + [frozenset({(a+1, b), (a, b+1), (a+1, b+1)}) for a, b in dns]
def add(u, d, j=1): return (u[0]+j*d[0], u[1]+j*d[1])
def sigs(P, X, Y, side):
    tris = {k[1] for k in P if k[0] == side}; out = {}
    for u in {u for T in tris for u in T}:
        P6 = pos(u)
        if all(T in tris for T in P6): out[u] = tuple((X.get((side, T), Z), Y.get((side, T), Z)) for T in P6)
    return out
def predict(P, X, Y):
    Ip = set(); Us = []
    for side in (0, 1):
        sg = sigs(P, X, Y, side); d = D[side]
        for u, s in sg.items():
            if s in S[(side, 'single')]: Ip |= {(side, T) for T in pos(u)}; Us.append((side, 1))
            if s in S[(side, 'first')]:
                j = 1
                while sg.get(add(u, d, j)) in S[(side, 'inner')]: j += 1
                if sg.get(add(u, d, j)) in S[(side, 'last')]:
                    for i in range(j+1): Ip |= {(side, T) for T in pos(add(u, d, i))}
                    Us.append((side, j+1))
    return Ip, Us
st = Counter(); seenstates = 0; ones = []
for key, b0 in census():
    for b, DV in states(b0, cap=CAP):
        seenstates += 1
        cs = corners(set(b.nest[0]['pos']))
        for x in cs:
            for y in cs:
                if x == y: continue
                for kind, (tx, ty) in (('BC', 'BC'), ('BB', 'BB')):
                    if kind == 'BB' and not x < y: continue
                    b1 = mv(b, x, tx, DV); b2 = mv(b, y, ty, DV); b12 = mv(b1, y, ty, DV); b21 = mv(b2, x, tx, DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    com = tiles_of(b12) == tiles_of(b21)
                    X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                    K = set(X) | set(Y) | set(Yx) | set(Xy)
                    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                    Ip, Us = predict(P, X, Y)
                    if kind == 'BC':
                        st[('BC com', com)] += 1
                        st[('BC SCEL I', Ip == I)] += 1
                        ok = all(tuple(a+c for a, c in zip(X.get(k, Z), Yx.get(k, Z))) == SV[k[0]] and
                                 tuple(a+c for a, c in zip(Y.get(k, Z), Xy.get(k, Z))) == SV[k[0]] for k in I)
                        st[('BC on-I formula', ok)] += 1
                        st[('BC |U| per side', tuple(sorted(Us)))] += 1
                        if Ip != I or not ok: ones.append((key, x, y, sorted(I), sorted(Ip)))
                    else:
                        st[('BB com', com)] += 1
                        st[('BB pred nonempty', bool(Ip))] += 1
print('cap', CAP, 'states', seenstates)
for k in sorted(st, key=str): print(k, st[k])
print('failures', len(ones)); pickle.dump(ones, open(f'widen-fail-{CAP}.pkl', 'wb'))
