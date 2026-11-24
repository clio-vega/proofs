"""Fiber lemma (computed): for each census state M0 and each vertex x of Delta1 (Lambda1 projection): F(x) = the
Lambda0-triangles of M0 whose Lambda1 part is x; k_i = number of squares of M0 whose Lambda1 edge projects to the edge
(x, x+d_i), d_i the 6 Lambda1 units in ccw order.  Check: (1) the stacks close: sum k_i n_i = 0 with n_i the 6 Lambda0
units in the matching cyclic order (both orientations tried; report which works); (2) |F(x)| = (k0+k1+k2)^2 - k0^2 - k2^2
- k4^2 (unit triangles of the hexagon with those sides); (3) rhombi over edges (nests) counted separately.
Same with the roles of the lattices swapped.   usage: python3 fiber.py SIZE"""
import sys, pickle
from collections import Counter
SIZE = int(sys.argv[1])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
U = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
st = Counter(); seen = set()
def run(M0, k):
    fx = (lambda P: P[:2]) if k == 1 else (lambda P: P[2:])     # coordinate held fixed on a k-triangle
    fr = (lambda P: P[2:]) if k == 1 else (lambda P: P[:2])
    F = Counter(); stack = Counter(); rhs = Counter(); verts = set()
    for t in M0:
        V = list(t[1])
        if t[0] == 'tri':
            if len({fx(P) for P in V}) == 1: verts |= {fr(P) for P in V}
            elif len({fr(P) for P in V}) == 1: F[fr(V[0])] += 1
        else:
            E = frozenset(fr(P) for P in V)
            if len(E) == 2: (stack if t[0] == 'sq' else rhs)[E] += 1
    for x in verts:
        ks = [stack[frozenset((x, (x[0] + d[0], x[1] + d[1])))] for d in U]
        nr = sum(rhs[frozenset((x, (x[0] + d[0], x[1] + d[1])))] for d in U)
        ok = []
        for sgn in (1, -1):
            for rot in range(6):
                s = [0, 0]
                for i, kk in enumerate(ks):
                    n = U[(sgn * i + rot) % 6]; s[0] += kk * n[0]; s[1] += kk * n[1]
                ok.append(s == [0, 0])
        closes = all(ok)  # closure for every assignment iff ...; report any/all separately
        anyc = any(ok)
        area = (ks[0] + ks[1] + ks[2]) ** 2 - ks[0] ** 2 - ks[2] ** 2 - ks[4] ** 2
        st[('lattice-fixed', k, 'rhombi at x', nr > 0, 'closure some orientation', anyc, 'area formula', area == F[x])] += 1
for P in PAIRS:
    key = frozenset(P['M0'])
    if key in seen: continue
    seen.add(key)
    for k in (0, 1): run(P['M0'], k)
for kk in sorted(st, key=str): print(st[kk], kk)
