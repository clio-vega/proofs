"""Bridge: every distinct census M0 equals (up to a Z^4 translation) build_mosaic(n, lab) for the KT puzzle lab read off M0
(labels: a puzzle edge psi(P)psi(Q) is 0 if Q-P is a Lambda0 unit step, 1 if a Lambda1 step, 10 = the short diagonal of a
square), and lab uses only Purbhoo-chirality pieces (mosaic.OK).  So Lemma D's hypothesis covers the census.
usage: python3 bridge.py SIZE"""
import sys, pickle
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
SIZE = int(sys.argv[1])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
U = {(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)}
def psi(P): return (P[0] + P[2], P[1] + P[3])
st = Counter(); seen = set()
for P in PAIRS:
    key = frozenset(P['M0'])
    if key in seen: continue
    seen.add(key); M0 = P['M0']
    allp = {psi(Q) for t in M0 if t[0] != 'rh' for Q in t[1]}
    a0 = min(p[0] for p in allp); b0 = min(p[1] for p in allp); n = max(p[0] + p[1] for p in allp) - a0 - b0
    sh = lambda p: (p[0] - a0, p[1] - b0)
    lab = {}; okl = True
    for t in M0:
        if t[0] == 'rh': continue
        V = list(t[1])
        for x in V:
            for y in V:
                D = tuple(b - a for a, b in zip(x, y))
                if (D[2:] == (0, 0) and D[:2] in U): l = 0
                elif (D[:2] == (0, 0) and D[2:] in U): l = 1
                else: continue
                e = MZ.E(sh(psi(x)), sh(psi(y)))
                if lab.get(e, l) != l: okl = False
                lab[e] = l
        if t[0] == 'sq':
            Vs = [sh(psi(x)) for x in V]
            for x in Vs:
                for y in Vs:
                    if x < y and (y[0] - x[0], y[1] - x[1]) in U and MZ.E(x, y) not in lab:
                        # the short diagonal: the unit pair that is not a side
                        lab[MZ.E(x, y)] = 2
    # every cell of the grid: ccw labels in MZ.OK
    cells_ok = all(tuple(lab.get(MZ.E(c[i], c[(i + 1) % 3])) for i in range(3)) in MZ.OK for c in MZ.tri_cells(n))
    st[('labels consistent', okl)] += 1; st[('every cell a Purbhoo-chirality KT piece', cells_ok)] += 1
    if not cells_ok: continue
    M = MZ.build_mosaic(n, lab)
    built = {(k, frozenset(v)) for k, v in M['tiles']} | {(k, frozenset(v)) for nest in M['nests'] for k, v in nest}
    # translate: match a vertex over psi = (0,0)
    Pm = [Q for t in M0 for Q in t[1] if sh(psi(Q)) == (0, 0) and t[0] != 'rh'][0]
    T0 = M['L'][(0, 0)]; d = tuple(a - b for a, b in zip(Pm, T0))
    built = {(k, frozenset(tuple(a + b for a, b in zip(Q, d)) for Q in v)) for k, v in built}
    st[('M0 == build_mosaic(puzzle) + translation', built == set(M0))] += 1
print('SIZE', SIZE, 'distinct M0', len(seen))
for k in sorted(st, key=str): print(st[k], k)
