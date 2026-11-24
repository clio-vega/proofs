"""(D) lesion: a cochain that does not close on one tile has no lift.  For every KT puzzle of size <= NMAX and every
interior 0-edge, add (0,0,1,0) to that edge's vector; the path-integration lift must report an inconsistency.
Control: unperturbed lift succeeds.   usage: python3 perturb.py NMAX"""
import sys
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
def lift(n, lab, bump=None):
    adj = {}
    for e, l in lab.items():
        if l == 2: continue
        p, q = tuple(e); adj.setdefault(p, []).append((q, l, e)); adj.setdefault(q, []).append((p, l, e))
    L = {(0, 0): (0, 0, 0, 0)}; st = [(0, 0)]
    while st:
        p = st.pop()
        for q, l, e in adj[p]:
            d = (q[0] - p[0], q[1] - p[1]); v = (d[0], d[1], 0, 0) if l == 0 else (0, 0, d[0], d[1])
            if e == bump:
                s = 1 if p < q else -1; v = (v[0], v[1], v[2] + s, v[3])
            Q = tuple(a + b for a, b in zip(L[p], v))
            if q in L:
                if L[q] != Q: return None
            else: L[q] = Q; st.append(q)
    return L
st = Counter()
for n in range(2, int(sys.argv[1]) + 1):
    bnd = set(MZ.E(*e) for s in MZ.side_edges(n) for e in s)
    for lab in MZ.puzzles(n, {}):
        st[('control lift ok', lift(n, lab) is not None)] += 1
        for e, l in lab.items():
            if l == 0 and e not in bnd: st[('perturbed lift exists', lift(n, lab, e) is not None)] += 1
for k in sorted(st, key=str): print(st[k], k)
