"""Boundary row of Theorem E1Z.  For every census state (nest5/nest6.pkl, distinct M0): nest index 2 = C (check frame),
1 = B.  Would-be cells outside the frame: C: (-1, j); B: (i, -1), for 0 <= i, j <= 12.  Does the rhombus of such a cell
meet the interior of ANY M0 tile (i.e. lie inside the region)?  Control: cells (0, j) / (i, 0) inside the nest DO coincide
with M0 nest rhombi.   usage: python3 boundrow.py"""
import sys, pickle
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import overlap
def add(P, Q, k=1): return tuple(a + k * b for a, b in zip(P, Q))
st = Counter(); seen = set()
for SIZE in (5, 6):
    for P in pickle.load(open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d.pkl' % SIZE, 'rb')):
        key = frozenset(P['M0'])
        if key in seen: continue
        seen.add(key); M0 = list(P['M0'])
        for idx, name in ((1, 'B'), (2, 'C')):
            N = P['nests'][idx]; X, u, w = N['X'], N['u'], N['w']
            st[(name, 'frame u,w', u, w)] += 1
            for t in range(0, 13):
                i, j = (-1, t) if name == 'C' else (t, -1)
                q = add(add(X, u, i), w, j)
                rh = ('rh', frozenset({q, add(q, u), add(q, w), add(add(q, u), w)}))
                inside = any(overlap(rh, m) for m in M0)
                st[(SIZE, name, 'outside-frame cell meets region', inside)] += 1
            for (i, j) in N['pos']:
                q = add(add(X, u, i), w, j); rh = ('rh', frozenset({q, add(q, u), add(q, w), add(add(q, u), w)}))
                st[(SIZE, name, 'CONTROL nest cell rhombus is an M0 tile', rh in P['M0'])] += 1
for k in sorted(st, key=str): print(st[k], k)
