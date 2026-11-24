"""Lift structure of census boards M0.  pi0 = Lambda0 part (a0,b0), pi1 = Lambda1 part (a1,b1).
kind of triangle: L1-tri if all vertices share the Lambda0 part, L0-tri if all share the Lambda1 part.
(H1) pi1 injective on L1-triangles?  pi0 injective on L0-triangles?
(H2) the image set pi1(L1-tris) (resp. pi0(L0-tris)) is the same for all states of the size?
Also on boards during journeys (after every step of every J_o).   usage: python3 proj.py SIZE"""
import sys, os, pickle
from collections import Counter
SIZE = int(sys.argv[1])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
def kindtri(t):
    V = list(t[1])
    if len({P[:2] for P in V}) == 1: return 1
    if len({P[2:] for P in V}) == 1: return 0
    return None
def proj(t, k): return frozenset(P[2:] if k == 1 else P[:2] for P in t[1])
st = Counter(); images = {0: Counter(), 1: Counter()}; seen = set()
def check(board, lab):
    for k in (0, 1):
        tr = [t for t in board if t[0] == 'tri' and kindtri(t) == k]
        im = [proj(t, k) for t in tr]
        st[(lab, 'L%d-tri pi%d injective' % (k, k), len(im) == len(set(im)))] += 1
        if lab == 'M0': images[k][frozenset(im)] += 1
    st[(lab, 'triangles of mixed kind', sum(1 for t in board if t[0] == 'tri' and kindtri(t) is None))] += 1
for P in PAIRS:
    key = frozenset(P['M0'])
    if key in seen: continue
    seen.add(key); M0 = set(P['M0']); check(M0, 'M0')
    for (tag, A, Bq, Jo) in P['roles']:
        bd = set(M0)
        for s in Jo:
            bd = (bd - s[1]) | s[2]; check(bd, 'during journey')
for k in sorted(st, key=str): print(st[k], k)
for k in (0, 1): print('distinct pi%d image sets of L%d-triangles over %d states:' % (k, k, len(seen)), len(images[k]), 'sizes', sorted({len(s) for s in images[k]}))
