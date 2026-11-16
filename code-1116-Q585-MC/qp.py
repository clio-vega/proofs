"""(q,p) coordinates of a mosaic state and the potential Phi.  Exact integers only.
q[tau] = 1-part (vertex of Lambda_1) of the 0-triangle whose 0-part is tau; p[sig] likewise.
Phi: on 0-triangles (-x2, x3) and on 1-triangles (y2, -y3), barycentric x=(a,b,-a-b)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1115c2-Q585-G1')
from rec import *
def qp(b):
    q, p = {}, {}
    for kind, vs in b.tiles.values():
        if kind != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: q[p0] = next(iter(p1))
        else: p[p1] = next(iter(p0))
    return q, p
def phi(b):
    q, p = qp(b); out = {}
    for t, (a, c) in q.items(): out[('q2', t)] = -c; out[('q3', t)] = -a-c
    for s, (a, c) in p.items(): out[('p2', s)] = c; out[('p3', s)] = a+c
    return out
def leq(P1, P2): return all(P1[k] <= P2[k] for k in P1)
def nestdata(b): return tuple(frozenset(b.nest[i]['pos']) for i in range(3))
def fullkey(b): return (nestdata(b), tuple(tiles_of(b)))
