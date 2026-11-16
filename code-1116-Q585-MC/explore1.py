"""(q,p) coordinates: q[tau] = 1-part of the 0-triangle with 0-part tau; p[sig] = 0-part of 1-triangle with 1-part sig.
Tests: (1) (q,p) determines the tiling (injective over reachable states);
(2) per-step effect of type a / type b flips on (q,p);
(3) net effect of a single migration on (q,p)."""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1115c2-Q585-G1')
from rec import *
from collections import Counter, defaultdict
def qp(b):
    q, p = {}, {}
    for kind, vs in b.tiles.values():
        if kind != 'tri': continue
        p0 = frozenset(v[:2] for v in vs); p1 = frozenset(v[2:] for v in vs)
        if len(p1) == 1: q[p0] = next(iter(p1))
        else: p[p1] = next(iter(p0))
    return q, p
def qpkey(b):
    q, p = qp(b); return (frozenset(q.items()), frozenset(p.items()))
st = Counter(); seen = {}; n = 0
ex = []
for key, b0 in census():
    for b, DV in states(b0, cap=25):
        k = qpkey(b); full = (tuple(tiles_of(b)), frozenset(b.nest[0]['pos']), frozenset(b.nest[1]['pos']), frozenset(b.nest[2]['pos']))
        kk = (key, k)
        if kk in seen: st['qp collision'] += 1; st['collision same tiling'] += seen[kk][0] == full[0]; st['collision same full'] += seen[kk] == full
        else: seen[kk] = full
        st['states'] += 1
        # net effect of each single migration
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                q0, p0 = qp(b); bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); q1, p1 = qp(bx)
                dq = {tau: tuple(a-c for a, c in zip(q1[tau], q0[tau])) for tau in q0 if q1[tau] != q0[tau]}
                dp = {s: tuple(a-c for a, c in zip(p1[s], p0[s])) for s in p0 if p1[s] != p0[s]}
                st[(t, 'dq vectors')] = st.get((t,'dq vectors'), Counter()) if False else 0
                for v in dq.values(): st[(t, 'dq', v)] += 1
                for v in dp.values(): st[(t, 'dp', v)] += 1
                st[(t, 'migrations')] += 1
                if len(ex) < 3 and dq: ex.append((key, x, t, dq, dp))
    n += 1
    if n >= 40: break
for k in sorted(st, key=str): print(k, st[k])
for e in ex: print(e)
