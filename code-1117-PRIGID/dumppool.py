"""Dump the reachable mosaic pool of 11-16c2 (ponly_hex.py logic, unchanged) to pool.pkl for control (i)."""
import sys, pickle, copy
sys.path.insert(0, '/home/clio/projects/proofs/code-1116-Q585-MC')
from qp import *
from collections import defaultdict
pool = defaultdict(dict)
for key, b0 in census():
    lam, mu, nu, bb = key; tn = trim(nu); nt = (max(tn[0], 1)+1, len(tn)+1)
    A = b0.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    seen = set(); stack = [b0]
    while stack:
        b = stack.pop(); T = tuple(tiles_of(b))
        if T in seen: continue
        seen.add(T)
        q, p = qp(b)
        pool[nt][T] = (dict(q), dict(p))
        for x in corners(set(b.nest[0]['pos'])):
            for t in 'BC':
                bx = copy.deepcopy(b); mig1(bx, 0, x, DV[t]); stack.append(bx)
pickle.dump(dict(pool), open('pool.pkl', 'wb'))
for nt, D in pool.items():
    T = next(iter(D)); V = set(v for k, vs in T for v in vs)
    A0 = sorted(set(v[:2] for v in V)); A1 = sorted(set(v[2:] for v in V))
    kinds = defaultdict(int)
    for k, vs in T: kinds[(k, len(vs))] += 1
    print(nt, len(D), 'A0', len(A0), A0[:12], 'A1', len(A1), A1[:12], dict(kinds))
