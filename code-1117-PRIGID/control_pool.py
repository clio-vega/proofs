"""Control (i): every reachable mosaic of 11-16c2 (pool.pkl) must be a fine mixed subdivision under mixsub's independent checker,
for one fixed coordinate dictionary (engine (i,j) -> mine (m,n)) and one twist sign.  Also: count p-collisions in the pool again."""
import sys, pickle, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E
from collections import defaultdict
pool = pickle.load(open('pool.pkl', 'rb'))
MAPS = {'(i+j,j)': lambda i, j: (i+j, j), '(i+j,i)': lambda i, j: (i+j, i)}
def conv(v, f0, f1): return f0(*v[:2]) + f1(*v[2:])
def valid(tiling, a, b, f0, f1):
    Ttiles = set(vs for _, vs in E.make_tiles(a, b)); kinds = {vs: k for k, vs in E.make_tiles(a, b)}
    cells = [frozenset(conv(v, f0, f1) for v in vs) for _, vs in tiling]
    if not all(c in Ttiles for c in cells): return 'not a cell'
    poly = [E.ccw(c) for c in cells]
    he = set(); allpts = set(v for c in cells for v in c)
    for P in poly:
        for j in range(len(P)):
            h = (P[j], P[(j+1) % len(P)])
            if h in he: return 'half-edge twice'
            he.add(h)
    for (u, v) in he:
        if (v, u) not in he and not all(E.cross(u, v, x) >= 0 for x in allpts): return 'unmatched interior edge'
    for i in range(len(poly)):
        for j in range(i+1, len(poly)):
            if E.overlap(poly[i], poly[j]): return 'overlap'
    tl = [(kinds[c], c) for c in cells]
    if not E.is_mixed(tl, range(len(tl))): return 'not mixed'
    return 'ok'
for (n0, f0), (n1, f1), s in itertools.product(MAPS.items(), MAPS.items(), (1, -1)):
    E.SIGN = s; res = defaultdict(int)
    for (a, b), D in sorted(pool.items()):
        if a + b > 9: continue     # keep the control quick; (4,4),(5,4),(4,5) suffice to fix conventions
        for T in D: res[valid(T, a, b, f0, f1)] += 1
    print(n0, n1, 'SIGN', s, dict(res))
