"""Faster slide search on the census pool (a+b<=9): for each reachable mosaic, apply every p-preserving 0-triangle slide
(tau x w + parallelograms on one side of tau+f  ->  tau x w' + the others) and validate the result independently:
half-edge matching, no overlap of the new tiles with any tile, labelled-intersection (is_mixed).  Exact arithmetic."""
import sys, pickle, itertools
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E, control_pool as CP
E.SIGN = 1; f0 = f1 = CP.MAPS['(i+j,i)']
pool = pickle.load(open('pool.pkl', 'rb')); st = Counter(); example = None
for (a, b), D in sorted(pool.items()):
    if a + b > 9: continue
    kinds = {vs: k for k, vs in E.make_tiles(a, b)}; E1 = E.unit_edges(b)
    allpts = set(v for vs in kinds for v in vs)
    for T in D:
        cells = set(frozenset(CP.conv(v, f0, f1) for v in vs) for _, vs in T); nfound = 0
        for c in cells:
            B0, B1 = E.parts(c)
            if len(B1) != 1 or len(B0) != 3: continue
            w = next(iter(B1))
            for f in E1:
                if w not in f: continue
                (w2,) = tuple(f - {w}); edges = [frozenset(e) for e in itertools.combinations(sorted(B0), 2)]
                par = {e: frozenset(x + y for x in e for y in f) for e in edges}
                present = [e for e in edges if par[e] in cells]
                if not present: continue
                new = [frozenset(x + w2 for x in B0)] + [par[e] for e in edges if e not in present]
                out = (cells - {c} - {par[e] for e in present}) | set(new)
                # independent validation
                poly = {x: E.ccw(x) for x in out}; he = Counter()
                for P in poly.values():
                    for j in range(len(P)): he[(P[j], P[(j+1) % len(P)])] += 1
                ok = all(n == 1 for n in he.values()) and all((v, u) in he or all(E.cross(u, v, x) >= 0 for x in allpts) for (u, v) in he)
                ok = ok and all(not E.overlap(poly[x], poly[y]) for x in new for y in out if y != x)
                ok = ok and E.is_mixed([(kinds[x], x) for x in out], range(len(out)))
                st[((a, b), ok)] += 1
                if ok:
                    nfound += 1
                    if example is None: example = ((a, b), sorted(map(sorted, cells)), sorted(map(sorted, out)))
        st[((a, b), 'mosaics with a valid p-preserving slide' if nfound else 'mosaics with none')] += 1
for k in sorted(st, key=str): print(k, st[k])
pickle.dump(example, open('slide_example.pkl', 'wb'))
