import sys, pickle
from collections import defaultdict
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E, windcheck as W
a, b, M = pickle.load(open('all_2_2.pkl', 'rb')); M = list(M)
rec = defaultdict(list)
for T in M:
    pure_other = [vs for k, vs in T if k == 't1']
    pure_own = {}
    for k, vs in T:
        if k == 't0': B0, B1 = E.parts(vs); pure_own[B0] = (vs, next(iter(B1)))
    rib = W.midlines(T, 0)
    for e, fs in rib.items():
        adj = defaultdict(list)
        for f in fs:
            u, v = tuple(f); adj[u].append(v); adj[v].append(u)
        endsW = [x for x in adj if len(adj[x]) == 1]
        if len(endsW) != 2: continue
        path = [endsW[0]]; prev = None
        while True:
            nx = [y for y in adj[path[-1]] if y != prev]
            if not nx: break
            prev = path[-1]; path.append(nx[0])
        pts = [W.ptsum([x + w for x in e], 3) for w in path]
        endtiles = frozenset(vs for own, (vs, w) in pure_own.items() if e <= own and w in (path[0], path[-1]))
        for O in pure_other:
            o = W.ptsum(list(O), 2)
            ind = W.winding(pts, o)
            rec[(e, frozenset((path[0], path[-1])), endtiles, O)].append((ind, path, T))
for key, L in rec.items():
    if len(set(x[0] for x in L)) > 1:
        e, ends, endtiles, O = key
        print('e', sorted(e), 'path ends', sorted(ends), 'end tiles', [sorted(t) for t in endtiles], 'O', sorted(O))
        seen = set()
        for ind, path, T in L:
            if (ind, tuple(path)) in seen: continue
            seen.add((ind, tuple(path))); print('   ind', ind, 'path', path)
        break
