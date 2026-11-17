"""Independent test of Claim W (the lemma the proof rests on), exact arithmetic in Z[sqrt3].
For a 0-edge e, the midline alpha_S of R_e runs between fixed points a1,a2 (midpoints of the e-edges of the end tiles).
For a 1-triangle tile O, ind_S(e,O) := winding number of alpha_S . rho^{-1} around the centroid of O, rho = straight segment a2->a1.
Claim W says: among all S that contain the same two end tiles of R_e and the tile O, ind_S(e,O) is constant.
(Also the symmetric statement for 1-edges f and 0-triangle tiles.)  Coordinates scaled by 6 to stay integral."""
import sys, pickle
from collections import defaultdict, Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1117-PRIGID')
import mixsub as E
E.SIGN = 1
def P(v): return E.planar(v)
def scal(k, X): return (k*X[0], k*X[1])
def ptsum(vs, k):  # k * (sum of planar(v)) as a point in Z[sqrt3]^2
    X, Y = (0, 0), (0, 0)
    for v in vs: Xv, Yv = P(v); X = E.add(X, Xv); Y = E.add(Y, Yv)
    return (scal(k, X), scal(k, Y))
def orient(p, q, r):
    ax, ay = E.sub(q[0], p[0]), E.sub(q[1], p[1]); bx, by = E.sub(r[0], p[0]), E.sub(r[1], p[1])
    return E.sgn(E.sub(E.mul(ax, by), E.mul(ay, bx)))
def winding(poly, o):
    w = 0
    for i in range(len(poly)):
        p, q = poly[i], poly[(i+1) % len(poly)]
        if orient(p, q, o) == 0 and min(E.sgn(E.sub(p[1], o[1])), E.sgn(E.sub(q[1], o[1]))) <= 0 <= max(E.sgn(E.sub(p[1], o[1])), E.sgn(E.sub(q[1], o[1]))):
            # o collinear with an edge and within its y-range: check it is not ON the edge
            pass
        py, qy = E.sgn(E.sub(p[1], o[1])), E.sgn(E.sub(q[1], o[1]))
        if py <= 0 < qy and orient(p, q, o) > 0: w += 1
        elif qy <= 0 < py and orient(p, q, o) < 0: w -= 1
    return w
def midlines(T, side):
    """side 0: for each 0-edge e, ordered midpoints (scaled x3) of the e-edges along R_e, plus end tiles.  side 1 symmetric."""
    rib = defaultdict(list); ends = defaultdict(list)
    for k, vs in T:
        B0, B1 = E.parts(vs)
        if k == 'par':
            e, f = (B0, B1) if side == 0 else (B1, B0)
            rib[e].append(f)
    return rib
def run(fn, side, maxS=None):
    a, b, M = pickle.load(open(fn, 'rb')); M = list(M)[:maxS] if maxS else list(M)
    groups = defaultdict(set); routes = defaultdict(set); bad = 0; tested = 0
    for T in M:
        tilesT = set(vs for _, vs in T)
        pure_other = [vs for k, vs in T if k == ('t1' if side == 0 else 't0')]
        pure_own = {}
        for k, vs in T:
            if k == ('t0' if side == 0 else 't1'):
                B0, B1 = E.parts(vs); own, oth = (B0, B1) if side == 0 else (B1, B0); pure_own[own] = (vs, next(iter(oth)))
        rib = midlines(T, side)
        for e, fs in rib.items():
            # order the ribbon: walk the path in the other lattice
            adj = defaultdict(list)
            for f in fs:
                u, v = tuple(f); adj[u].append(v); adj[v].append(u)
            endsW = [x for x in adj if len(adj[x]) == 1]
            if len(endsW) != 2: continue
            path = [min(endsW)]; prev = None   # canonical orientation (fix 11-17: arbitrary start negated the index)
            while True:
                nx = [y for y in adj[path[-1]] if y != prev]
                if not nx: break
                prev = path[-1]; path.append(nx[0])
            def lab(x, w): return (x + w) if side == 0 else (w + x)
            pts = [ptsum([lab(x, w) for x in e], 3) for w in path]   # 3*(sum of 2 endpoints) = 6*midpoint
            endtiles = frozenset(vs for own, (vs, w) in pure_own.items() if e <= own and w in (path[0], path[-1]))
            for O in pure_other:
                o = ptsum(list(O), 2)  # 2*(sum of 3 vertices) = 6*centroid
                if any(orient(pts[i], pts[i+1], o) == 0 for i in range(len(pts)-1)) or orient(pts[-1], pts[0], o) == 0: continue
                ind = winding(pts, o)   # closed by the straight segment back
                groups[(e, frozenset(path[::len(path)-1]), endtiles, O)].add(ind); routes[(e, frozenset(path[::len(path)-1]), endtiles, O)].add(tuple(path)); tested += 1
    viol = sum(len(s) > 1 for s in groups.values())
    nontriv = sum(1 for s in groups.values() if s != {0})
    multi = sum(1 for k in groups if len(routes[k]) > 1); multi_nz = sum(1 for k in groups if len(routes[k]) > 1 and groups[k] != {0})
    print('   groups with >=2 distinct routes', multi, 'of which index nonzero', multi_nz)
    print(fn, 'side', side, 'tests', tested, 'groups', len(groups), 'groups with >1 index (VIOLATIONS of W)', viol,
          '| groups where the index vs the straight chord is nonzero (non-vacuity)', nontriv)
for fn in sys.argv[1:]:
    for side in (0, 1): run(fn, side, maxS=20000)
