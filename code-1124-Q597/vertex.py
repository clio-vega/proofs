"""Angles at the vertex o = min R* on J_o's board at e+1 (abstract canonical configuration): which tiles of CJ+RJ contain
o and the angular interval each covers; the uncovered gap; also the angular sectors of the two universal triangles Z1,Z2
and of the straight configuration S.   usage: python3 vertex.py"""
import sys, cmath, math
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import phi
STRT = {'B': 'a1', 'C': 'a0'}
Z = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
     'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
def sector(t, v):
    """angular sector (deg, ccw from start) of tile t at vertex v, or None"""
    pts = [P for P in t[1] if tuple(P) != tuple(v)]
    if len(pts) == len(t[1]): return None
    zv = phi(v); zs = sorted(phi(P) for P in pts) if False else [phi(P) for P in pts]
    # neighbours of v along edges = the two points at unit distance with the tile's edges: take the two whose angle is extreme
    angs = sorted((math.degrees(cmath.phase(z - zv)) % 360, z) for z in zs if abs(abs(z - zv) - 1) < 1e-6)
    a = [x for x, _ in angs]
    # the interior sector is the arc between the two edge directions that contains the centroid direction
    cdir = math.degrees(cmath.phase(sum(phi(P) for P in t[1]) / len(t[1]) - zv)) % 360
    for s, e_ in ((a[0], a[1]), (a[1], a[0])):
        w = (e_ - s) % 360
        if (cdir - s) % 360 < w: return (round(s, 1), round(w, 1))
for Tp in 'CB':
    Rs = CF[Tp][(k0, 'b')][0]; HJ, RJ, CJ = step(Tp, Rs, 'b'); S, _, _ = step(Tp, RJ, STRT[Tp])
    o = min(Rs[1]); print('==== J_o toward', Tp, 'o =', o)
    vtx = o if Tp == 'C' else tuple(a + b for a, b in zip(o, (0, 1, 0, 1)))
    print(' vertex used:', vtx)
    for lab, tiles in (('CJ+RJ', CJ | {RJ}), ('HJ (before)', HJ), ('S straight', S), ('Z', {('tri', frozenset(tuple(x + y for x, y in zip(P, o)) for P in z)) for z in Z[Tp]})):
        for t in tiles:
            s = sector(t, vtx)
            if s: print('  ', lab, t[0], 'sector start %.1f width %.1f' % s)
