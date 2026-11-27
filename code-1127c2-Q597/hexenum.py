"""(Hex), exhaustively: every edge-connected, interior-disjoint cluster {rhombus, square, triangle, triangle} of Z^4 unit
cells (three forms) that the engine's acceptance test would accept (boundary with 6 turning corners; corner set invariant
under the 180 deg (a) or 120 deg (b) rotation about its centre) has: boundary = 6 unit edges, the rotation maps the full
boundary edge set to itself, and the hexagon is convex.  Up to translation, every rhombus orientation, every growth order."""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
from migrate import prot
G = MZ.Geo(-1)
U = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
def add(P, Q): return tuple(a + b for a, b in zip(P, Q))
def sub(P, Q): return tuple(a - b for a, b in zip(P, Q))
def L0(d): return (d[0], d[1], 0, 0)
def L1(d): return (0, 0, d[0], d[1])
def ph(P): return G.phi(P)
def ccw(cyc):
    z = [ph(P) for P in cyc]; a = sum((z[i].conjugate() * z[(i + 1) % len(z)]).imag for i in range(len(z)))
    return tuple(cyc) if a > 0 else tuple(reversed(cyc))
def cell_ef(p, e, f):
    return ccw([p, add(p, L0(e)), add(add(p, L0(e)), L1(f)), add(p, L1(f))])
def kind(cyc):
    if len(cyc) == 3: return 'tri'
    z = [ph(P) for P in cyc]; v1 = z[1] - z[0]; v2 = z[3] - z[0]
    return 'sq' if abs((v1.conjugate() * v2).real) < 1e-9 else 'rh'
def tiles_on_edge(P, Q):
    d = sub(Q, P); out = []
    if d[2:] == (0, 0):
        d0 = d[:2]; i = U.index(d0)
        for r in (U[(i + 1) % 6], U[(i - 1) % 6]):          # the two unit triangles on segment P0Q0
            third = add(P, L0(r)) if add(d0, (0, 0)) != r else None
            for c in (add(P, L0(U[(i + 1) % 6])), add(P, L0(U[(i - 1) % 6]))): pass
        out += [ccw([P, Q, add(P, L0(U[(i + 1) % 6]))]), ccw([P, Q, add(P, L0(U[(i - 1) % 6]))])]
        for f in U: out.append(cell_ef(P, d0, f))
    else:
        d1 = d[2:]; i = U.index(d1)
        out += [ccw([P, Q, add(P, L1(U[(i + 1) % 6]))]), ccw([P, Q, add(P, L1(U[(i - 1) % 6]))])]
        for e in U: out.append(cell_ef(P, e, d1))
    return out
def edges(cyc): return [frozenset((cyc[i], cyc[(i + 1) % len(cyc)])) for i in range(len(cyc))]
def overlap(a, b):
    pa = [ph(P) for P in a]; pb = [ph(P) for P in b]
    for poly in (pa, pb):
        for i in range(len(poly)):
            e = poly[(i + 1) % len(poly)] - poly[i]; n = complex(-e.imag, e.real)
            ra = [(z * n.conjugate()).real for z in pa]; rb = [(z * n.conjugate()).real for z in pb]
            if max(ra) <= min(rb) + 1e-9 or max(rb) <= min(ra) + 1e-9: return False
    return True
def key(c): return frozenset(c)
def boundary(cl):
    from collections import Counter
    cnt = Counter(e for c in cl for e in edges(c)); return {e for e, k in cnt.items() if k == 1}
def cycle_of(bd):
    adj = {}
    for e in bd:
        p, q = tuple(e); adj.setdefault(p, []).append(q); adj.setdefault(q, []).append(p)
    if any(len(v) != 2 for v in adj.values()): return None
    start = next(iter(adj)); cyc = [start]; prev = None; cur = start
    while True:
        nxt = [q for q in adj[cur] if q != prev][0] if prev is not None else adj[cur][0]
        if nxt == start: break
        cyc.append(nxt); prev, cur = cur, nxt
    return cyc if len(cyc) == len(bd) else None
def rotc(P, c, k): return tuple(a + b for a, b in zip(c, prot(tuple(Fr(p) - q for p, q in zip(P, c)), k)))
# all rhombus cells at origin
RH = {key(c): c for e in U for f in U for c in [cell_ef((0, 0, 0, 0), e, f)] if kind(c) == 'rh'}
seen = set(); acc = 0; bad = []; stats = {}
frontier = [(c,) for c in RH.values()]
for size in (2, 3, 4):
    new = {}
    for cl in frontier:
        for e in boundary(cl):
            P, Q = tuple(e)
            for t in tiles_on_edge(P, Q):
                if key(t) in {key(c) for c in cl}: continue
                if kind(t) == 'rh': continue
                if any(overlap(t, c) for c in cl): continue
                ncl = cl + (t,)
                kinds = sorted(kind(c) for c in ncl)
                if kinds.count('sq') > 1 or kinds.count('tri') > 2: continue
                new[frozenset(key(c) for c in ncl)] = ncl
    frontier = list(new.values())
for cl in frontier:
    if sorted(kind(c) for c in cl) != ['rh', 'sq', 'tri', 'tri']: continue
    bd = boundary(cl); cyc = cycle_of(bd)
    if cyc is None: continue
    zz = [ph(P) for P in cyc]; m = len(zz)
    corners = [cyc[i] for i in range(m) if abs(((zz[i] - zz[i - 1]).conjugate() * (zz[(i + 1) % m] - zz[i])).imag) > 1e-9]
    if len(corners) != 6: continue
    c = tuple(sum(Fr(P[k]) for P in corners) / 6 for k in range(4)); V = set(tuple(Fr(x) for x in P) for P in corners)
    typ = 'a' if {rotc(P, c, 3) for P in V} == V else ('b' if {rotc(P, c, 2) for P in V} == V else None)
    if typ is None: continue
    acc += 1
    ks = [3] if typ == 'a' else [2, 4]
    bdF = {frozenset(tuple(Fr(x) for x in P) for P in e) for e in bd}
    ok_bd = all({frozenset(rotc(P, c, k) for P in e) for e in bdF} == bdF for k in ks)
    zc = [ph(P) for P in corners]
    turns = [((zc[i] - zc[i - 1]).conjugate() * (zc[(i + 1) % 6] - zc[i])).imag for i in range(6)]
    convex = all(t > 0 for t in turns) or all(t < 0 for t in turns)
    s = (typ, len(bd) == 6, ok_bd, convex); stats[s] = stats.get(s, 0) + 1
print('accepted clusters (with rhombus at origin, all orientations):', acc)
for k, v in sorted(stats.items()): print(v, '(type, 6 unit edges, rotation maps boundary to itself, convex) =', k)
