"""PROVE 11-20: INDEPENDENT enumeration of step hexagons containing R=[0,u]x[0,u], NOT using Lemma hex's normal forms.
BFS over edge-adjacent lattice tiles (0-triangles, 1-triangles, squares; no other rhombus), 4 tiles incl. R, pairwise
interior-disjoint (separating axes), union a convex hexagon (6 corners, area = sum of tile areas), centrally or 3-fold symmetric.
Compares the resulting set of configurations with localenum.classify."""
import itertools, cmath
from fractions import Fraction as Fr
from collections import Counter
import localenum as LE
from claimU import interior_overlap
U = LE.UNITS
def tiles_on_edge(e):
    P, Q = sorted(e); d = LE.sub(Q, P); out = []
    if d[2:] == (0, 0):     # 0-edge: 0-triangles (third vertex in Lambda_0) and parallelograms with a 1-edge
        a, b = P[:2], Q[:2]; q = P[2:]
        for w in U:
            c = LE.add(a, w)
            if LE.sub(b, c) in U: out.append(('tri', frozenset(LE.L0(x) and (x + q) for x in (a, b, c))))
        for f in U:
            V = frozenset({P, Q, LE.add(P, LE.L1(f)), LE.add(Q, LE.L1(f))}); out.append(('par', V))
    else:
        a, b = P[2:], Q[2:]; p = P[:2]
        for w in U:
            c = LE.add(a, w)
            if LE.sub(b, c) in U: out.append(('tri', frozenset(p + x for x in (a, b, c))))
        for f in U:
            V = frozenset({P, Q, LE.add(P, LE.L0(f)), LE.add(Q, LE.L0(f))}); out.append(('par', V))
    res = []
    for k, V in out:
        if k == 'par':
            ar = LE.area(('x', V))
            if abs(ar - 1) < 1e-7: res.append(('sq', V))
            elif abs(ar - 0.5) < 1e-7: res.append(('rh', V))
        else: res.append(('tri', V))
    return res
def enum(u):
    R = ('rh', LE.P_((0, 0, 0, 0), [LE.L0(u), LE.L1(u)]))
    found = set(); frontier = {frozenset([R])}
    for size in (2, 3, 4):
        new = set()
        for S in frontier:
            for t in S:
                for e in LE.edges_of(t):
                    for nb in tiles_on_edge(e):
                        if nb in S or nb[0] == 'rh': continue
                        if any(interior_overlap(nb, s) for s in S): continue
                        new.add(S | {nb})
        frontier = new
    for S in frontier:
        if Counter(t[0] for t in S) != Counter(rh=1, sq=1, tri=2): continue
        V = set().union(*[t[1] for t in S]); zs = LE.convex_hull([LE.phi(P) for P in V])
        corners = [P for P in V if LE.is_corner(LE.phi(P), zs)]
        if len(corners) != 6: continue
        if abs(sum(LE.area(t) for t in S) - LE.polyarea(zs)) > 1e-7: continue
        c = tuple(sum(Fr(P[i]) for P in corners)/6 for i in range(4))
        Vc = set(corners)
        if {LE.rot_about(P, c, 3) for P in Vc} == Vc: typ = 'a'
        elif {LE.rot_about(P, c, 2) for P in Vc} == Vc: typ = 'b'
        else: continue
        found.add((typ, frozenset(S)))
    return R, found
if __name__ == '__main__':
    tot = Counter()
    for u in U:
        R, F2 = enum(u)
        _, Hs = LE.classify('B', u)
        F1 = {(typ[0], H) for typ, H in Hs}
        tot[('u', u, '#indep', len(F2), '#localenum', len(F1), 'equal', F1 == F2)] += 1
    for k in sorted(tot, key=str): print(tot[k], k)
