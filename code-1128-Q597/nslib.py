"""PROVE 2026-11-28: shared instrument for (NS).  cands() is a verbatim copy of the candidate part of
migrate.Board.step (S1); across_edge() is the (S2) across edge.  gash() computes the puzzle picture of a board."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import migrate
from fractions import Fraction as Fr
from collections import Counter
G = migrate.G; prot = migrate.prot; fr = migrate.fr
def psi(P): return (P[0] + P[2], P[1] + P[3])
def cands(self, R, d):
    rot = d.conjugate()
    def x(P): return (G.phi(P)*rot).real
    xmin = min(x(P) for P in self.tiles[R][1])
    out = []; frontier = {frozenset([R])}
    for size in (2, 3, 4):
        new = set()
        for S in frontier:
            for t in S:
                for nb in self.neighbours(t):
                    if nb in S or self.tiles[nb][0] == 'rh': continue
                    new.add(S | {nb})
        frontier = new
        for S in frontier:
            kinds = Counter(self.tiles[t][0] for t in S)
            if not (size == 4 and kinds == Counter(rh=1, tri=2, sq=1)): continue
            cyc = self.boundary_cycle(S)
            if cyc is None: continue
            zz = [G.phi(P) for P in cyc]; m_ = len(zz)
            cyc = [cyc[i] for i in range(m_) if abs(((zz[i]-zz[i-1]).conjugate()*(zz[(i+1) % m_]-zz[i])).imag) > 1e-9]
            if len(cyc) != 6: continue
            if min(x(P) for t in S for P in self.tiles[t][1]) < xmin - 1e-9: continue
            c = tuple(sum(Fr(P[k]) for P in cyc) / 6 for k in range(4))
            V = set(fr(P) for P in cyc)
            def rotc(P, k): return tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k)))
            if set(rotc(P, 3) for P in V) == V: typ = 'a'
            elif set(rotc(P, 2) for P in V) == V: typ = 'b'
            else: continue
            out.append((typ, S, c))
    return out
def across_edge(self, R, d):
    """(S2): the edge [v,v'] of R at its unique x-minimal vertex, parallel to d; None if undefined."""
    rot = d.conjugate(); vs = self.tiles[R][1]; z = [G.phi(P)*rot for P in vs]
    xs = [w.real for w in z]; lm = min(range(4), key=lambda i: xs[i])
    if sum(1 for w in xs if abs(w-xs[lm]) < 1e-9) != 1: return None
    for nb_i in ((lm+1) % 4, (lm-1) % 4):
        if abs((z[nb_i]-z[lm]).imag) < 1e-9: return frozenset((vs[lm], vs[nb_i]))
    return None
def classify(self, R, d):
    """list of candidates with verdict good/bad (bad = across edge defined and on the boundary of H)"""
    e = across_edge(self, R, d); res = []
    for typ, S, c in cands(self, R, d):
        if e is None: res.append(('good', typ, S, c)); continue
        across = self.edges.get(e, set()) - {R}
        res.append(('good' if across & S else 'bad', typ, S, c))
    return e, res
# ---------------------------------------------------------------- puzzle picture
def vsub(a, b): return tuple(x - y for x, y in zip(a, b))
def lab(v):   # Z^4 edge vector -> 0 (Lambda_0), 1 (Lambda_1), None (mixed)
    if not any(v[2:]): return 0
    if not any(v[:2]): return 1
    return None
def unit_tris(kind, vs):
    """psi-image of a tri/sq tile as unit lattice triangles (frozensets of 3 points), plus its labelled boundary
    edges {(frozenset of 2 puzzle points): label}"""
    pts = [psi(P) for P in vs]; edges = {}
    for i in range(len(vs)):
        edges[frozenset((pts[i], pts[(i+1) % len(vs)]))] = lab(vsub(vs[(i+1) % len(vs)], vs[i]))
    if kind == 'tri': return [frozenset(pts)], edges
    # square -> puzzle rhombus: split along its short diagonal (the two 60-degree... ) -- the diagonal of unit length
    a, b, c, e = pts
    def d2(p, q): x, y = p[0]-q[0], p[1]-q[1]; return x*x + x*y + y*y   # squared length in the triangular lattice
    if d2(a, c) == 1: return [frozenset((a, b, c)), frozenset((a, c, e))], edges
    assert d2(b, e) == 1, (pts,)
    return [frozenset((b, c, e)), frozenset((b, e, a))], edges
def gash(self, R):
    """Return (ok_segment, P, M, Q, u, R-edge images) for the moving rhombus."""
    vs = self.tiles[R][1]
    # base b: the vertex with both edges going 'out' with the same coordinate vector
    for i in range(4):
        p, q, r = vs[i-1], vs[i], vs[(i+1) % 4]
        e1, e2 = vsub(p, q), vsub(r, q)
        if psi(e1) == psi(e2) and {lab(e1), lab(e2)} == {0, 1}:
            u = psi(e1); P = psi(q)
            return True, P, (P[0]+u[0], P[1]+u[1]), (P[0]+2*u[0], P[1]+2*u[1]), u, q
    return False, None, None, None, None, None
def cands_nox(self, R):
    """cands() with the x-rule switched off (acceptance = 4 tiles, 6-corner boundary, rotation symmetry)."""
    import cmath
    return cands(self, R, 1+0j) if False else _cands_generic(self, R)
def _cands_generic(self, R):
    out = []; frontier = {frozenset([R])}
    for size in (2, 3, 4):
        new = set()
        for S in frontier:
            for t in S:
                for nb in self.neighbours(t):
                    if nb in S or self.tiles[nb][0] == 'rh': continue
                    new.add(S | {nb})
        frontier = new
        for S in frontier:
            kinds = Counter(self.tiles[t][0] for t in S)
            if not (size == 4 and kinds == Counter(rh=1, tri=2, sq=1)): continue
            cyc = self.boundary_cycle(S)
            if cyc is None: continue
            zz = [G.phi(P) for P in cyc]; m_ = len(zz)
            cyc = [cyc[i] for i in range(m_) if abs(((zz[i]-zz[i-1]).conjugate()*(zz[(i+1) % m_]-zz[i])).imag) > 1e-9]
            if len(cyc) != 6: continue
            c = tuple(sum(Fr(P[k]) for P in cyc) / 6 for k in range(4))
            V = set(fr(P) for P in cyc)
            def rotc(P, k): return tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k)))
            if set(rotc(P, 3) for P in V) == V: typ = 'a'
            elif set(rotc(P, 2) for P in V) == V: typ = 'b'
            else: continue
            out.append((typ, S, c))
    return out
