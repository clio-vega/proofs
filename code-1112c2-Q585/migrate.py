"""Migration on mosaics (Purbhoo 0705.1184 §3.1), exact Z^4 geometry; floats only for left/right/horizontal tests."""
from mosaic import *
from fractions import Fraction as Fr
from collections import Counter
SIGN, RULE = -1, 0          # the unique geometrically valid convention (geotest.py: 41/41 convex, area-exact)
G = Geo(SIGN)
def fr(P): return tuple(Fr(x) for x in P)
def prot(P, k):             # rotate by k*60 degrees (exact, works on Fractions)
    a0, b0, a1, b1 = P
    for _ in range(k % 6): a0, b0 = -b0, a0+b0; a1, b1 = -b1, a1+b1
    return (a0, b0, a1, b1)
def norm(P): return tuple(int(x) if (not isinstance(x, Fr) or x.denominator == 1) else x for x in P)

Global_ties = [0, 0]   # [ties met, ties left unresolved]
class Board:
    def __init__(self, M):
        self.tiles = {}; self.edges = {}; self.nid = 0
        for kind, vs in M['tiles']: self.add(kind, vs)
        self.nest = []    # per nest: frame + positions
        for X, rh in zip(M['corners'], M['nests']):
            u = tuple(-x for x in X['s_first']); w = X['s_second']
            info = dict(X=X['X'], u=u, w=w, pos={})
            for _, vs in rh:
                t = self.add('rh', vs); ij = self.packed_at(info, vs)
                assert ij is not None and ij not in info['pos']; info['pos'][ij] = t
            self.nest.append(info)
        for info in self.nest: assert self.is_ideal(set(info['pos']))
    def add(self, kind, vs):
        t = self.nid; self.nid += 1; vs = tuple(norm(v) for v in vs); self.tiles[t] = (kind, vs)
        for i in range(len(vs)): self.edges.setdefault(frozenset((vs[i], vs[(i+1) % len(vs)])), set()).add(t)
        return t
    def remove(self, t):
        kind, vs = self.tiles.pop(t)
        for i in range(len(vs)):
            e = frozenset((vs[i], vs[(i+1) % len(vs)])); self.edges[e].discard(t)
            if not self.edges[e]: del self.edges[e]
        return kind, vs
    @staticmethod
    def coef(D, v):
        """D = k*v for a Z^4 unit step v (one component zero); return k (Fraction) or None"""
        k = None
        for a, b in zip(D, v):
            if b == 0:
                if a != 0: return None
            else:
                kk = Fr(a) / b
                if k is None: k = kk
                elif kk != k: return None
        return k
    def packed_at(self, info, vs):
        X, u, w = info['X'], info['u'], info['w']
        S = set(vs); best = None
        for p in vs:
            if {p, padd(p, u), padd(p, w), padd(padd(p, u), w)} == S:
                D = psub(p, X)
                Du = (D[0], D[1], 0, 0) if u[2:] == (0, 0) else (0, 0, D[2], D[3])
                Dw = psub(D, Du)
                i = self.coef(Du, u) if any(Du) else Fr(0); j = self.coef(Dw, w) if any(Dw) else Fr(0)
                if i is None or j is None: continue
                if i >= 0 and j >= 0 and i.denominator == 1 and j.denominator == 1: return (int(i), int(j))
        return None
    @staticmethod
    def is_ideal(S):
        return all(((i-1, j) in S or i == 0) and ((i, j-1) in S or j == 0) for (i, j) in S)
    def neighbours(self, t):
        kind, vs = self.tiles[t]; out = set()
        for i in range(len(vs)):
            out |= self.edges.get(frozenset((vs[i], vs[(i+1) % len(vs)])), set())
        out.discard(t); return out
    # ------------------------------------------------------------- one migration step
    def step(self, R, d):
        """d: unit complex = migration direction.  Rotate the minimal hexagon; return new tid of the rhombus."""
        rot = d.conjugate()
        def x(P): return (G.phi(P)*rot).real
        xmin = min(x(P) for P in self.tiles[R][1])
        # connected subsets containing R of size 3,4 with exactly one rhombus
        cands = []
        frontier = {frozenset([R])}
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
                cyc = [cyc[i] for i in range(m_)
                       if abs(((zz[i]-zz[i-1]).conjugate()*(zz[(i+1) % m_]-zz[i])).imag) > 1e-9]
                if len(cyc) != 6: continue
                if min(x(P) for t in S for P in self.tiles[t][1]) < xmin - 1e-9: continue
                c = tuple(sum(Fr(P[k]) for P in cyc) / 6 for k in range(4))
                V = set(fr(P) for P in cyc)
                def rotc(P, k): return tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k)))
                if set(rotc(P, 3) for P in V) == V: typ = 'a'
                elif set(rotc(P, 2) for P in V) == V: typ = 'b'
                else: continue
                cands.append((size, typ, S, c))
        if not cands: return None
        # tie-break (interpretation, flagged): the rhombus must not move left ("moves roughly to the right,
        # although sometimes strictly vertically", §3.1).  Count how often it is needed.
        def newx(z):
            size, typ, S, c = z
            ks = [3] if typ == 'a' else [2, 4]
            vals = []
            for k in ks:
                vs = [tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k))) for P in self.tiles[R][1]]
                vals.append(sum(x(P) for P in vs)/4)
            return max(vals)
        x0 = sum(x(P) for P in self.tiles[R][1])/4
        if len(cands) > 1:
            # tie-break (interpretation of Figs (a),(b), flagged): the hexagon contains the tile across the
            # HORIZONTAL edge of the rhombus that starts at its leftmost vertex.
            Global_ties[0] += 1
            vs = self.tiles[R][1]; z = [G.phi(P)*rot for P in vs]
            xs = [w.real for w in z]; lm = min(range(4), key=lambda i: xs[i])
            if sum(1 for w in xs if abs(w-xs[lm]) < 1e-9) == 1:
                for nb_i in ((lm+1) % 4, (lm-1) % 4):
                    if abs((z[nb_i]-z[lm]).imag) < 1e-9:
                        e = frozenset((vs[lm], vs[nb_i])); across = self.edges.get(e, set()) - {R}
                        cands = [zz for zz in cands if across & zz[2]]
            if len(cands) != 1: Global_ties[1] += 1
        if not cands: return None
        m = min(cands, key=lambda z: z[0])[0]
        best = [z for z in cands if z[0] == m]
        if len(best) != 1:
            self.tie = [(z[1], sorted(round((G.phi(P)*rot).real,3) for t in z[2] for P in self.tiles[t][1])[:2], sorted(self.tiles[t][0] for t in z[2]), z[2]) for z in best]
            raise RuntimeError('non-unique minimal hexagon')
        size, typ, S, c = best[0]
        def rotc(P, k): return norm(tuple(a+b for a, b in zip(c, prot(tuple(Fr(p)-q for p, q in zip(P, c)), k))))
        if typ == 'a': ks = [3]
        else:
            ks = []
            for k in (2, 4):
                vs = [rotc(P, k) for P in self.tiles[R][1]]
                z = [G.phi(P)*rot for P in vs]
                if any(abs((z[i]-z[(i+1) % 4]).imag) < 1e-9 or abs((z[i]-z[(i+1) % 4]).real) < 1e-9 for i in range(4)):
                    ks.append(k)
            assert len(ks) == 1, ('120-rotation rule ambiguous', ks)
        k = ks[0]; newR = None
        old = [(t, self.remove(t)) for t in S]
        for t, (kind, vs) in old:
            nt = self.add(kind, tuple(rotc(P, k) for P in vs))
            if t == R: newR = nt
        self.last = typ
        return newR
    def boundary_cycle(self, S):
        cnt = Counter()
        for t in S:
            vs = self.tiles[t][1]
            for i in range(len(vs)): cnt[frozenset((vs[i], vs[(i+1) % len(vs)]))] += 1
        bd = [e for e, c in cnt.items() if c == 1]
        adj = {}
        for e in bd:
            p, q = tuple(e); adj.setdefault(p, []).append(q); adj.setdefault(q, []).append(p)
        if any(len(v) != 2 for v in adj.values()): return None
        start = next(iter(adj)); cyc = [start]; prev = None; cur = start
        while True:
            nxt = [q for q in adj[cur] if q != prev][0] if prev is not None else adj[cur][0]
            if nxt == start: break
            cyc.append(nxt); prev, cur = cur, nxt
            if len(cyc) > 20: return None
        return cyc if len(cyc) == len(bd) else None
    # ------------------------------------------------------------- one rhombus journey
    def journey(self, R, src, d, maxsteps=2000):
        steps = Counter()
        for _ in range(maxsteps):
            vs = self.tiles[R][1]
            for X, info in enumerate(self.nest):
                if X == src: continue
                ij = self.packed_at(info, vs)
                if ij is not None and ij not in info['pos'] and self.is_ideal(set(info['pos']) | {ij}):
                    info['pos'][ij] = R; return R, X, ij, steps
            R2 = self.step(R, d)
            if R2 is None: raise RuntimeError('stuck: no hexagon to the right')
            steps[self.last] += 1; R = R2
        raise RuntimeError('journey too long')
