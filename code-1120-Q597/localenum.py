"""PROVE 11-20: ABSTRACT enumeration of every local step configuration (no boards).
A moving rhombus R = s0 x s1 with s0=[A,A+u] in Lambda_0, s1=[Bp,Bp+u] in Lambda_1 (same coordinate vector u; this
class is preserved by a-steps (translation) and b-steps (both lattices rotated by the same angle)).  Translation-invariance:
A = Bp = 0.  For each journey type T in {B,C}, each u in the six units, each step hexagon containing R as a tile
(type a0, a1: zonotope tilings; type b: tilings of tau+sigma by tau+q, p+sigma, R, square), compute:
  xrule  : every vertex of H has x >= xmin(R) (x = component along d_T)          [engine candidates()]
  across : the engine tie-break edge of R (horizontal edge at the unique leftmost vertex) is INTERIOR to H, or undefined
  k      : rotation (3 for a; for b the k in {2,4} whose image of R has an edge parallel/perp to d, None if not unique)
  DL     : every moved i-triangle moves from min h_i(pi(H)) to min+1           (direction law)
  MU     : mu_i = min h_i over pi(R) changes by 0 or -1, for i=0,1               (monotonicity)
  W1     : the new rhombus class is allowed (eta(u') != 0)                        (F1 preserved)
Output: a table of all configurations with xrule True, split by allowed/forbidden u."""
import cmath, math, itertools
from fractions import Fraction as Fr
from collections import Counter
ZETA = cmath.exp(1j*math.pi/3); z0 = cmath.exp(1j*math.pi/12); z1 = cmath.exp(-1j*math.pi/12)
def tri(a, b): return a + b*ZETA
def phi(P): return z0*tri(*P[:2]) + z1*tri(*P[2:])
UNITS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
def rot60(v): return (-v[1], v[0]+v[1])
def rotk(v, k):
    for _ in range(k % 6): v = rot60(v)
    return v
def add(P, Q): return tuple(a+b for a, b in zip(P, Q))
def sub(P, Q): return tuple(a-b for a, b in zip(P, Q))
def L0(v): return (v[0], v[1], 0, 0)
def L1(v): return (0, 0, v[0], v[1])
DIR = {'B': phi((0, 0, 0, 1)), 'C': -phi((0, 1, 0, 0))}
DIR = {k: v/abs(v) for k, v in DIR.items()}
ETA = {'B': lambda v: v[1], 'C': lambda v: v[0]+v[1]}
EPS = 1e-9
def par(A, B): return P_(A, [B])
def P_(base, gens):
    """parallelogram / zonotope tile vertex set"""
    pts = {base}
    for g in gens: pts |= {add(p, g) for p in pts}
    return frozenset(pts)
def x(P, T): return (phi(P)*DIR[T].conjugate()).real
def edges_of(tile):
    kind, V = tile; V = list(V); out = []
    for P, Q in itertools.combinations(V, 2):
        d = sub(Q, P)
        if (d[2:] == (0, 0) and d[:2] in UNITS) or (d[:2] == (0, 0) and d[2:] in UNITS): out.append(frozenset((P, Q)))
    return out
def boundary_edges(tiles):
    c = Counter(e for t in tiles for e in edges_of(t)); return {e for e, n in c.items() if n == 1}, {e for e, n in c.items() if n == 2}
def area(tile):
    V = list(tile[1]); zs = [phi(P) for P in V]; c = sum(zs)/len(zs)
    zs = sorted(zs, key=lambda z: cmath.phase(z-c))
    return abs(sum((zs[i].conjugate()*zs[(i+1) % len(zs)]).imag for i in range(len(zs))))/2
def is_tiling(tiles):
    """4 tiles, every edge used once or twice, boundary = convex hexagon (6 corners), area sum = hull area"""
    bd, inner = boundary_edges(tiles)
    pts = set().union(*[t[1] for t in tiles])
    zs = [phi(P) for P in pts]
    # convex hull area
    hull = convex_hull(zs)
    if abs(sum(area(t) for t in tiles) - polyarea(hull)) > 1e-7: return False
    corners = [h for i, h in enumerate(hull) if abs(((h-hull[i-1]).conjugate()*(hull[(i+1) % len(hull)]-h)).imag) > 1e-9]
    if len(corners) != 6: return False
    # every boundary edge lies on the hull boundary
    for e in bd:
        P, Q = tuple(e); m = (phi(P)+phi(Q))/2
        if not on_hull(m, hull): return False
    return True
def convex_hull(zs):
    pts = sorted(set((round(z.real, 9), round(z.imag, 9)) for z in zs))
    def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lo, hi = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 1e-12: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 1e-12: hi.pop()
        hi.append(p)
    return [complex(*p) for p in lo[:-1] + hi[:-1]]
def polyarea(h): return abs(sum((h[i].conjugate()*h[(i+1) % len(h)]).imag for i in range(len(h))))/2
def on_hull(m, h):
    for i in range(len(h)):
        a, b = h[i], h[(i+1) % len(h)]
        if abs(((b-a).conjugate()*(m-a)).imag) < 1e-7 and -1e-9 <= ((m-a)/(b-a)).real <= 1+1e-9: return True
    return False
def tri0(cell, q): return ('tri', frozenset(add(L0(c), L1(q)) for c in cell))
def tri1(p, cell): return ('tri', frozenset(add(L0(p), L1(c)) for c in cell))
def cells_at(seg):
    a, b = seg; u = sub(b, a); out = []
    for w in (rot60(u), rotk(u, 5)):
        c = add(a, w)
        if sub(b, c) in UNITS: out.append((a, b, c))
    return out
def kindof(t):
    if t[0] != 'tri': return t[0]
    V = t[1]; return 'tri0' if len({P[2:] for P in V}) == 1 else 'tri1'
def classify(T, u):
    """all step hexagons containing R = [0,u]_0 x [0,u]_1"""
    R = ('rh', P_((0, 0, 0, 0), [L0(u), L1(u)]))
    out = []
    # type a0 / a1 : zonotopes Z(L0 u, L1 u, g), g in the 12 units, R a tile at offset 0 or -g relative to base
    for lat in (0, 1):
        for g in UNITS:
            G = L0(g) if lat == 0 else L1(g)
            for off in ((0, 0, 0, 0), tuple(-c for c in G)):
                # tiles: R at 0 (i.e. zonotope base = off... ) -> build zonotope base so that R = P(base', U0,U1)
                base = off      # zonotope base point
                U0, U1 = L0(u), L1(u)
                # two tilings of Z(U0,U1,G); keep those containing R
                for tiling in zon_tilings(base, U0, U1, G):
                    if R not in tiling: continue
                    if not is_tiling(tiling): continue
                    ks = [t for t in tiling if t[0] == 'tri']
                    out.append(('a%d' % lat, frozenset(tiling)))
    # type b: tau a cell of s0, sigma a cell of s1; q vertex of sigma, p vertex of tau, square = edge(tau) x edge(sigma)
    s0 = ((0, 0), u); s1 = ((0, 0), u)
    for tau in cells_at(s0):
        for sig in cells_at(s1):
            for q in sig:
                for p in tau:
                    for e0 in itertools.combinations(tau, 2):
                        for e1 in itertools.combinations(sig, 2):
                            S = ('sq', P_(add(L0(e0[0]), L1(e1[0])), [L0(sub(e0[1], e0[0])), L1(sub(e1[1], e1[0]))]))
                            if abs(area(S) - 1) > 1e-7: continue
                            tiling = [R, tri0(tau, q), tri1(p, sig), S]
                            if len(set(tiling)) == 4 and is_tiling(tiling): out.append(('b', frozenset(tiling)))
    return R, list(dict.fromkeys(out))
def zon_tilings(base, g1, g2, g3):
    """the two tilings of Z(g1,g2,g3) based at `base` (whatever angular order): tiles P(gi,gj)+o_k, offsets
    either (0,0,0)-pattern-with-exceptional or complement; we brute force offset choices and test tiling-ness"""
    gs = [g1, g2, g3]; res = []
    for offs in itertools.product((0, 1), repeat=3):
        tiles = []
        for k in range(3):
            i, j = [m for m in range(3) if m != k]
            o = base if offs[k] == 0 else add(base, gs[k])
            V = P_(o, [gs[i], gs[j]])
            kind = None
            tiles.append(V)
        tl = []
        for V in tiles:
            pts = list(V)
            a = {P[:2] for P in V}; b = {P[2:] for P in V}
            if len(a) == 1 or len(b) == 1:   # degenerate in one lattice -> a lozenge made of two triangles
                tl.append(('loz', V))
            else:
                tl.append(('par', V))
        res.append(tl)
    out = []
    for tl in res:
        tiles = []
        ok = True
        for kind, V in tl:
            if kind == 'par':
                ar = area(('x', V)); tiles.append(('sq' if abs(ar-1) < 1e-7 else 'rh', V))
            else:
                # split lozenge into its two unit triangles
                a = {P[:2] for P in V}; b = {P[2:] for P in V}
                if len(b) == 1:
                    q = next(iter(b)); A = sorted(a); cellsA = lozenge_cells(A)
                    if cellsA is None: ok = False; break
                    tiles += [tri0(c, q) for c in cellsA]
                else:
                    p = next(iter(a)); Bs = sorted(b); cellsB = lozenge_cells(Bs)
                    if cellsB is None: ok = False; break
                    tiles += [tri1(p, c) for c in cellsB]
        if ok and len(tiles) == 4 and Counter(t[0] for t in tiles) == Counter(rh=1, sq=1, tri=2): out.append(tiles)
    return out
def lozenge_cells(pts):
    """pts: 4 lattice points of a unit lozenge -> its two unit triangles (split along the SHORT diagonal)"""
    for P, Q in itertools.combinations(pts, 2):
        if sub(Q, P) in UNITS:
            others = [r for r in pts if r not in (P, Q)]
            if all(sub(r, P) in UNITS and sub(r, Q) in UNITS for r in others) is False: continue
            if len(others) == 2 and all(sub(r, P) in UNITS and sub(r, Q) in UNITS for r in others):
                return [(P, Q, others[0]), (P, Q, others[1])]
    return None
def rot_about(P, c, k):
    v = sub(P, c); w = rotk(v[:2], k) + rotk(v[2:], k); return add(c, w)
def analyse(T, R, typ, H):
    V = set().union(*[t[1] for t in H])
    xm = min(x(P, T) for P in R[1])
    xrule = all(x(P, T) >= xm - EPS for P in V)
    # hexagon corners and centre
    zs = convex_hull([phi(P) for P in V])
    corners = [P for P in V if any(abs(phi(P)-h) < 1e-7 for h in zs) and is_corner(phi(P), zs)]
    c = tuple(sum(Fr(P[i]) for P in corners)/len(corners) for i in range(4)) if len(corners) == 6 else None
    if typ.startswith('a'): k = 3
    else:
        ks = []
        for kk in (2, 4):
            Rv = [rot_about(P, c, kk) for P in R[1]]
            ok = False
            for P, Q in itertools.combinations(Rv, 2):
                d = sub(Q, P)
                if (d[2:] == (0, 0) and d[:2] in UNITS) or (d[:2] == (0, 0) and d[2:] in UNITS):
                    w = (phi(Q)-phi(P))*DIR[T].conjugate()
                    if abs(w.imag) < 1e-9 or abs(w.real) < 1e-9: ok = True
            if ok: ks.append(kk)
        k = ks[0] if len(ks) == 1 else None
    # across edge
    acr = across_edge(T, R)
    inner = boundary_edges(H)[1]
    across = None if acr is None else (acr in inner)
    res = dict(xrule=xrule, k=k, across=across)
    if k is None: return res
    new = {t: (t[0], frozenset(rot_about(P, c, k) for P in t[1])) for t in H}
    eta = ETA[T]
    h = {0: (lambda q: -eta(q)), 1: (lambda p: eta(p))}
    DL = True
    for t in H:
        if t[0] != 'tri': continue
        i = 0 if kindof(t) == 'tri0' else 1
        piH = {(P[2:] if i == 0 else P[:2]) for P in V}
        lo = min(h[i](z) for z in piH)
        b = next(iter(t[1])); a_ = next(iter(new[t][1]))
        before = h[i](b[2:] if i == 0 else b[:2]); after = h[i](a_[2:] if i == 0 else a_[:2])
        if not (before == lo and after == lo+1): DL = False
    Rn = new[R]
    MU = []
    for i in (0, 1):
        pr = lambda S: {(P[2:] if i == 0 else P[:2]) for P in S[1]}
        MU.append(min(h[i](z) for z in pr(Rn)) - min(h[i](z) for z in pr(R)))
    un = sub(*sorted({P[:2] for P in Rn[1]})[::-1])
    res.update(DL=DL, MU=tuple(MU), W1=(eta(un) != 0), sqmove=None)
    return res
def is_corner(z, hull):
    for i, h in enumerate(hull):
        if abs(z-h) < 1e-7:
            return abs(((h-hull[i-1]).conjugate()*(hull[(i+1) % len(hull)]-h)).imag) > 1e-9
    return False
def across_edge(T, R):
    """engine stepset.across: unique leftmost vertex (rotated frame) and the edge of R from it that is horizontal"""
    V = list(R[1]); xs = [x(P, T) for P in V]; m = min(xs)
    lm = [P for P, xx in zip(V, xs) if abs(xx-m) < 1e-9]
    if len(lm) != 1: return None
    P = lm[0]
    for Q in V:
        if Q == P: continue
        d = sub(Q, P)
        if not ((d[2:] == (0, 0) and d[:2] in UNITS) or (d[:2] == (0, 0) and d[2:] in UNITS)): continue
        w = (phi(Q)-phi(P))*DIR[T].conjugate()
        if abs(w.imag) < 1e-9: return frozenset((P, Q))
    return None
if __name__ == '__main__':
    tab = Counter()
    for T in 'BC':
        for u in UNITS:
            R, Hs = classify(T, u)
            allowed = ETA[T](u) != 0
            for typ, H in Hs:
                r = analyse(T, R, typ, set(H))
                key = (T, 'allowed' if allowed else 'FORBIDDEN', 'u', u, typ, 'xrule', r['xrule'], 'across', r['across'], 'k', r['k'],
                       'DL', r.get('DL'), 'MU', r.get('MU'), 'W1', r.get('W1'))
                tab[key] += 1
    for k in sorted(tab, key=str):
        if k[7] if False else True: print(tab[k], k)
