"""Independent enumerator of ALL fine mixed subdivisions of phi(P0+P1), P0 = a*Delta in Lambda0, P1 = b*Delta in Lambda1,
phi(x0,x1) = x0 + zeta*x1, zeta = exp(SIGN*i*pi/6).  Exact arithmetic in Z[sqrt3] (no floats anywhere).

Lattice points x = m + n*omega are encoded (m,n).  A labelled vertex is (m0,n0,m1,n1).
Tiles: tau x {w} (tau a unit triangle of P0), {v} x sigma, e0 x e1 (unit edges).
A tiling = set of tiles such that every ccw half-edge of a tile is matched by the reverse half-edge of exactly one other tile,
or is a boundary half-edge of the region; interiors pairwise disjoint.  (All tiles are positively oriented, so the covering has
degree 1: this is exactly an edge-to-edge tiling, i.e. a mixed subdivision since labels are read off planar position by lem:inj.)
"""
import sys, itertools
from functools import lru_cache

# ---------- Z[sqrt3] arithmetic: pairs (r,s) = r + s*sqrt3 ----------
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def sub(x, y): return (x[0]-y[0], x[1]-y[1])
def mul(x, y): return (x[0]*y[0] + 3*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def sgn(x):
    r, s = x
    if r >= 0 and s >= 0: return 0 if (r == 0 and s == 0) else 1
    if r <= 0 and s <= 0: return -1
    # opposite signs: compare r^2 vs 3 s^2
    d = r*r - 3*s*s
    return (1 if r > 0 else -1) if d > 0 else (-1 if r > 0 else 1)

SIGN = +1
def planar(v):
    """2*phi(v) as (X,Y) in Z[sqrt3]^2.  omega = -1/2 + i sqrt3/2; zeta = sqrt3/2 + SIGN i/2."""
    m0, n0, m1, n1 = v
    # 2*x0 = (2m0 - n0) + i*sqrt3*n0
    X0, Y0 = (2*m0 - n0, 0), (0, n0)
    # 2*zeta*x1: x1 = u + i w with 2u = 2m1-n1, 2w = sqrt3 n1 ; zeta = (sqrt3 + SIGN i)/2
    # 2*zeta*x1 = (sqrt3 + SIGN i) * (u + i w) = (sqrt3 u - SIGN w) + i (SIGN u + sqrt3 w)
    # with u = (2m1-n1)/2, w = sqrt3 n1/2:  real = sqrt3(2m1-n1)/2 - SIGN sqrt3 n1/2 ; imag = SIGN(2m1-n1)/2 + 3 n1/2
    # multiply everything by 2 once more to stay integral -> we return 4*phi
    X = (2*(2*m0 - n0), 2*m1 - n1 - SIGN*n1)
    Y = (SIGN*(2*m1 - n1) + 3*n1, 2*n0)
    return (X, Y)

def cross(o, a, b):
    """sign of (a-o) x (b-o) for labelled vertices."""
    O, A, B = planar(o), planar(a), planar(b)
    ax, ay = sub(A[0], O[0]), sub(A[1], O[1]); bx, by = sub(B[0], O[0]), sub(B[1], O[1])
    return sgn(sub(mul(ax, by), mul(ay, bx)))

def tri_pts(a):  # lattice points of a*Delta = conv(0, a, a*(1+omega))?  -> use conv(0, a, -a*omega^2) = conv(0,a,a+a*omega)
    return [(m, n) for m in range(a+1) for n in range(m+1)]   # points m + n*omega with 0<=n<=m<=a: triangle conv(0, a, a(1+omega))

def unit_tris(a):
    P = set(tri_pts(a)); out = []
    for (m, n) in P:
        for T in (((m, n), (m+1, n), (m+1, n+1)), ((m, n), (m, n+1), (m+1, n+1))):
            if all(t in P for t in T): out.append(frozenset(T))
    return out

def unit_edges(a):
    E = set()
    for T in unit_tris(a):
        for e in itertools.combinations(sorted(T), 2): E.add(frozenset(e))
    return sorted(E, key=sorted)

def ccw(vs):
    """order the vertices of a convex polygon (given as labelled vertices) counterclockwise."""
    vs = list(vs)
    # gift wrap from the lowest by exact comparisons; n <= 4 so brute force over permutations
    for perm in itertools.permutations(vs[1:]):
        cyc = [vs[0]] + list(perm)
        if all(cross(cyc[i], cyc[(i+1) % len(cyc)], cyc[(i+2) % len(cyc)]) > 0 for i in range(len(cyc))):
            return tuple(cyc)
    raise ValueError('degenerate tile %r' % (vs,))

def make_tiles(a, b):
    T0, T1 = unit_tris(a), unit_tris(b); E0, E1 = unit_edges(a), unit_edges(b)
    A0, A1 = tri_pts(a), tri_pts(b)
    tiles = []
    for t in T0:
        for w in A1: tiles.append(('t0', frozenset(x + w for x in t)))
    for s in T1:
        for v in A0: tiles.append(('t1', frozenset(v + y for y in s)))
    for e in E0:
        for f in E1: tiles.append(('par', frozenset(x + y for x in e for y in f)))
    return tiles

def overlap(P, Q):
    """interiors of convex ccw polygons P,Q intersect?  SAT on edge normals."""
    for A, B in ((P, Q), (Q, P)):
        for i in range(len(A)):
            o, c = A[i], A[(i+1) % len(A)]
            if all(cross(o, c, x) <= 0 for x in B): return False
    return True

def boundary_halfedges(a, b):
    """ccw boundary half-edges of the region: the hexagon's sides, subdivided at labelled lattice points.
    Found as: half-edges of tiles that lie on a supporting line of the whole point set."""
    return None  # computed in enumerate() from tile half-edges

def enumerate_all(a, b, verbose=False):
    tiles = make_tiles(a, b)
    poly = [ccw(vs) for _, vs in tiles]
    allpts = set(v for _, vs in tiles for v in vs)
    n = len(tiles)
    # half-edges
    HE = {}  # halfedge -> list of tile indices having it (ccw)
    for i, P in enumerate(poly):
        for j in range(len(P)):
            HE.setdefault((P[j], P[(j+1) % len(P)]), []).append(i)
    # boundary half-edges: all points weakly left of the directed line, i.e. no point strictly right
    bd = set()
    for h in HE:
        u, v = h
        if all(cross(u, v, x) >= 0 for x in allpts): bd.add(h)
    # conflict bitsets
    conf = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if overlap(poly[i], poly[j]): conf[i] |= 1 << j; conf[j] |= 1 << i
    if verbose: print('tiles', n, 'boundary half-edges', len(bd), file=sys.stderr)
    results = []
    # state: needed = set of half-edges that must be supplied by a new tile; used = half-edges already supplied
    def rec(needed, used, chosen, forbid):
        if not needed:
            results.append(frozenset(chosen)); return
        # most constrained needed half-edge
        best, bestc = None, None
        for h in needed:
            c = [i for i in HE.get(h, []) if not (forbid >> i) & 1]
            if best is None or len(c) < len(bestc):
                best, bestc = h, c
                if len(c) <= 1: break
        for i in bestc:
            P = poly[i]; ok = True; nn = set(needed); uu = set(used)
            for j in range(len(P)):
                h = (P[j], P[(j+1) % len(P)])
                if h in uu: ok = False; break
                uu.add(h)
                if h in nn: nn.discard(h)
                elif h in bd: pass  # cannot happen unless boundary already covered... handled by uu
                else:
                    r = (h[1], h[0])
                    if r in uu: ok = False; break  # reverse already supplied but was not needed -> inconsistent
                    nn.add(r)
            if not ok: continue
            chosen.append(i)
            rec(nn, uu, chosen, forbid | conf[i] | (1 << i))
            chosen.pop()
    # boundary half-edges are "needed" from the start; a boundary half-edge h is supplied by a tile having h.
    rec(set(bd), set(), [], 0)
    return tiles, poly, results

if __name__ == '__main__':
    a, b = int(sys.argv[1]), int(sys.argv[2])
    if len(sys.argv) > 3: SIGN = int(sys.argv[3])
    tiles, poly, R = enumerate_all(a, b, verbose=True)
    print(a, b, 'tilings', len(R))

# ---------- the mixed-subdivision condition (labelled intersections) ----------
def parts(vs):
    return frozenset(v[:2] for v in vs), frozenset(v[2:] for v in vs)
def is_mixed(tiles, T):
    """For every pair of cells B=B0+B1, B'=B0'+B1': common labelled vertices == (B0 cap B0') x (B1 cap B1')."""
    cells = [tiles[i][1] for i in T]; P = [parts(c) for c in cells]
    for i in range(len(cells)):
        for j in range(i+1, len(cells)):
            common = cells[i] & cells[j]
            prod = frozenset(x + y for x in P[i][0] & P[j][0] for y in P[i][1] & P[j][1])
            if common != prod: return False
    return True

# ---------- independent check: Cayley tetrahedra form a triangulation (pseudomanifold + volume) ----------
def cay(pt, level):
    """3D Cayley point, scaled: (4*phi_level(pt), level*K).  Uses planar() with the other part zero."""
    if level == 0: X, Y = planar(pt + (0, 0))
    else: X, Y = planar((0, 0) + pt)
    return (X, Y, ((4 if level else 0), 0))
def det3(u, v, w):
    a = sub(mul(u[1], w[2]) if False else mul(v[1], w[2]), mul(v[2], w[1]))
    b = sub(mul(v[0], w[2]), mul(v[2], w[0]))
    c = sub(mul(v[0], w[1]), mul(v[1], w[0]))
    return add(sub(mul(u[0], a), mul(u[1], b)), mul(u[2], c))
def vsub(p, q): return tuple(sub(p[k], q[k]) for k in range(3))
def tets(tiles, T):
    out = []
    for i in T:
        B0, B1 = parts(tiles[i][1])
        out.append(tuple([cay(x, 0) for x in sorted(B0)] + [cay(y, 1) for y in sorted(B1)]))
    return out
def vol6(t):
    d = det3(vsub(t[1], t[0]), vsub(t[2], t[0]), vsub(t[3], t[0]))
    return d if sgn(d) >= 0 else (-d[0], -d[1])
def cayley_ok(tiles, T, a, b):
    """(1) every triangle face of a tetrahedron is either on the boundary of conv(Cayley) or shared by exactly two tets on
    opposite sides; (2) total volume = volume of conv(Cayley) (computed from a reference: any one triangulation's volume
    is passed in via the global REFVOL; we check equality across all tilings instead)."""
    TT = tets(tiles, T); faces = {}
    for t in TT:
        for k in range(4):
            f = tuple(sorted(t[:k] + t[k+1:])); apex = t[k]
            s = sgn(det3(vsub(f[1], f[0]), vsub(f[2], f[0]), vsub(apex, f[0])))
            if s == 0: return None, 'degenerate tet'
            faces.setdefault(f, []).append(s)
    bad = [f for f, ss in faces.items() if not (len(ss) == 1 or (len(ss) == 2 and ss[0] != ss[1]))]
    if bad: return None, 'face condition'
    V = (0, 0)
    for t in TT: V = add(V, vol6(t))
    return V, 'ok'
