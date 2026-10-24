r"""
Q569 -- Pylyavskyy-Yang hexagonal puzzle tiles: boundary piece or bulk piece?

Tile extraction is VERBATIM from the tikz of
  library/2026-10-23-browse-arxiv-raw/ep-1801.07667/ktiles.tex
using the paper's own edge macros (lines 232-236 of ktiles.tex):

    \r{i} = "to node {0} ++(60*i:1)"   -> unit edge, direction 60i deg, label 0
    \b{i} = "to node {1} ++(60*i:1)"   -> unit edge, direction 60i deg, label 1
    \bb{i}= "to node {2} ++(60*i:1)"   -> unit edge, direction 60i deg, label 2

Coordinates: \xcoord=cos60, \ycoord=sin60 (ktiles.tex:4-5), so the tikz
coordinate system IS the triangular-lattice basis
    e_0 = (1,0) (dir 0),  e_1 = (1/2, sqrt3/2) (dir 60).
We work in INTEGER (a,b) coords throughout; the shear is affine so
point-in-polygon is unaffected.
"""
from functools import lru_cache

# direction d  ->  step in (a,b) coords
D = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

# ---------------------------------------------------------------- tiles
# Each tile is a CCW list of (direction, label) unit steps, verbatim from tex.
TILES = {
    # KTW (fig LR-tiles, ktiles.tex:578-597)
    "zero_A":  [(0,0),(2,0),(4,0)],                       # \r0\r2\r4
    "zero_B":  [(1,0),(3,0),(5,0)],                       # \r1\r3\r5
    "one_A":   [(0,1),(2,1),(4,1)],                       # \b0\b2\b4
    "one_B":   [(1,1),(3,1),(5,1)],                       # \b1\b3\b5
    "rh_014":  [(0,1),(1,0),(3,1),(4,0)],                 # \b0\r1\b3\r4
    "rh_230":  [(2,1),(3,0),(5,1),(0,0)],                 # \b2\r3\b5\r0
    "rh_125":  [(1,1),(2,0),(4,1),(5,0)],                 # \b1\r2\b4\r5
    # the four K-theoretic pieces (fig K-tiles, ktiles.tex:768-792)
    "K":       [(1,1),(1,0),(3,1),(3,0),(5,1),(5,0)],     # \b1\r1\b3\r3\b5\r5
    "tri":     [(0,1),(0,0),(2,1),(2,0),(4,1),(4,0)],     # \b0\r0\b2\r2\b4\r4
    "hex":     [(1,0),(2,1),(3,0),(4,1),(5,0),(0,1)],     # \r1\b2\r3\b4\r5\b0
    "hexR":    [(0,0),(1,1),(2,0),(3,1),(4,0),(5,1)],     # \r0\b1\r2\b3\r4\b5
    # section ss:protrude complete set (ktiles.tex:1601-1625)
    "one_A2":  [(0,1),(2,1),(4,2)],                       # \b0\b2\bb4
    "one_B2":  [(1,2),(3,1),(5,1)],                       # \bb1\b3\b5
    "rh_125a": [(1,2),(2,0),(4,1),(5,0)],                 # \bb1\r2\b4\r5
    "rh_125b": [(1,2),(2,0),(4,2),(5,0)],                 # \bb1\r2\bb4\r5
    "hex2":    [(1,0),(2,1),(3,0),(4,2),(5,0),(0,1)],     # \r1\b2\r3\bb4\r5\b0
    "trap":    [(1,0),(2,1),(4,0),(4,1),(0,1)],           # \r1\b2\r4\b4\b0
}

def closes(steps):
    a = sum(D[d][0] for d, _ in steps)
    b = sum(D[d][1] for d, _ in steps)
    return (a, b) == (0, 0)

def ccw_area2(steps):
    """Twice the signed area in (a,b) coords (positive = CCW)."""
    p = (0, 0); s = 0
    for d, _ in steps:
        q = (p[0] + D[d][0], p[1] + D[d][1])
        s += p[0] * q[1] - q[0] * p[1]
        p = q
    return s

def verts(steps, anchor=(0, 0)):
    p = anchor; out = [p]
    for d, _ in steps:
        p = (p[0] + D[d][0], p[1] + D[d][1]); out.append(p)
    return out[:-1]

def edge_key(p, d):
    """Canonical undirected lattice edge."""
    if d < 3:
        return (p, d)
    q = (p[0] + D[d][0], p[1] + D[d][1])
    return (q, d - 3)

def tile_edges(steps, anchor):
    """[(edge_key, label)] for the tile's boundary unit edges."""
    p = anchor; out = []
    for d, lab in steps:
        out.append((edge_key(p, d), lab))
        p = (p[0] + D[d][0], p[1] + D[d][1])
    return out

# ----------------------------------------- unit triangles & covering
def inside(poly, pt3):
    """Ray-casting; poly in (a,b) ints, pt3 = 3*(point) to avoid vertex hits."""
    x, y = pt3; n = len(poly); c = False
    for i in range(n):
        x1, y1 = poly[i][0]*3, poly[i][1]*3
        x2, y2 = poly[(i+1) % n][0]*3, poly[(i+1) % n][1]*3
        if (y1 > y) != (y2 > y):
            xint = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < xint:
                c = not c
    return c

def tri_centroid3(tri):
    """3 * centroid of a unit triangle given as 3 (a,b) vertices -> ints."""
    return (sum(v[0] for v in tri), sum(v[1] for v in tri))

def up_tri(a, b):   return ((a,b), (a+1,b), (a,b+1))
def dn_tri(a, b):   return ((a+1,b), (a+1,b+1), (a,b+1))

def tile_cover(steps, anchor, n):
    """Unit triangles of the side-n region covered by the tile. None if it
    protrudes (covers area outside the region)."""
    poly = verts(steps, anchor)
    area2 = ccw_area2(steps)          # = (# unit triangles) since unit tri has area2=1
    got = []
    for a in range(-2, n + 2):
        for b in range(-2, n + 2):
            for tag, tri in (("U", up_tri(a,b)), ("D", dn_tri(a,b))):
                if inside(poly, tri_centroid3(tri)):
                    got.append((tag, a, b))
    if len(got) != area2:
        raise AssertionError("cover mismatch %d vs %d" % (len(got), area2))
    for tag, a, b in got:
        if a < 0 or b < 0: return None
        if tag == "U" and a + b > n - 1: return None
        if tag == "D" and a + b > n - 2: return None
    return got

def region_tris(n):
    out = []
    for a in range(n):
        for b in range(n):
            if a + b <= n - 1: out.append(("U", a, b))
            if a + b <= n - 2: out.append(("D", a, b))
    return out

# ----------------------------------------------- partitions <-> strings
def sigma(lam, n, k):
    """Binary string of partition lam in the k x (n-k) ambient rectangle:
    lower-right boundary read top-right to bottom-left, 0=horizontal,1=vertical.
    (ktiles.tex:526-536 and fig 'binary string')."""
    lam = list(lam) + [0] * (k - len(lam))
    assert len(lam) == k and all(0 <= x <= n - k for x in lam), (lam, n, k)
    assert all(lam[i] >= lam[i+1] for i in range(k-1)), lam
    s = []; col = n - k
    for i in range(k):
        s += [0] * (col - lam[i]); s += [1]; col = lam[i]
    s += [0] * col
    assert len(s) == n and sum(s) == k
    return s

def transpose(lam):
    lam = [x for x in lam if x > 0]
    if not lam: return []
    return [sum(1 for x in lam if x > j) for j in range(lam[0])]

# ----------------------------------------------- the solver
def boundary_labels(n, lam, mu, nu, k):
    """edge_key -> label on the three sides of the side-n triangle.
    bottom (dir0, (0,0)->(n,0))    = sigma(nu)
    left   (dir1, (0,0)->(0,n))    = sigma(lam)
    right  (dir5, (0,n)->(n,0))    = sigma(mu)
    All three verified against fig 'bdry' (ktiles.tex:610-617)."""
    bl = {}
    sn, sl, sm = sigma(nu, n, k), sigma(lam, n, k), sigma(mu, n, k)
    for i in range(n): bl[edge_key((i, 0), 0)] = sn[i]
    for j in range(n): bl[edge_key((0, j), 1)] = sl[j]
    for m in range(n): bl[edge_key((m, n - m), 5)] = sm[m]
    return bl

def count_puzzles(n, lam, mu, nu, k, tileset, limit=None, collect=False):
    """Exact-cover count of tilings of the side-n triangle by lattice
    translates of the given tiles, edge labels matching, boundary = Delta."""
    bl = boundary_labels(n, lam, mu, nu, k)
    # precompute placements
    placements = []           # (frozenset(tris), [(edge,label)])
    for name in tileset:
        steps = TILES[name]
        for a in range(n + 1):
            for b in range(n + 1):
                try:
                    cov = tile_cover(steps, (a, b), n)
                except AssertionError:
                    raise
                if cov is None: continue
                edges = tile_edges(steps, (a, b))
                ok = True
                for e, lab in edges:
                    if e in bl and bl[e] != lab: ok = False; break
                if ok: placements.append((name, (a,b), set(cov), edges))
    tris = region_tris(n)
    idx = {t: i for i, t in enumerate(tris)}
    by_tri = {t: [] for t in tris}
    for p in placements:
        first = min(p[2], key=lambda t: idx[t])
        by_tri[first].append(p)
    order = tris
    labels = {}          # edge -> label currently asserted
    used = set()
    sols = []
    cnt = 0
    def rec(i):
        nonlocal cnt
        if limit is not None and cnt >= limit: return
        while i < len(order) and order[i] in used: i += 1
        if i == len(order):
            cnt += 1
            if collect: sols.append(list(cur))
            return
        t = order[i]
        for name, anc, cov, edges in by_tri[t]:
            if cov & used: continue
            newe = []
            bad = False
            for e, lab in edges:
                if e in labels:
                    if labels[e] != lab: bad = True; break
                else: newe.append((e, lab))
            if bad: continue
            for e, lab in newe: labels[e] = lab
            used.update(cov); cur.append((name, anc))
            rec(i + 1)
            cur.pop(); used.difference_update(cov)
            for e, _ in newe: del labels[e]
    cur = []
    rec(0)
    return (cnt, sols) if collect else cnt
