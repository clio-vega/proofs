"""Purbhoo mosaics and migration, built from KT puzzles via the R^4 lift.

Points are exact: P = (a0,b0,a1,b1) in Z^4 = (0-part, 1-part), each part a triangular-lattice vector
a + b*w6 (w6 = e^{i pi/3}).  Planar embedding phi(P) = z0*tri(a0,b0) + z1*tri(a1,b1), z0,z1 = e^{-/+ i pi/12}
(chirality chosen so that KT rhombus pieces become squares; asserted).
"""
import cmath, math, itertools
from fractions import Fraction as Fr

W6 = cmath.exp(1j*math.pi/3)
def tri(a, b): return a + b*W6

class Geo:
    def __init__(self, sign):
        self.z0 = cmath.exp(-sign*1j*math.pi/12); self.z1 = cmath.exp(sign*1j*math.pi/12)
    def phi(self, P): return self.z0*tri(float(P[0]), float(P[1])) + self.z1*tri(float(P[2]), float(P[3]))

def padd(P, Q): return tuple(p+q for p, q in zip(P, Q))
def psub(P, Q): return tuple(p-q for p, q in zip(P, Q))
def rot60(a, b): return (-b, a+b)          # multiplication by w6 on a+b*w6
def protk(P, k):
    a0, b0, a1, b1 = P
    for _ in range(k % 6):
        a0, b0 = rot60(a0, b0); a1, b1 = rot60(a1, b1)
    return (a0, b0, a1, b1)

# ---------------- KT puzzles (Purbhoo's chirality: ccw labels a rotation of (0,0,0),(1,1,1),(2,0,1), 2 = "10")
OK = set()
for t in [(0,0,0),(1,1,1),(2,0,1)]:
    for r in range(3): OK.add(t[r:]+t[:r])
def tri_cells(n):
    """triangles as ccw vertex triples in (a,b) lattice coords"""
    T = []
    for b in range(n):
        for a in range(n-b):
            T.append(((a,b),(a+1,b),(a,b+1)))
            if a+b <= n-2: T.append(((a+1,b),(a+1,b+1),(a,b+1)))
    return T
def E(p, q): return frozenset((p, q))
def side_edges(n):
    """boundary edges in clockwise traversal from lower-left: left side up, right side down, bottom leftwards"""
    left = [((0,k),(0,k+1)) for k in range(n)]
    right = [((k,n-k),(k+1,n-k-1)) for k in range(n)]
    bottom = [((n-k,0),(n-k-1,0)) for k in range(n)]
    return left, right, bottom
def puzzles(n, fixed, limit=None):
    """fixed: dict edge->label for (some) boundary edges.  Returns list of edge-label dicts."""
    T = tri_cells(n)
    # order triangles: BFS from cells touching most fixed edges (simple: sort by distance to the fixed set)
    fe = set(fixed)
    def score(t):
        es = [E(t[i], t[(i+1)%3]) for i in range(3)]
        return -sum(e in fe for e in es)
    order = []; known = set(fe); rest = list(T)
    while rest:
        rest.sort(key=lambda t: -sum(E(t[i], t[(i+1)%3]) in known for i in range(3)))
        t = rest.pop(0); order.append(t)
        for i in range(3): known.add(E(t[i], t[(i+1)%3]))
    bnd = set(E(*e) for s in side_edges(n) for e in s)
    out = []; lab = dict(fixed)
    def rec(k):
        if limit and len(out) >= limit: return
        if k == len(order): out.append(dict(lab)); return
        t = order[k]; es = [E(t[i], t[(i+1)%3]) for i in range(3)]
        for c in OK:
            new = []; good = True
            for e, l in zip(es, c):
                if e in lab:
                    if lab[e] != l: good = False; break
                else:
                    if l == 2 and e in bnd: good = False; break
                    new.append(e)
            if not good: continue
            for e, l in zip(es, c):
                if e not in lab: lab[e] = l
            rec(k+1)
            for e in new: del lab[e]
    rec(0); return out

# ---------------- mosaic from puzzle
class Mosaic:
    pass
def lift(n, lab):
    """vertex -> Z^4 lift; asserts consistency"""
    adj = {}
    for e, l in lab.items():
        if l == 2: continue
        p, q = tuple(e)
        adj.setdefault(p, []).append((q, l)); adj.setdefault(q, []).append((p, l))
    L = {(0,0): (0,0,0,0)}; st = [(0,0)]
    while st:
        p = st.pop()
        for q, l in adj[p]:
            d = (q[0]-p[0], q[1]-p[1])
            Q = padd(L[p], (d[0], d[1], 0, 0) if l == 0 else (0, 0, d[0], d[1]))
            if q in L: assert L[q] == Q, 'lift inconsistent'
            else: L[q] = Q; st.append(q)
    assert len(L) == (n+1)*(n+2)//2
    return L
def build_mosaic(n, lab, geo_sign=+1, corner_rule=0):
    L = lift(n, lab); tiles = []
    T = tri_cells(n); twos = {}
    for t in T:
        es = [E(t[i], t[(i+1)%3]) for i in range(3)]
        ls = [lab[e] for e in es]
        if 2 in ls:
            twos.setdefault(es[ls.index(2)], []).append(t)
        else:
            tiles.append(('tri', tuple(L[v] for v in t)))
    for e, ts in twos.items():
        assert len(ts) == 2
        t1, t2 = ts
        # quadrilateral ccw: start at vertex of t1 not on e, go ccw in t1 then the far vertex of t2
        i = [v for v in t1 if v not in e][0]; k = t1.index(i)
        v0, v1, v2 = t1[k], t1[(k+1)%3], t1[(k+2)%3]
        far = [v for v in t2 if v not in e][0]
        tiles.append(('sq', (L[v0], L[v1], L[far], L[v2])))
    # nests: per side, bubble the boundary string to the corner string
    sides = side_edges(n); corners = []; nests = []
    for s, edges in enumerate(sides):
        labels = [lab[E(*e)] for e in edges]
        vec = (edges[0][1][0]-edges[0][0][0], edges[0][1][1]-edges[0][0][1])
        def step(l): return (vec[0], vec[1], 0, 0) if l == 0 else (0, 0, vec[0], vec[1])
        start = L[edges[0][0]]
        target = sorted(labels, reverse=(corner_rule == 1))
        cur = list(labels); rh = []
        changed = True
        while changed:
            changed = False
            for k in range(n-1):
                if cur[k] != target[0] and cur[k+1] == target[0]:
                    p = start
                    for l in cur[:k]: p = padd(p, step(l))
                    x, y = cur[k], cur[k+1]
                    rh.append(('rh', (p, padd(p, step(x)), padd(padd(p, step(x)), step(y)), padd(p, step(y)))))
                    cur[k], cur[k+1] = y, x; changed = True
        assert cur == target
        p = start; cnt = sum(1 for l in target if l == target[0])
        for l in target[:cnt]: p = padd(p, step(l))
        corners.append(dict(Xp=start, X=p, Yp=L[edges[-1][1]],
                            s_first=step(target[0]), s_second=step(target[-1]), count_first=cnt))
        nests.append(rh)
    return dict(n=n, tiles=tiles, nests=nests, corners=corners, L=L, lab=lab)
