"""Step 0 (PROVE 11-24 c2).  For every distinct census mosaic M0 (tiles in Z^4 = Lambda0 x Lambda1):
 (a) psi(P) = Lambda0-part + Lambda1-part (the 'puzzle projection', Purbhoo 2.4 / KT) sends every tri and sq of M0 onto
     a unit triangle / unit 60-rhombus, and these tile the size-n triangle exactly once: M0 minus rhombi IS a lifted puzzle.
 (b) every rhombus is a product e x f with e, f parallel (same planar lattice direction), and that direction is a side
     direction of the puzzle triangle (a nest rhombus).
 (c) pi_k of every rhombus lies on the boundary of Delta_k (k = 0, 1): no nest rhombus over an interior edge.
 (d) side length of Delta_k vs number of k-letters on a side of the puzzle.
Lesion (control on the instrument): move one rhombus by an interior Lambda0 vector so it sits over an interior edge;
(c) must flag it.     usage: python3 nestcheck.py SIZE"""
import sys, pickle
from collections import Counter
SIZE = int(sys.argv[1])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))

def p0(P): return (P[0], P[1])
def p1(P): return (P[2], P[3])
def psi(P): return (P[0] + P[2], P[1] + P[3])
def pk(P, k): return p0(P) if k == 0 else p1(P)

def hull_triangle(pts):
    """smallest up-triangle {a>=a0, b>=b0, a+b<=s} containing pts (lattice coords a + b w6)."""
    a0 = min(p[0] for p in pts); b0 = min(p[1] for p in pts); s = max(p[0] + p[1] for p in pts)
    return a0, b0, s

def on_boundary(seg, tri):
    a0, b0, s = tri
    return all(p[0] == a0 for p in seg) or all(p[1] == b0 for p in seg) or all(p[0] + p[1] == s for p in seg)

def analyse(M0, st, tag):
    tris = [u for u in M0 if u[0] == 'tri']; sqs = [u for u in M0 if u[0] == 'sq']; rhs = [u for u in M0 if u[0] == 'rh']
    # (a) puzzle projection
    cells = Counter()
    for u in tris:
        cells[frozenset(psi(P) for P in u[1])] += 1
    for u in sqs:
        V = [psi(P) for P in u[1]]
        # a 60-rhombus = two unit triangles; find its short diagonal: the two vertices at distance 1 that are not edges...
        Vs = set(V); assert len(Vs) == 4
        # split: unit triangles inside the 4 points
        for a in Vs:
            for b in Vs:
                for c in Vs:
                    if len({a, b, c}) == 3 and a < b < c:
                        T = sorted([a, b, c]);
                        # unit triangle test: pairwise differences are unit lattice vectors
                        U = {(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)}
                        if all((q[0] - p[0], q[1] - p[1]) in U for p in T for q in T if p != q):
                            cells[frozenset(T)] += 1
    allp = {psi(P) for u in M0 for P in u[1]}
    a0, b0, s = hull_triangle(allp); n = s - a0 - b0
    want = set()
    for b in range(n):
        for a in range(n - b):
            want.add(frozenset({(a0 + a, b0 + b), (a0 + a + 1, b0 + b), (a0 + a, b0 + b + 1)}))
            if a + b <= n - 2: want.add(frozenset({(a0 + a + 1, b0 + b), (a0 + a + 1, b0 + b + 1), (a0 + a, b0 + b + 1)}))
    st[(tag, '(a) tri+sq tile the puzzle triangle exactly once', set(cells) == want and all(v == 1 for v in cells.values()))] += 1
    # (b) rhombi are parallel products along a side direction
    sides = {(0, 1), (1, -1), (-1, 0)}
    for u in rhs:
        V = sorted(u[1]); A = V[0]
        d0 = {(p0(P)[0] - p0(A)[0], p0(P)[1] - p0(A)[1]) for P in V} - {(0, 0)}
        d1 = {(p1(P)[0] - p1(A)[0], p1(P)[1] - p1(A)[1]) for P in V} - {(0, 0)}
        ok = len(d0) == 1 and len(d1) == 1
        if ok:
            e = next(iter(d0)); f = next(iter(d1))
            ok = (e == f or e == (-f[0], -f[1])) and (e in sides or (-e[0], -e[1]) in sides)
        st[(tag, '(b) rhombus = e x f, e || f || a side of the puzzle', ok)] += 1
    # (c) pi_k of rhombi on boundary of Delta_k;  (d) side lengths
    for k in (0, 1):
        ktris = [u for u in tris if len({pk(P, 1 - k) for P in u[1]}) == 1]   # Lambda_k-triangles
        img = {pk(P, k) for u in M0 for P in u[1]}
        D = hull_triangle(img); m = D[2] - D[0] - D[1]
        for u in rhs:
            st[(tag, 'k', k, '(c) rhombus image on boundary of Delta_k', on_boundary({pk(P, k) for P in u[1]}, D))] += 1
        st[(tag, 'k', k, '(d) #Lambda_k-tris == m^2', len(ktris) == m * m)] += 1
        st[(tag, 'k', k, 'side m', m, 'n', n)] += 1

st = Counter(); seen = set(); lesioned = 0
for P in PAIRS:
    key = frozenset(P['M0'])
    if key in seen: continue
    seen.add(key)
    M0 = set(P['M0']); analyse(M0, st, 'clean')
    # lesion: translate one rhombus by an interior direction of Lambda0 (perpendicular-ish to its side)
    rh = sorted((u for u in M0 if u[0] == 'rh'), key=str)
    if rh and lesioned < 50:
        u = rh[0]; V = sorted(u[1]); A = V[0]
        d0 = [(p0(Q)[0] - p0(A)[0], p0(Q)[1] - p0(A)[1]) for Q in V if p0(Q) != p0(A)][0]
        # rotate side direction by 60 degrees: (a,b) -> (-b, a+b): points into / out of the triangle
        for rot in [(-d0[1], d0[0] + d0[1]), (d0[0] + d0[1], -d0[0])]:
            v = (rot[0], rot[1], 0, 0)
            u2 = ('rh', frozenset(tuple(x + y for x, y in zip(Q, v)) for Q in u[1]))
            L = (M0 - {u}) | {u2}; st2 = Counter(); analyse(L, st2, 'lesion')
            st[('lesion moved-rhombus', 'k=0 any rhombus flagged off-boundary',
                st2[('lesion', 'k', 0, '(c) rhombus image on boundary of Delta_k', False)] > 0)] += 1
        lesioned += 1
print('SIZE', SIZE, 'distinct M0:', len(seen))
for kk in sorted(st, key=str): print(st[kk], kk)
