"""Finite check (36 cases): for puzzle unit vectors e (0-edge) and f (1-edge), the Z^4 cell (e,0)x(0,f) under the engine's
phi is a square / 30-150 rhombus / segment; compare with: the puzzle parallelogram on e,f is a Purbhoo-chirality KT rhombus
(both triangles ccw-labelled in mosaic.OK, short diagonal labelled 10=2); and record the orientation sign of psi o phi^-1."""
import sys, cmath, math
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
G = MZ.Geo(-1)
U = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
def cz(v): return MZ.tri(float(v[0]), float(v[1]))
def det(a, b): return (a.conjugate() * b).imag
def add(p, q): return (p[0] + q[0], p[1] + q[1])
rows = []
for e in U:
    for f in U:
        pe = G.phi((e[0], e[1], 0, 0)); pf = G.phi((0, 0, f[0], f[1]))
        ang = round(math.degrees(cmath.phase(pf / pe))) % 360
        shape = 'square' if ang in (90, 270) else ('rhombus' if ang in (30, 150, 210, 330) else 'other')
        dpsi = det(cz(e), cz(f)); dphi = det(pe, pf)
        if abs(dpsi) < 1e-9: kt = 'degenerate'; orient = 'collapsed'
        else:
            orient = 'preserved' if dpsi * dphi > 0 else 'REVERSED'
            # short diagonal
            if add(e, (-f[0], -f[1])) in U: tris = [((0, 0), e, f), (e, add(e, f), f)]; diag = (e, f)
            else: tris = [((0, 0), e, add(e, f)), ((0, 0), add(e, f), f)]; diag = ((0, 0), add(e, f))
            def lab(p, q):
                d = (q[0] - p[0], q[1] - p[1])
                if {p, q} == set(diag): return 2
                if d in (e, (-e[0], -e[1])): return 0
                if d in (f, (-f[0], -f[1])): return 1
                raise ValueError
            good = True
            for t in tris:
                if det(cz(t[1]) - cz(t[0]), cz(t[2]) - cz(t[0])) < 0: t = (t[0], t[2], t[1])
                good &= tuple(lab(t[i], t[(i + 1) % 3]) for i in range(3)) in MZ.OK
            kt = 'OK-chirality' if good else 'wrong-chirality'
        rows.append((e, f, shape, kt, orient))
from collections import Counter
c = Counter((r[2], r[3], r[4]) for r in rows)
for k, v in sorted(c.items()): print(v, k)
# tiles of the other two forms: triangles tau x w and v x sigma: phi is a rotation by +-15 deg, psi the identity: preserved.
