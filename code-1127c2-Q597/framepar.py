"""The three nest-frame parallelograms P_s = A'_s + [0,m0](D_s,0) + [0,m1](0,D_s) (Def M0 frames) have pairwise disjoint
interiors and lie in H, for all 1<=m0,m1<=15 (separating-axis test, exact up to 1e-9)."""
import sys; sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
G = MZ.Geo(-1)
Ds = [(0, 1), (1, -1), (-1, 0)]   # side directions of T_n (D tex Lemma D(2)); clockwise order as in side_edges
def P(v): return G.phi(v)
def para(Ap, D, m0, m1):
    a = (D[0], D[1], 0, 0); b = (0, 0, D[0], D[1])
    pts = [Ap, tuple(x + m0 * y for x, y in zip(Ap, a)), tuple(x + m0 * y + m1 * z for x, y, z in zip(Ap, a, b)), tuple(x + m1 * z for x, z in zip(Ap, b))]
    return [P(q) for q in pts], tuple(x + m0 * y + m1 * z for x, y, z in zip(Ap, a, b))
def sep(p, q):
    for poly in (p, q):
        for i in range(len(poly)):
            e = poly[(i + 1) % len(poly)] - poly[i]; n = complex(-e.imag, e.real)
            pr = [(z * n.conjugate()).real for z in p]; qr = [(z * n.conjugate()).real for z in q]
            if max(pr) <= min(qr) + 1e-9 or max(qr) <= min(pr) + 1e-9: return True
    return False
bad = 0; tot = 0
for m0 in range(1, 16):
    for m1 in range(1, 16):
        Ap = (0, 0, 0, 0); polys = []
        for D in Ds:
            poly, nxt = para(Ap, D, m0, m1); polys.append(poly); Ap = nxt
        closed = Ap == (0, 0, 0, 0)
        for i in range(3):
            for j in range(i + 1, 3):
                tot += 1; bad += (not sep(polys[i], polys[j])) or not closed
print('pairs', tot, 'overlapping or not closed', bad)
