"""(FB) via the hexagon H = z0*Delta0 + z1*Delta1 (Minkowski sum; geom.py's phi, z0 = e^{i pi/12}, z1 = e^{-i pi/12}).
For every distinct census M0 (code-1123c2-Q597/nest{5,6}.pkl, which carries the nest frames X, u, w):
 (1) every vertex P of every M0 tile has pi_k(P) in Delta_k (k = 0, 1);
 (2) every nest frame has u = -(Lambda0 unit step), w = +(Lambda1 unit step);
 (3) the run X + R w (resp. X + R u) lies on a supporting line of H: phi(M0 vertices) all weakly on one side, the
     outside-frame cells (-1, j) (resp. (i, -1)), -3 <= i, j <= 15, weakly on the other side.
Lesion: swap z0 <-> z1 (the other chirality of the planar embedding): (3) must fail.   usage: python3 fbhex.py"""
import pickle, cmath, math
from collections import Counter
W6 = cmath.exp(1j * math.pi / 3)
def mk(z0, z1):
    return lambda P: z0 * (P[0] + P[1] * W6) + z1 * (P[2] + P[3] * W6)
def add(P, Q, k=1): return tuple(a + k * b for a, b in zip(P, Q))
def tri(pts):
    a0 = min(p[0] for p in pts); b0 = min(p[1] for p in pts); s = max(p[0] + p[1] for p in pts); return a0, b0, s
def inT(p, T): return p[0] >= T[0] and p[1] >= T[1] and p[0] + p[1] <= T[2]
def side(phi, X, d, Q):  # sign of cross(phi(X+d)-phi(X), phi(Q)-phi(X))
    a = phi(add(X, d)) - phi(X); b = phi(Q) - phi(X); return a.real * b.imag - a.imag * b.real
st = Counter(); seen = set()
for SIZE in (5, 6):
    for P in pickle.load(open('/home/clio/projects/proofs/code-1123c2-Q597/nest%d.pkl' % SIZE, 'rb')):
        key = frozenset(P['M0'])
        if key in seen: continue
        seen.add(key); M0 = P['M0']; V = {Q for t in M0 for Q in t[1]}
        # Delta_k from the Lambda_k-triangles (G1: they tile Delta_k)
        for k in (0, 1):
            img = {(Q[0], Q[1]) if k == 0 else (Q[2], Q[3]) for t in M0 if t[0] == 'tri'
                   and len({(R[2], R[3]) if k == 0 else (R[0], R[1]) for R in t[1]}) == 1 for Q in t[1]}
            T = tri(img)
            st[(SIZE, 'k', k, '(1) all M0 vertices have pi_k in Delta_k',
                all(inT((Q[0], Q[1]) if k == 0 else (Q[2], Q[3]), T) for Q in V))] += 1
        for name, (z0, z1) in (('clean', (cmath.exp(1j * math.pi / 12), cmath.exp(-1j * math.pi / 12))),
                               ('LESION swapped chirality', (cmath.exp(-1j * math.pi / 12), cmath.exp(1j * math.pi / 12)))):
            phi = mk(z0, z1)
            for idx, N in enumerate(P['nests']):
                X, u, w = N['X'], N['u'], N['w']
                if name == 'clean':
                    st[(SIZE, '(2) u = -Lambda0 step, w = +Lambda1 step', u[2:] == (0, 0) and w[:2] == (0, 0))] += 1
                for run, other, rows in ((w, u, 'cells (-1,j)'), (u, w, 'cells (i,-1)')):
                    sv = [side(phi, X, run, Q) for Q in V]
                    if rows == 'cells (-1,j)':
                        cells = [add(add(X, u, -1), w, j) for j in range(-3, 16)]
                    else:
                        cells = [add(add(X, w, -1), u, i) for i in range(-3, 16)]
                    cv = [side(phi, X, run, add(add(q, a), b)) for q in cells for a in ((0,) * 4, u) for b in ((0,) * 4, w)]
                    ok = (max(sv) <= 1e-9 and min(cv) >= -1e-9) or (min(sv) >= -1e-9 and max(cv) <= 1e-9)
                    st[(name, '(3) run is a supporting line separating H from', rows, ok)] += 1
print('distinct M0:', len(seen))
for k in sorted(st, key=str): print(st[k], k)
