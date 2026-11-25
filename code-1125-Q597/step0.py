"""PROVE 11-25 Step 0: positive control for the c2 mirror-lesion split (G1 582/582, G2 fails 277).
Family: every labelling of T_n, n<=4, by KT pieces AND mirror pieces (ccw rotation of (10,1,0)), that lifts (= g12fresh mirror).
For each lifted labelling lab, and for transforms
   id      : lab itself
   R       : reflect (a,b)->(b,a)                         (orientation reversing; KT piece -> mirror piece)
   RS      : reflect + swap labels 0<->1                  (KT -> KT)
compute (G1); (G2) with the KT normals (proof's n_i) and with the MIRRORED normals (class i+1 for k=0, i-1 for k=1);
(FB) for all three nest frames with z = (e^{i pi/12}, e^{-i pi/12}) and swapped z.
Chirality class of lab: 'KT' (all rhombi KT), 'mirror' (all rhombi mirror), 'mixed', 'none' (no rhombi).
Prediction: pure-KT: G2-KT ok, FB-z ok.  pure-mirror: G2-mirror ok (and G2-KT fails iff some interior column has a square).
Combined RS on pure-KT: still pure KT, all three checks pass (0 failures).  usage: python3 step0.py NMAX"""
import sys, cmath, math
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
KT = set(MZ.OK); MIR = set()
for r in range(3):
    t = (2, 1, 0); MIR.add(t[r:] + t[:r])
MZ.OK |= MIR                       # enumerate the mirror family (as g12fresh mirror)
NMAX = int(sys.argv[1])
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
CLS = {0: (1, 0), 1: (0, 1), 2: (-1, 1)}
def left(d, v):
    dx, dy = d[0] + d[1] / 2, d[1] * math.sqrt(3) / 2; vx, vy = v[0] + v[1] / 2, v[1] * math.sqrt(3) / 2
    return dx * vy - dy * vx > 1e-9
def normal(k, i, mirrored):
    sh = (-1 if k == 0 else 1) * (-1 if mirrored else 1)
    c = CLS[(i + sh) % 3]
    return c if left(DIRS[i], c) else (-c[0], -c[1])
def pk(P, k): return (P[0], P[1]) if k == 0 else (P[2], P[3])
def chir(n, lab):
    kinds = set()
    for t in MZ.tri_cells(n):
        ls = tuple(lab[MZ.E(t[i], t[(i + 1) % 3])] for i in range(3))
        if 2 in ls: kinds.add('KT' if ls in KT else 'mirror')
    return 'none' if not kinds else (kinds.pop() if len(kinds) == 1 else 'mixed')
def g12(n, lab, M):
    tiles = M['tiles']; rh = [r for nest in M['nests'] for r in nest]; out = {}
    for k in (0, 1):
        ktri = [t for t in tiles if t[0] == 'tri' and len({pk(P, 1 - k) for P in t[1]}) == 1]
        N = Counter(frozenset(pk(P, k) for P in t[1]) for t in ktri)
        H = {frozenset(pk(P, k) for P in t[1]): pk(t[1][0], 1 - k) for t in ktri}
        o = pk(M['L'][(0, 0)], k)
        m = sum(1 for e in MZ.side_edges(n)[0] if lab[MZ.E(*e)] == k)
        want = set()
        for b in range(m):
            for a in range(m - b):
                want.add(frozenset({(o[0] + a, o[1] + b), (o[0] + a + 1, o[1] + b), (o[0] + a, o[1] + b + 1)}))
                if a + b <= m - 2: want.add(frozenset({(o[0] + a + 1, o[1] + b), (o[0] + a + 1, o[1] + b + 1), (o[0] + a, o[1] + b + 1)}))
        out[('G1', k)] = set(N) == want and all(v == 1 for v in N.values())
        SQE = {}
        for t in tiles:
            if t[0] != 'sq': continue
            E = frozenset(pk(P, k) for P in t[1])
            if len(E) == 2: SQE.setdefault(E, []).append(frozenset(pk(P, 1 - k) for P in t[1]))
        RH = Counter(frozenset(pk(P, k) for P in r[1]) for r in rh)
        for mirrored in (False, True):
            bad = 0
            for tau in want:
                for x in tau:
                    for i in range(3):
                        d = DIRS[i]; y = (x[0] + d[0], x[1] + d[1]); E = frozenset((x, y))
                        tm = frozenset((x, y, (x[0] + DIRS[i - 1][0], x[1] + DIRS[i - 1][1])))
                        tp = frozenset((x, y, (x[0] + DIRS[i + 1][0], x[1] + DIRS[i + 1][1])))
                        if tm not in want or tp not in want: continue
                        if RH[E] or tm not in H or tp not in H: bad += 1; continue
                        h = H[tm]; nn = normal(k, i, mirrored); K = len(SQE.get(E, []))
                        chain = {frozenset(((h[0] + r * nn[0], h[1] + r * nn[1]), (h[0] + (r + 1) * nn[0], h[1] + (r + 1) * nn[1]))) for r in range(K)}
                        bad += not (set(SQE.get(E, [])) == chain and H[tp] == (h[0] + K * nn[0], h[1] + K * nn[1]))
            out[('G2', 'mirrored' if mirrored else 'KT', k)] = bad == 0
    return out
W6 = cmath.exp(1j * math.pi / 3)
def mkphi(z0, z1): return lambda P: z0 * (P[0] + P[1] * W6) + z1 * (P[2] + P[3] * W6)
def add(P, Q, c=1): return tuple(a + c * b for a, b in zip(P, Q))
def side(phi, X, d, Q):
    a = phi(add(X, d)) - phi(X); b = phi(Q) - phi(X); return a.real * b.imag - a.imag * b.real
def fb(M, z0, z1):
    phi = mkphi(z0, z1)
    V = {Q for t in M['tiles'] for Q in t[1]} | {Q for nest in M['nests'] for r in nest for Q in r[1]}
    for c in M['corners']:
        X = c['X']; u = tuple(-a for a in c['s_first']); w = c['s_second']
        for run, rows in ((w, 'c'), (u, 'r')):
            sv = [side(phi, X, run, Q) for Q in V]
            cells = [add(add(X, u, -1), w, j) for j in range(-3, 8)] if rows == 'c' else [add(add(X, w, -1), u, i) for i in range(-3, 8)]
            cv = [side(phi, X, run, add(add(q, a), b)) for q in cells for a in ((0,) * 4, u) for b in ((0,) * 4, w)]
            if not ((max(sv) <= 1e-9 and min(cv) >= -1e-9) or (min(sv) >= -1e-9 and max(cv) <= 1e-9)): return False
    return True
ZP = (cmath.exp(1j * math.pi / 12), cmath.exp(-1j * math.pi / 12))
def transform(n, lab, name):
    if name == 'id': return lab
    sw = (lambda l: l) if name == 'R' else (lambda l: {0: 1, 1: 0, 2: 2}[l])
    return {frozenset((p[1], p[0]) for p in e): sw(l) for e, l in lab.items()}
st = Counter()
for n in range(1, NMAX + 1):
    for lab in MZ.puzzles(n, {}):
        try: MZ.build_mosaic(n, lab)
        except AssertionError: continue
        c0 = chir(n, lab)
        for name in ('id', 'R', 'RS'):
            lab2 = transform(n, lab, name)
            try: M = MZ.build_mosaic(n, lab2)
            except AssertionError: st[(c0, name, 'no lift')] += 1; continue
            c2 = chir(n, lab2)
            r = g12(n, lab2, M)
            m0 = sum(1 for e in MZ.side_edges(n)[0] if lab2[MZ.E(*e)] == 0)
            key = (c0, name, 'image chirality', c2)
            st[key + ('count',)] += 1
            st[key + ('G1 ok',)] += all(r[('G1', k)] for k in (0, 1))
            st[key + ('G2 KT-normals ok',)] += all(r[('G2', 'KT', k)] for k in (0, 1))
            st[key + ('G2 mirrored-normals ok',)] += all(r[('G2', 'mirrored', k)] for k in (0, 1))
            if 0 < m0 < n:
                st[key + ('FB-eligible',)] += 1
                st[key + ('FB z ok',)] += fb(M, *ZP)
                st[key + ('FB z-swapped ok',)] += fb(M, ZP[1], ZP[0])
for kk in sorted(st, key=str): print(st[kk], kk)
