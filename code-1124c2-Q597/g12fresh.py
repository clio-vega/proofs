"""Independent check of Lemma D / (G1) / (G2) on FRESH puzzles (not the census): every KT puzzle of size n (all boundaries),
built by mosaic.py.  For k in {0,1}:
 N(t) = #Lambda_k-triangles over unit triangle t: must be 1 on Delta_k (side m_k = #k-letters per side), 0 elsewhere;
 (G2) every interior edge E_i(x) of Delta_k: squares over E have Lambda_{1-k}-edges = unit chain h(t-) -> h(t+) step n_i,
      with n_i PREDICTED by the proof (k=0: class i-1, k=1: class i+1, sign: left of DIRS[i]); no nest rhombus over E.
Lesion: allow the MIRROR KT piece (ccw labels 10,1,0) as well; count lift failures and (G1)/(G2) failures among lifts.
usage: python3 g12fresh.py NMAX [mirror]"""
import sys
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
import mosaic as MZ
NMAX = int(sys.argv[1]); MIRROR = len(sys.argv) > 2
if MIRROR:
    for r in range(3):
        t = (2, 1, 0); MZ.OK.add(t[r:] + t[:r])
DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
CLS = {0: (1, 0), 1: (0, 1), 2: (-1, 1)}
def left(d, v):  # v strictly to the left of d (planar, a + b w6): cross product sign in real coords
    import math
    dx, dy = d[0] + d[1] / 2, d[1] * math.sqrt(3) / 2; vx, vy = v[0] + v[1] / 2, v[1] * math.sqrt(3) / 2
    return dx * vy - dy * vx > 1e-9
def predicted_n(k, i):
    c = CLS[(i - 1) % 3 if k == 0 else (i + 1) % 3]
    return c if left(DIRS[i], c) else (-c[0], -c[1])
def pk(P, k): return (P[0], P[1]) if k == 0 else (P[2], P[3])
def run(n, lab, st):
    try:
        M = MZ.build_mosaic(n, lab)
    except AssertionError:
        st['lift fails'] += 1; return
    st['lifted'] += 1
    tiles = M['tiles']; rh = [r for nest in M['nests'] for r in nest]
    for k in (0, 1):
        ktri = [t for t in tiles if t[0] == 'tri' and len({pk(P, 1 - k) for P in t[1]}) == 1]
        N = Counter(frozenset(pk(P, k) for P in t[1]) for t in ktri)
        H = {frozenset(pk(P, k) for P in t[1]): pk(t[1][0], 1 - k) for t in ktri}
        B = M['L'][(0, 0)]; o = pk(B, k)
        m = sum(1 for e in MZ.side_edges(n)[0] if lab[MZ.E(*e)] == k)
        want = set()
        for b in range(m):
            for a in range(m - b):
                want.add(frozenset({(o[0] + a, o[1] + b), (o[0] + a + 1, o[1] + b), (o[0] + a, o[1] + b + 1)}))
                if a + b <= m - 2: want.add(frozenset({(o[0] + a + 1, o[1] + b), (o[0] + a + 1, o[1] + b + 1), (o[0] + a, o[1] + b + 1)}))
        st[('k', k, 'G1: N=1 on Delta_k, 0 off', set(N) == want and all(v == 1 for v in N.values()))] += 1
        SQ = Counter(); SQE = {}
        for t in tiles:
            if t[0] != 'sq': continue
            E = frozenset(pk(P, k) for P in t[1])
            if len(E) == 2: SQE.setdefault(E, []).append(frozenset(pk(P, 1 - k) for P in t[1]))
        RH = Counter(frozenset(pk(P, k) for P in r[1]) for r in rh)
        bad = 0; tot = 0
        for tau in want:
            for x in tau:
                for i in range(3):
                    d = DIRS[i]; y = (x[0] + d[0], x[1] + d[1]); E = frozenset((x, y))
                    tm = frozenset((x, y, (x[0] + DIRS[i - 1][0], x[1] + DIRS[i - 1][1])))
                    tp = frozenset((x, y, (x[0] + DIRS[i + 1][0], x[1] + DIRS[i + 1][1])))
                    if tm not in want or tp not in want: continue
                    tot += 1
                    if RH[E]: bad += 1; continue
                    if tm not in H or tp not in H: bad += 1; continue
                    h = H[tm]; nn = predicted_n(k, i); K = len(SQE.get(E, []))
                    chain = {frozenset(((h[0] + r * nn[0], h[1] + r * nn[1]), (h[0] + (r + 1) * nn[0], h[1] + (r + 1) * nn[1]))) for r in range(K)}
                    ok = set(SQE.get(E, [])) == chain and H[tp] == (h[0] + K * nn[0], h[1] + K * nn[1])
                    bad += not ok
        st[('k', k, 'G2 predicted-chain all interior edges', bad == 0)] += 1
        st[('k', k, 'interior edges checked')] += tot // 1
st = Counter()
for n in range(1, NMAX + 1):
    for lab in MZ.puzzles(n, {}):
        st[('n', n)] += 1; run(n, lab, st)
for kk in sorted(st, key=str): print(st[kk], kk)
