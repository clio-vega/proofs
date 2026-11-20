"""PROVE 11-20: Claim U as a FINITE LOCAL CHECK (no boards).
For each journey type T, each allowed class u, each GOOD selected config G at R (xrule True, across in {True,None}),
let H' = image of G under its rotation, containing the new rhombus R'.  For each BAD config X at R' (class u', xrule True,
across False): X and H' are both sets of tiles of the same board, so every tile of X and every tile of H' must be EQUAL or
INTERIOR-DISJOINT.  Claim: compatible  <=>  X == H'.  Also: which edges of R are interior to each bad config (for k = 0)."""
import itertools, cmath
from collections import Counter
from fractions import Fraction as Fr
import localenum as LE
def poly(t):
    zs = [LE.phi(P) for P in t[1]]; c = sum(zs)/len(zs); return sorted(zs, key=lambda z: cmath.phase(z-c))
def interior_overlap(s, t):
    """separating axis test for convex polygons; touching (zero-area) counts as disjoint"""
    A, B = poly(s), poly(t)
    for P in (A, B):
        for i in range(len(P)):
            e = P[(i+1) % len(P)] - P[i]; n = complex(-e.imag, e.real)
            pa = [(z*n.conjugate()).real for z in A]; pb = [(z*n.conjugate()).real for z in B]
            if max(pa) <= min(pb) + 1e-9 or max(pb) <= min(pa) + 1e-9: return False
    return True
def shift(H, v): return frozenset((t[0], frozenset(LE.add(P, v) for P in t[1])) for t in H)
def canon_R(R):
    pts0 = sorted({P[:2] for P in R[1]}); pts1 = sorted({P[2:] for P in R[1]})
    u0 = LE.sub(pts0[1], pts0[0]); u1 = LE.sub(pts1[1], pts1[0])
    if u0 == u1: return u0, pts0[0] + pts1[0]
    return u0, pts0[0] + pts1[1]
CONF = {}
for T in 'BC':
    for u in LE.UNITS:
        if LE.ETA[T](u) == 0: continue
        R, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, R, typ, set(H))
            if r['xrule']: CONF.setdefault((T, u), []).append((typ, frozenset(H), r, R))
st = Counter()
for (T, u), lst in CONF.items():
    for typ, G, r, R in lst:
        if r['across'] is False: continue          # G must be selectable
        # image of G
        V = set().union(*[t[1] for t in G]); zs = LE.convex_hull([LE.phi(P) for P in V])
        corners = [P for P in V if LE.is_corner(LE.phi(P), zs)]
        c = tuple(sum(Fr(P[i]) for P in corners)/6 for i in range(4))
        Hn = frozenset((t[0], frozenset(LE.rot_about(P, c, r['k']) for P in t[1])) for t in G)
        Rn = next(t for t in Hn if t[0] == 'rh' and True)
        # the image of R is the rhombus of Hn (only one rhombus)
        un, base = canon_R(Rn)
        base = tuple(int(x) for x in base)
        for typ2, X, r2, R2 in CONF.get((T, un), []):
            if r2['across'] is not False: continue
            Xs = shift(X, base)
            ok = all(s == t or not interior_overlap(s, t) for s in Xs for t in Hn)
            st[(T, 'G', u, typ, 'k', r['k'], '-> bad X', un, typ2, 'compatible', ok, 'X == image(G)', Xs == Hn)] += 1
for k in sorted(st, key=str): print(st[k], k)
print('---- interior edges of R in each BAD config (canonical R=[0,v]x[0,v]); E_u={v}x[0,v] (1-edge), E_w=[0,v]x{0} (0-edge)')
for (T, u), lst in sorted(CONF.items()):
    for typ, X, r, R in lst:
        if r['across'] is not False: continue
        inner = LE.boundary_edges(X)[1]
        names = []
        for e in inner:
            if e <= R[1]:
                P, Q = sorted(e); d = LE.sub(Q, P)
                names.append(('1-edge at pi0=' + str(P[:2])) if d[:2] == (0, 0) else ('0-edge at pi1=' + str(P[2:])))
        print(T, 'u', u, typ, 'R edges interior to X:', sorted(names))
