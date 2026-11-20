"""PROVE 11-20: the reflection Sigma(P0,P1) = (zbar c(P1), zbar c(P0)), c(a,b)=(a+b,-b) [complex conjugation on Z[zeta]],
zbar = rotation by -60 deg.  Claims: (1) phi(Sigma P) = e^{-i60} conj(phi P) (an isometric reflection, swapping Lambda_0/Lambda_1);
(2) x_C(Sigma P) = x_B(P); (3) h0^C(Sigma P) = h1^B(P), h1^C(Sigma P) = h0^B(P);
(4) Sigma maps the set of B-configurations (with all analysis data) onto the C-configurations, swapping a0<->a1 and the MU components."""
import cmath, math
from collections import Counter
import localenum as LE
def c(v): return (v[0]+v[1], -v[1])
def zbar(v): return (v[0]+v[1], -v[0])
def Sig(P): return zbar(c(P[2:])) + zbar(c(P[:2]))
def SigT(t):
    V = frozenset(Sig(P) for P in t[1]); return (t[0], V)
ok = Counter()
import itertools
for P in itertools.product(range(-2, 3), repeat=4):
    ok['phi'] += abs(LE.phi(Sig(P)) - cmath.exp(-1j*math.pi/3)*LE.phi(P).conjugate()) < 1e-9
    ok['x'] += abs(LE.x(Sig(P), 'C') - LE.x(P, 'B')) < 1e-9
    ok['h0'] += (-LE.ETA['C'](Sig(P)[2:])) == LE.ETA['B'](P[:2])
    ok['h1'] += LE.ETA['C'](Sig(P)[:2]) == -LE.ETA['B'](P[2:])
print('pointwise checks over 625 points:', dict(ok))
# configurations
def confs(T):
    out = {}
    for u in LE.UNITS:
        R, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, R, typ, set(H)); out[H] = (typ, r)
    return out
CB, CC = confs('B'), confs('C')
st = Counter()
for H, (typ, r) in CB.items():
    HS = frozenset(SigT(t) for t in H)
    # translate so that the image rhombus sits at canonical position: compare up to translation by matching any C config
    R = next(t for t in HS if t[0] == 'rh'); pts0 = sorted({P[:2] for P in R[1]}); pts1 = sorted({P[2:] for P in R[1]})
    base = pts0[0] + pts1[0]
    Hc = frozenset((t[0], frozenset(tuple(a-b for a, b in zip(P, base)) for P in t[1])) for t in HS)
    if Hc not in CC: st['image not found'] += 1; continue
    t2, r2 = CC[Hc]
    sw = {'a0': 'a1', 'a1': 'a0', 'b': 'b'}[typ]
    same = (t2 == sw and r2['xrule'] == r['xrule'] and r2['across'] == r['across'] and r2.get('DL') == r.get('DL')
            and (r.get('MU') is None or r2.get('MU') == (r['MU'][1], r['MU'][0])) and ({2: 4, 4: 2, 3: 3, None: None}[r['k']] == r2['k']))
    st[('image matches with a0<->a1, MU swapped, k 2<->4', same)] += 1
print(dict(st), 'B configs', len(CB), 'C configs', len(CC))
