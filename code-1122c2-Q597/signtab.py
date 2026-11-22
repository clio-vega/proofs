"""PROVE 11-22 c2 Step 0: is the role a sign?  Offline on rolesN.pkl (rolesdump.py).
Per diverging role: p, c0 (first Jo step whose removed tile A[p] reads in W), types, g = a0+b0 of the lexmin vertex of R_p
(Lemma phi, 11-21 c2), psi_T(R_p) (no-repeat tex: h0=-eta_T(P1), h1=eta_T(P0), psi=(min h0, min h1)).
Tabulates candidate signs against I's side; role-swap lesion; B/B control.   usage: python3 signtab.py SIZE"""
import sys, pickle
from collections import Counter, defaultdict
from fractions import Fraction as Fr
SIZE = int(sys.argv[1])
sys.argv = [sys.argv[0], str(SIZE), '1', '1']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
eta = {'B': lambda a, b: b, 'C': lambda a, b: a + b}
def lexmin(vs): return min(tuple(Fr(c) for c in v) for v in vs)
def g(R): m = lexmin(R); return m[0] + m[1]
def psi(R, T):
    e = eta[T]; return (min(-e(v[2], v[3]) for v in R), min(e(v[0], v[1]) for v in R))
def sgn(v): return (v > 0) - (v < 0)
def role(M0, A, Bq, Jo, T):
    mA, mB, p, q = middles(A, Bq)
    if len(mA) < 1 or len(mB) < 1: return None
    W = set(M0) - set(apply(M0, Jo)); rem = {}
    for j, s in enumerate(Jo):
        for u in s[1]: rem.setdefault(u, j)
    reads = [sorted({rem[u] for u in s[1] if u in W}) for s in A]
    fc = next((k for k, r in enumerate(reads) if r), None)
    TO, TA = typed(Jo), typed(A)
    Hp = {u for u in A[p][1] if u[0] != 'rh'}
    sheetp = not any(TO[rem[u]][0] == 1 and TO[rem[u]][1] in ('a0', 'a1') for u in Hp if u in W)
    R = rh(A[p][1])
    return dict(p=p, fc=fc, c0=reads[p][0] if reads[p] else None, creads=reads[p], tA=TA[p], tB=typed(Bq)[p],
                good=(sheetp and len(Hp & W) == 2), S=len(Hp & W), g=g(R), psi=psi(R, T), lenA=len(A), lenJ=len(Jo),
                gstart=g(rh(A[0][1])), mlen=len(mA))
st = Counter(); rows = []
for P in PAIRS:
    if P['kind'] == 'BC':
        sides = {k[0] for k in P['I']}
        side = 'empty' if not sides else ('two' if len(sides) > 1 else int(next(iter(sides))))
        Ts = ('B', 'C')
    else:
        side = 'BB'; Ts = ('B', 'B')
    r = [role(P['M0'], A, Bq, Jo, T) for (t, A, Bq, Jo), T in zip(P['roles'], Ts)]
    st[('PRE', P['kind'], 'both roles diverge', all(x is not None for x in r))] += 1
    if not all(r): continue
    st[('PRE', P['kind'], 'fc==p both', all(x['fc'] == x['p'] for x in r))] += 1
    st[('PRE', P['kind'], 'side', side)] += 1
    # mutual reading: role0's contact reads Jo=role1's journey at c0; is that role1's contact step p?
    st[('MUTUAL', P['kind'], side, 'c0_0 - p_1', r[0]['c0'] - r[1]['p'], 'c0_1 - p_0', r[1]['c0'] - r[0]['p'])] += 1
    rows.append((P, side, r))
    goodrole = tuple(x['good'] for x in r)
    st[('GOOD', P['kind'], side, goodrole)] += 1
    feats = {
        'dg = g0-g1': r[0]['g'] - r[1]['g'],
        'dp = p0-p1': r[0]['p'] - r[1]['p'],
        'dc0 = c0_0-c0_1': r[0]['c0'] - r[1]['c0'],
        'p0-c0_1': r[0]['p'] - r[1]['c0'],
        'dlen = lenA0-lenA1': r[0]['lenA'] - r[1]['lenA'],
        'dgstart': r[0]['gstart'] - r[1]['gstart'],
        'class0 of R_p (role0,role1)': (r[0]['tA'][0], r[1]['tA'][0]),
        'dpsi0 (psi_B role0 - psi_C role1, comp0)': r[0]['psi'][0] - r[1]['psi'][0],
        'dpsi1': r[0]['psi'][1] - r[1]['psi'][1],
        'dmlen': r[0]['mlen'] - r[1]['mlen'],
    }
    for f, v in feats.items():
        st[('FEAT', P['kind'], f, side, sgn(v) if isinstance(v, (int, Fr)) else v)] += 1
for k in sorted(st, key=str): print(st[k], k)
pickle.dump([(P['kind'], P['key'], P['x'], P['y'], side, r) for P, side, r in rows], open('/home/clio/projects/proofs/code-1122c2-Q597/signrows%d.pkl' % SIZE, 'wb'))
