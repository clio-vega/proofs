"""Common potentials.  For a rhombus R and target T put
   Pi_B(R) = (lam0, lam1, -mu_B1, -mu_B0 + lam1),   Pi_C(R) = (lam0, lam1, -(mu_C1 - lam0), -mu_C0)
(lam = (min a0, min a1), mu_T = psi_T(R)).  Claim 1: each component of Pi_T is non-decreasing along every good T-step
(read off the 12 configurations of the AW table, and on every journey of the census dumps).  Claim 2: Pi_B and Pi_C
shift by the SAME vector (D0, D2, -D1, D2+D3) under translation by D, so Pi_C(R') - Pi_B(R) depends only on R'-R and the
two shapes.  Per role (survivor s of Theorem local, target T, other target T'):
   e(s) = [Pi_T'(R^o_c0) - Pi_T(R_p)] - [Pi_T(R_p+1) - Pi_T(R_p)].
Failure of (Share) forces e(x) + e(y) >= 0 componentwise (sum of the two monotonicity inequalities).  Share forces
   e'(x) + e'(y) = 0, e' = Pi_T'(R^o_c0) - Pi_T(R_p)   (control)."""
import os, pickle, random
exec(open('/home/clio/projects/proofs/code-1123-Q597/pairloc.py').read().split("\nZ = (0, 0, 0, 0)\n")[0])
def Pi(R, T):
    l0, l1 = lam(R); m0, m1 = mu(R, T)
    return (l0, l1, -m1, -m0 + l1) if T == 'B' else (l0, l1, -(m1 - l0), -m0)
def sub(a, b): return tuple(p - q for p, q in zip(a, b))
def add(a, b): return tuple(p + q for p, q in zip(a, b))
st = Counter()
# Claim 2 on random translates of every rhombus shape in the table
random.seed(2)
for (T_, ck, typ), (Rc, rem, Rn, cre) in CF.items():
    for R in (Rc, Rn):
        for _ in range(30):
            D = tuple(random.randint(-6, 6) for _ in range(4)); sh = (D[0], D[2], -D[1], D[2] + D[3])
            for T in 'BC':
                assert Pi(tr(R, D), T) == add(Pi(R, T), sh)
st['claim2 shift (D0,D2,-D1,D2+D3) verified on table rhombi x 30 translates x 2 targets'] += 1
# Claim 1 on the table
for (T, ck, typ), (Rc, rem, Rn, cre) in CF.items():
    d = sub(Pi(Rn, T), Pi(Rc, T))
    st[('TABLE', T, cls(T, Rc), typ, 'dPi', d, 'nondecreasing', all(v >= 0 for v in d))] += 1
# Claim 1 on the census dumps (every journey of every pair, both sizes), and a lesion: the WRONG target's Pi
for SIZE in (5, 6):
    PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
    for P in PAIRS:
        for (t, A, Bq, Jo) in P['roles']:
            T = 'B' if (P['kind'] == 'BB' or t.endswith('B')) else 'C'; To = {'B': 'C', 'C': 'B'}[T] if P['kind'] == 'BC' else T
            for J, TT in ((A, T), (Bq, T), (Jo, To)):
                for s in J:
                    R0 = next(u for u in s[1] if u[0] == 'rh'); R1 = next(u for u in s[2] if u[0] == 'rh')
                    R0 = it(R0); R1 = it(R1)
                    d = sub(Pi(R1, TT), Pi(R0, TT))
                    st[('DATA size', SIZE, 'Pi_T nondecreasing', all(v >= 0 for v in d))] += 1
                    dl = sub(Pi(R1, {'B': 'C', 'C': 'B'}[TT]), Pi(R0, {'B': 'C', 'C': 'B'}[TT]))
                    st[('LESION wrong-target Pi, size', SIZE, 'nondecreasing', all(v >= 0 for v in dl))] += 1
# per-survivor excess
E = {}
for s in F:
    T = s['T']; To = OTHER[T]
    RA = s['RA']; R1 = step(T, RA, s['t'])[1]; Rc = s['seg'][s['c0']][0]
    ep = sub(Pi(Rc, To), Pi(RA, T)); e = sub(ep, sub(Pi(R1, T), Pi(RA, T)))
    E[id(s)] = (e, ep)
    st[('SURVIVOR', s['key'][:12], 'e', e, "e'", ep)] += 1
for sx in SB:
    for sy in SC:
        e = add(E[id(sx)][0], E[id(sy)][0]); ep = add(E[id(sx)][1], E[id(sy)][1])
        st[('PAIR fail needs e(x)+e(y)>=0:', all(v >= 0 for v in e), 'neg comps', tuple(i for i in range(4) if e[i] < 0))] += 1
        x = role(sx, 'B', (0, 0, 0, 0)); y0 = role(sy, 'C', (0, 0, 0, 0))
        D = sub(base(x['seg'][0][0]), base(y0['RA']))
        y = role(sy, 'C', D)
        cons = y['RA'] == x['seg'][0][0] and x['seg'][0][1] == sy['t'] and y['seg'][0][0] == x['RA'] and y['seg'][0][1] == sx['t']
        st[('PAIR share-consistent (seg_x[0] = y contact AND seg_y[0] = x contact)', cons, "e'(x)+e'(y)=0", ep == (0, 0, 0, 0))] += 1
for k in sorted(st, key=str): print(st[k], k)
