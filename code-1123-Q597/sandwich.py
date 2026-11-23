"""Tile sandwich.  tau(u) = (min a0, min a1, -min b0, min(a1+b1)) is target-free and shifts like Pi.
For every good configuration (T, class, type) of the AW table and every removed non-rhombus tile u, record
   lo = Pi_T(R_k) - tau(u),  hi = Pi_T(R_k+1) - tau(u)
by tile SHAPE (kind + vertex set translated to base 0).  If for every shape  max lo <= min hi  (over both targets), then any
two journeys (any targets) that remove a common tile u at steps k (first) and j (second) satisfy
   Pi_T2(R^2_j) <= Pi_T1(R^1_k+1)."""
exec(open('/home/clio/projects/proofs/code-1123-Q597/commonpot.py').read().split("st = Counter()")[0])
def tau(u):
    return (min(P[0] for P in u[1]), min(P[2] for P in u[1]), -min(P[1] for P in u[1]), min(P[2] + P[3] for P in u[1]))
def shape(u): b = base(u); return (u[0], tuple(sorted(tuple(a - c for a, c in zip(P, b)) for P in u[1])))
LO = {}; HI = {}; rows = []
for (T, ck, typ), (Rc, rem, Rn, cre) in CF.items():
    for u in rem:
        if u[0] == 'rh': continue
        lo = sub(Pi(Rc, T), tau(u)); hi = sub(Pi(Rn, T), tau(u)); sh = shape(u)
        LO.setdefault(sh, []).append(lo); HI.setdefault(sh, []).append(hi); rows.append((sh, T, cls(T, Rc), typ, lo, hi))
bad = 0
for sh in sorted(LO):
    mlo = tuple(max(l[i] for l in LO[sh]) for i in range(4)); mhi = tuple(min(h[i] for h in HI[sh]) for i in range(4))
    ok = all(a <= b for a, b in zip(mlo, mhi)); bad += not ok
    print(sh[0], 'n', len(LO[sh]), 'max lo', mlo, 'min hi', mhi, 'SANDWICH', ok)
print('shapes', len(LO), 'failing', bad)
for r in sorted(rows, key=str): print('ROW', r)
# which S-tiles does J_o's step c0 read, in each Theorem-local survivor; per square, the e that the sandwich gives
import pickle
for tag in ('', '-same'):
    for s in pickle.load(open('/home/clio/projects/proofs/code-1123-Q597/ccsurv%s.pkl' % tag, 'rb')):
        T = s['T']; To = T if tag else OTHER[T]; seg = s['seg']; c0 = s['c0']
        read = [u for u in seg[c0][2] if u in set(s['S'])]
        Rc = seg[c0][0]; R1 = step(T, s['RA'], s['t'])[1]
        e = sub(Pi(Rc, To), Pi(R1, T))
        print('SURV', tag or 'BC', T, (s['cA'], s['t']), 'c0 reads', sorted(u[0] for u in read), 'e', e)
