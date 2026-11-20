"""PROVE 11-20: unified tile potential, checked on (1) every GOOD abstract configuration, (2) every real journey step.
psi(t) = (min h0 over pi_1(t), min h1 over pi_0(t)) for a non-rhombus tile t; mu(R) = same for the rhombus.
(REM) every removed non-rh tile t of H:  psi(t) <= mu(R) componentwise.
(CRE) every created non-rh tile t' of H': psi_i(t') >= mu_i(R') + 1 for SOME i.
(MU)  mu(R') <= mu(R) componentwise.
Then no non-rh tile is created at step j and removed at step k > j:  psi_i(t) >= mu_i(j+1)+1 >= mu_i(k)+1 > mu_i(k)."""
import sys, pickle, glob, itertools
from collections import Counter
from fractions import Fraction as Fr
import localenum as LE
def psi(t, T):
    eta = LE.ETA[T]
    return (min(-eta(P[2:]) for P in t[1]), min(eta(P[:2]) for P in t[1]))
def check(T, H, Hn):
    R = next(t for t in H if t[0] == 'rh'); Rn = next(t for t in Hn if t[0] == 'rh')
    mu, mun = psi(R, T), psi(Rn, T)
    rem = all(psi(t, T)[0] <= mu[0] and psi(t, T)[1] <= mu[1] for t in H if t[0] != 'rh')
    cre = all(psi(t, T)[0] >= mun[0]+1 or psi(t, T)[1] >= mun[1]+1 for t in Hn if t[0] != 'rh')
    mono = mun[0] <= mu[0] and mun[1] <= mu[1]
    return rem, cre, mono
st = Counter()
# (1) abstract
for T in 'BC':
    for u in LE.UNITS:
        if LE.ETA[T](u) == 0: continue
        R, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, R, typ, set(H))
            if not r['xrule'] or r['k'] is None: continue
            V = set().union(*[t[1] for t in H]); zs = LE.convex_hull([LE.phi(P) for P in V])
            corners = [P for P in V if LE.is_corner(LE.phi(P), zs)]
            c = tuple(sum(Fr(P[i]) for P in corners)/6 for i in range(4))
            Hn = frozenset((t[0], frozenset(LE.rot_about(P, c, r['k']) for P in t[1])) for t in H)
            good = r['across'] is not False
            st[('abstract', T, 'good' if good else 'bad', typ, 'REM,CRE,MU', check(T, H, Hn))] += 1
# (2) real journeys, + lesion (reversed journeys)
if len(sys.argv) > 1:
    for fn in sorted(glob.glob('dumpj5-*.pkl')):
        D = pickle.load(open(fn, 'rb'))
        for fam in ('BC', 'BB2'):
            for TT, M0, J in D[fam]:
                for lesion in (False, True):
                    JJ = [(ty, n, g) for ty, g, n in reversed(J)] if lesion else J
                    res = [check(TT, g, n) for ty, g, n in JJ]
                    st[('real', fam, 'LESION' if lesion else 'journey', 'all steps REM,CRE,MU', tuple(all(x[i] for x in res) for i in range(3)))] += 1
for k in sorted(st, key=str): print(st[k], k)
