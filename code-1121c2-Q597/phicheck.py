"""PROVE 11-21 c2: engine-side check of the facts behind Lemma phi, on the 11-20 dumps (3881 B/C journeys, 898 B/B-second),
independent of the configuration table: g(R) = a0 + b0 of the lex-min vertex.
 (i)  Delta g >= 0 at every step toward B, <= 0 toward C;  Delta g at class-0 b-steps: +1 (B), 0 (C);
 (ii) off(t) = phi(lexmin t) - g(R_k) for removed non-rhombus tiles: {0} for B, subset of {-1,0} for C;
 (iii) FIN: the final rhombus of a journey never has the nest-A class (Lambda_0 edge direction +-(0,1)).
Lesion: the same with phi replaced by eta_B on Lambda_0 (b0 alone) must fail (i) for some target.   usage: python3 phicheck.py"""
import pickle, glob
from fractions import Fraction as Fr
from collections import Counter
def lexmin(t): return min(tuple(Fr(x) for x in P) for P in t[1])
def rh(S): return next(u for u in S if u[0] == 'rh')
def l0dir(R):
    vs = [tuple(Fr(x) for x in P) for P in R[1]]
    for P in vs:
        for Q in vs:
            d = tuple(b - a for a, b in zip(P, Q))
            if any(d[:2]) and not any(d[2:]): return d[:2] if d[:2] > (0, 0) else tuple(-x for x in d[:2])
PHI = {'phi': lambda P: P[0] + P[1], 'LESION b0': lambda P: P[1]}
st = Counter()
for fn in sorted(glob.glob('/home/clio/projects/proofs/code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for grp in ('BC', 'BB2'):
        for TT, M0, J in D[grp]:
            for name, f in PHI.items():
                for s in J:
                    R = rh(s[1]); Rn = rh(s[2]); dg = f(lexmin(Rn)) - f(lexmin(R))
                    cls0 = l0dir(R) == (0, 1)
                    st[(grp, TT, name, 'Delta g', str(dg))] += 1
                    if cls0 and s[0] == 'b': st[(grp, TT, name, 'class-0 b-step Delta g', str(dg))] += 1
                    if name == 'phi':
                        for t in s[1]:
                            if t[0] != 'rh': st[(grp, TT, 'off(t)', str(f(lexmin(t)) - f(lexmin(R))))] += 1
            st[(grp, TT, 'FIN final rhombus nest-A class', l0dir(rh(J[-1][2])) == (0, 1))] += 1
            st[(grp, TT, 'journeys')] += 1
for k in sorted(st, key=str): print(st[k], k)
