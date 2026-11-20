"""PROVE 11-20 c2: the cross-class past. For every step j, removed non-rh u, ascending edge vector v of R_j, t = u+v,
and every k < j: is class(R_k) = class(R_j)?  For cross-class k: is psi(t) <= mu_k (psi fails to exclude) and is
om(t, R_k) in Q (4-offset fails to exclude)?  Also the number of b-steps in [k, j).   Restricted to k with t a board tile BEFORE step k (else removal at k is impossible).  usage: python3 cross2.py"""
import pickle, glob
from collections import Counter
from psilib import psi, rh, tr, edgevecs, ETA
from offsets import om
QS = {'B': {(-1,0,0,0),(0,-1,0,0),(0,0,0,0),(0,0,1,0),(-1,0,1,0)}, 'C': {(-1,0,0,0),(0,-1,0,0),(0,0,0,0),(0,0,0,1),(0,-1,0,1)}}
def cls(R):
    V = sorted(R[1]); d = set()
    for P in V:
        for Q in V:
            w = tuple(b-a for a, b in zip(P, Q))
            if any(w[:2]) and not any(w[2:]): d.add(frozenset({w[:2], tuple(-x for x in w[:2])}))
    return frozenset(d)
st = Counter()
for fn in sorted(glob.glob('../code-1120-Q597/dumpj5-*.pkl')):
    D = pickle.load(open(fn, 'rb'))
    for fam in ('BC',):
        for TT, M0, J in D[fam]:
            boards = [set(M0)]
            for _, g, n in J: boards.append((boards[-1] - set(g)) | set(n))
            Rs = [rh(g) for _, g, _ in J]; C = [cls(R) for R in Rs]; mu = [psi(R, TT) for R in Rs]
            for j, (ty, g, n) in enumerate(J):
                for u in g:
                    if u[0] == 'rh': continue
                    pu = psi(u, TT)
                    for v in edgevecs(Rs[j]):
                        c = 1 if any(v[:2]) else 0
                        t = tr(u, v); pt = psi(t, TT)
                        if pt[c] - pu[c] != 1: continue
                        for k in range(j):
                            if t not in boards[k]: continue
                            same = C[k] == C[j]
                            nb = sum(1 for s in J[k:j] if s[0] == 'b')
                            psile = pt[0] <= mu[k][0] and pt[1] <= mu[k][1]
                            inq = om(t, Rs[k], TT) in QS[TT]
                            st[('same class' if same else 'CROSS class', 'b-steps in [k,j) parity', nb % 2, 'psi<=mu_k', psile, 'om in Q', inq)] += 1
for k in sorted(st, key=str): print(st[k], k)
