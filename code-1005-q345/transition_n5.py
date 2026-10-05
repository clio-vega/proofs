"""Closes gap (1) of the paper for the converse direction: the transition
matrices at n=5.  Two separate claims, reported separately:
  (a) I_mat is unimodular over Z  -- the evidence for h_alpha being a Z-BASIS,
      which is the one step resting on a cited (not re-derived) classical fact;
  (b) eta = N . I^{-1} is integral, and eta computed ONLY from u=e reproduces
      N_{w/u}(p_beta) for all (u,w) -- universality.
(b) is run on a bounded sample at n=5, and the sample size is REPORTED, not hidden.
"""
from fractions import Fraction
from itertools import permutations, combinations_with_replacement
from schub import length
from chains import T_tensor
from ls_side import I_table
from transition import comps_le_delta, word, inv_matrix

n = 5
perms = list(permutations(range(1, n+1)))
e = tuple(range(1, n+1))
byl = {}
for v in perms: byl.setdefault(length(v), []).append(v)
maxm = n*(n-1)//2

for m in range(1, maxm+1):
    A = comps_le_delta(n, m); V = byl[m]
    assert len(A) == len(V), (m, len(A), len(V))
    Ie = {v: I_table(e, v, n) for v in V}
    I_mat = [[Ie[v].get(al, 0) for v in V] for al in A]
    Iinv = inv_matrix(I_mat)
    unimod = Iinv is not None and all(x.denominator == 1 for r in Iinv for x in r)
    BETAS = []
    for c in combinations_with_replacement(range(1, n), m):
        b = [0]*(n-1)
        for k in c: b[k-1] += 1
        BETAS.append(tuple(b))
    Te = {v: T_tensor(e, v, n) for v in V}
    N_mat = [[Te[v].get(word(be), 0) for v in V] for be in BETAS]
    eta = [[sum(Fraction(N_mat[i][k])*Iinv[k][j] for k in range(len(V)))
            for j in range(len(A))] for i in range(len(BETAS))]
    eta_int = all(x.denominator == 1 for r in eta for x in r)
    nz = sum(1 for r in eta for x in r if x) 
    multi = sum(1 for r in eta if sum(1 for x in r if x) > 1)
    # universality on a bounded sample
    bad = tot = 0; pairs = 0
    for u in perms:
        lu = length(u)
        if lu+m > maxm: continue
        for w in byl.get(lu+m, []):
            pairs += 1
            if pairs % 7:                     # 1-in-7 sample, stated
                continue
            T = T_tensor(u, w, n); I = I_table(u, w, n)
            if not T and not I: continue
            for i, be in enumerate(BETAS):
                lhs = T.get(word(be), 0)
                rhs = sum(eta[i][j]*I.get(A[j], 0) for j in range(len(A)))
                tot += 1
                if lhs != rhs:
                    bad += 1
                    if bad == 1: print("    UNIVERSALITY FAILS", u, w, be, lhs, rhs)
    print(f"n=5 m={m}: dim={len(A)}, |beta|={len(BETAS)}; I_mat unimodular over Z: {unimod}; "
          f"eta integral: {eta_int}; {multi} of {len(BETAS)} eta rows are multi-term; "
          f"universality {tot-bad}/{tot} on a 1-in-7 sample of the {pairs} eligible (u,w) pairs, {bad} bad",
          flush=True)
