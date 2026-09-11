"""
Check the PROOF of prop:reflect: for each (I)/(II) configuration, the partition
   B_lambda = B \ {b,c} u {b',e}          (same formula for both types)
is a partition of n, differs from B by exactly the two prescribed bead moves, and
the row lambda of M has EXACTLY TWO nonzero entries, with the predicted heights.
Verified against engine E2 (the matrix built from rim hooks), which knows nothing
about this argument.
"""
from graph import maya, partitions, vertices
from matrix_engine import matrix, contains, is_rim_hook
from itertools import combinations

def beads_to_partition(Bset, ell_pad):
    bs = sorted(Bset, reverse=True)
    lam = [bs[i] + (i+1) for i in range(len(bs))]
    return tuple(x for x in lam if x > 0)

def configs(mu):
    """all (I) and (II) quadruples; yields (type, b, c, bp, e, col1, col2, pred_exp)"""
    ell, bs, is_bead = maya(mu)
    n = sum(mu); N = ell + n + 5
    B = {(mu[i] if i < ell else 0) - (i+1) for i in range(N)}
    holes = [u for u in range(-ell, bs[0]+n+2) if u not in B]
    out = []
    for c, b in combinations(sorted(bs), 2):                     # (I)
        for bp in [u for u in holes if u < c]:
            e = b + c - bp
            if e not in B:
                h = sum(1 for v in range(c+1, b) if v in B)
                out.append(('I', b, c, bp, e, (b, bp), (c, bp), h))
    for b in bs:                                                  # (II)
        hb = [u for u in holes if u < b]
        for bp, e in combinations(hb, 2):
            c = bp + e - b
            if c in B:
                h = 1 + sum(1 for v in range(bp+1, e) if v in B)
                out.append(('II', b, c, bp, e, (b, bp), (b, e), h))
    return out, B, N

def pair_of_gamma(mu, gam, N):
    ell = len(mu); lg = len(gam)
    Bm = {(mu[i] if i < ell else 0)-(i+1) for i in range(N)}
    Bg = {(gam[i] if i < lg else 0)-(i+1) for i in range(N)}
    return ((Bm-Bg).pop(), (Bg-Bm).pop())

bad = 0; tot = 0; nlam = 0
for n in range(1, 9):
    for mu in partitions(n):
        cfg, B, N = configs(mu)
        rows, cols, M = matrix(mu)
        P = {g: pair_of_gamma(mu, g, N) for g in cols}
        for (ty, b, c, bp, e, col1, col2, hpred) in cfg:
            tot += 1
            Blam = (B - {b, c}) | {bp, e}
            lam = beads_to_partition(Blam, N)
            if sum(lam) != n: bad += 1; print('SIZE', mu, ty, lam); continue
            if lam == mu:     bad += 1; print('LAM=MU', mu, ty); continue
            nz = [(g, M[(lam, g)]) for g in cols if (lam, g) in M]
            if len(nz) != 2:
                bad += 1; print('NOT TWO', mu, ty, (b,c,bp,e), 'lam=',lam, 'nz=',nz); continue
            got = {P[g]: h for g, h in nz}
            if set(got) != {col1, col2}:
                bad += 1; print('WRONG COLS', mu, ty, set(got), {col1,col2}); continue
            if got[col1] - got[col2] != hpred:
                bad += 1; print('WRONG EXP', mu, ty, (b,c,bp,e), got[col1]-got[col2], hpred); continue
            nlam += 1
print(f'prop:reflect proof-witness: {nlam}/{tot} configurations (n<=8) give a row lambda != mu')
print(f'   with EXACTLY two nonzero entries, in the predicted columns, with the predicted')
print(f'   exponent difference.   failures = {bad}')
