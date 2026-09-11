"""Independent verification that w_N=(-1)^N gives multiplication by p_e.
Route: SSYT -> Kostka -> monomial basis -> invert Kostka.  Uses NO ribbon
combinatorics, no abacus, no Murnaghan-Nakayama."""
import sympy as sp
from itertools import product
from engine import partitions, apply_R

def ssyt_count(lam, alpha):
    """Kostka K_{lam,alpha}: # SSYT of shape lam, content alpha."""
    lam = list(lam)
    rows = len(lam)
    # fill entry by entry, row by row, enforcing weak-increase along rows,
    # strict increase down columns, and content alpha
    n = len(alpha)
    total = 0
    T = [[0]*lam[i] for i in range(rows)]
    def rec(i, j, remaining):
        nonlocal total
        if i == rows:
            total += 1 if all(r == 0 for r in remaining) else 0
            return
        if j == lam[i]:
            rec(i+1, 0, remaining); return
        lo = 1
        if j > 0: lo = max(lo, T[i][j-1])
        if i > 0: lo = max(lo, T[i-1][j] + 1)
        for v in range(lo, n+1):
            if remaining[v-1] > 0:
                T[i][j] = v
                remaining[v-1] -= 1
                rec(i, j+1, remaining)
                remaining[v-1] += 1
        T[i][j] = 0
    rec(0, 0, list(alpha))
    return total

def compositions(n, k):
    if k == 1:
        yield (n,); return
    for i in range(n+1):
        for rest in compositions(n-i, k-1):
            yield (i,) + rest

def schur_to_mon(lam, n, nvars):
    """dict: content-composition (as sorted partition) -> coefficient; s_lam in m-basis."""
    out = {}
    for mu in partitions(n):
        if len(mu) > nvars: continue
        out[mu] = ssyt_count(lam, tuple(mu) + (0,)*(nvars-len(mu)))
    return out

def mult_p_e(vec_m, e, n, nvars):
    """multiply an m-basis vector (partition->coeff) by p_e = sum x_i^e, return m-basis."""
    out = {}
    for mu, c in vec_m.items():
        if c == 0: continue
        # m_mu * p_e = sum over monomials x^alpha in m_mu of x^alpha * x_i^e
        # enumerate distinct rearrangements alpha of mu with <= nvars parts
        base = tuple(mu) + (0,)*(nvars-len(mu))
        seen = set()
        for perm in set(__import__('itertools').permutations(base)):
            if perm in seen: continue
            seen.add(perm)
            for i in range(nvars):
                new = list(perm); new[i] += e
                key = tuple(sorted([x for x in new if x > 0], reverse=True))
                out[key] = out.get(key, 0) + c
    # each partition key counted once per rearrangement -> divide by #rearrangements
    res = {}
    for key, v in out.items():
        r = len(set(__import__('itertools').permutations(tuple(key)+(0,)*(nvars-len(key)))))
        assert v % r == 0, (key, v, r)
        res[key] = v // r
    return res

def mon_to_schur(vec_m, n, nvars):
    """invert the Kostka matrix (unitriangular in dominance)."""
    parts = [p for p in partitions(n) if len(p) <= nvars]
    K = sp.Matrix([[ssyt_count(l, tuple(m)+(0,)*(nvars-len(m))) for m in parts] for l in parts])
    b = sp.Matrix([vec_m.get(m, 0) for m in parts])
    x = K.T.solve(b)         # s = sum_lam c_lam s_lam ; m-coeffs = K^T c
    return {parts[i]: x[i] for i in range(len(parts)) if x[i] != 0}

ok = bad = 0
for e in (1,2,3,4):
    for n in range(0, 6):
        for lam in partitions(n):
            nv = n + e
            vm = schur_to_mon(lam, n, nv)
            prod = mult_p_e(vm, e, n+e, nv)
            target = mon_to_schur(prod, n+e, nv)
            got = apply_R({tuple(lam): sp.Integer(1)}, e, lambda g,N: sp.Integer((-1)**N))
            got = {k: sp.Integer(v) for k, v in got.items()}
            if got == {k: sp.Integer(v) for k, v in target.items()}:
                ok += 1
            else:
                bad += 1
                if bad < 4: print("MN MISMATCH", lam, e, got, target)
print("R_e^{w_N=(-1)^N} s_lam == p_e * s_lam :  %d/%d pairs (e<=4, |lam|<=5)" % (ok, ok+bad))
