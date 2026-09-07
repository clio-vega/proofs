"""Q96 engine A: abacus / Maya-set representation.

Partitions as tuples (weakly decreasing positive ints).
Operators are dicts {lambda: coeff} -> {mu: coeff}, coeffs in QQ(t) via sympy.
"""
import sympy as sp
from functools import lru_cache

t = sp.Symbol('t')

# ---------- abacus ----------
def beta(lam, n):
    """beta-numbers with n beads: lam_j + n - j, j=1..n (lam padded with 0)."""
    lam = list(lam) + [0]*(n-len(lam))
    return tuple(lam[j] + n - 1 - j for j in range(n))   # j=0..n-1 -> lam_{j+1} + n - (j+1)

def from_beta(bs, n):
    """inverse: sorted descending beta -> partition."""
    bs = sorted(bs, reverse=True)
    lam = [bs[j] - (n - 1 - j) for j in range(n)]
    lam = [x for x in lam if x > 0]
    return tuple(lam)

def R_abacus(lam, e, tt, nbeads=None):
    """R_e(tt) s_lam  = sum over legal e-moves, weight tt^(#beads strictly between)."""
    n = (nbeads if nbeads is not None else len(lam) + sum(lam) + e + 2)
    bs = set(beta(lam, n))
    out = {}
    for b in list(bs):
        if b + e in bs:
            continue
        ht = sum(1 for x in bs if b < x < b + e)
        mu = from_beta((bs - {b}) | {b + e}, n)
        out[mu] = out.get(mu, 0) + tt**ht
    return out

# ---------- engine B: direct border strips on the Young diagram ----------
def cells(lam):
    return set((i, j) for i, r in enumerate(lam) for j in range(r))

def is_partition(rows):
    rows = [r for r in rows if r > 0]
    return all(rows[i] >= rows[i+1] for i in range(len(rows)-1))

def R_borderstrip(lam, e, tt):
    """Brute force: enumerate all mu with |mu|=|lam|+e, mu contains lam,
    mu/lam a connected border strip (skew shape, connected, no 2x2)."""
    out = {}
    L = list(lam)
    maxrows = len(L) + e
    # enumerate mu by adding a composition of e to the rows
    def rec(i, remaining, cur):
        if remaining == 0:
            full = cur + L[i:]
            mu = tuple(x for x in full if x > 0)
            if is_partition(list(full)):
                yield mu
            return
        if i >= maxrows:
            return
        base = L[i] if i < len(L) else 0
        for add in range(0, remaining + 1):
            yield from rec(i + 1, remaining - add, cur + [base + add])
    for mu in sorted(set(rec(0, e, []))):
        if len(mu) < len(lam):
            continue
        if any((mu[i] if i < len(mu) else 0) < (lam[i] if i < len(lam) else 0) for i in range(max(len(mu), len(lam)))):
            continue
        sk = cells(mu) - cells(lam)
        if len(sk) != e:
            continue
        # no 2x2
        if any((i,j) in sk and (i+1,j) in sk and (i,j+1) in sk and (i+1,j+1) in sk for (i,j) in sk):
            continue
        # connected (edge adjacency)
        start = next(iter(sk)); seen = {start}; stack = [start]
        while stack:
            (i,j) = stack.pop()
            for (a,b) in ((i+1,j),(i-1,j),(i,j+1),(i,j-1)):
                if (a,b) in sk and (a,b) not in seen:
                    seen.add((a,b)); stack.append((a,b))
        if len(seen) != e:
            continue
        rows = len(set(i for (i,j) in sk))
        out[mu] = out.get(mu, 0) + tt**(rows - 1)
    return out

# ---------- operator algebra on formal sums ----------
def apply(op, vec):
    """op: callable lam -> dict; vec: dict lam->coeff"""
    out = {}
    for lam, c in vec.items():
        for mu, d in op(lam).items():
            out[mu] = sp.expand(out.get(mu, 0) + c*d)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}

def comm(op1, op2):
    """returns callable lam -> dict for [op1,op2]"""
    def f(lam):
        a = apply(op1, op2(lam))
        b = apply(op2, op1(lam))
        out = dict(a)
        for k, v in b.items():
            out[k] = sp.expand(out.get(k, 0) - v)
        return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}
    return f

def Rop(e, tt):
    return lambda lam: R_abacus(lam, e, tt)

def Rop_bs(e, tt):
    return lambda lam: R_borderstrip(lam, e, tt)
