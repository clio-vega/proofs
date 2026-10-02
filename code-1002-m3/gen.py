"""General-m slice-sum machinery for cylindric Kostka numbers k(a,b)=K^c_{lam/mu,(a,b,d-a-b)}.

Everything here follows [c1, Prop 2.3 / Prop "reform"]:
    Box(mu) = prod_i [mu_i, mu_{i+1}-1],   mu_{m+1}=mu_1+n
    L_i(nu) = max(nu_i, lam_{i-1}+1),      lam_0 = lam_m - n
    R_i(nu) = min(nu_{i+1}-1, lam_i),      nu_{m+1} = nu_1 + n
    f_nu    = *_i 1_[L_i(nu)-nu_i, R_i(nu)-nu_i]
    k(a,b)  = sum_{nu in Box(mu), S_nu=|mu|+b} f_nu(a)
"""
from itertools import product
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cyl import is_shape, hstrip, size


def lam_prev(lam, i, n, m):
    """lam_{i-1} with 0-based i, using lam_0 = lam_m - n."""
    return lam[i - 1] if i >= 1 else lam[m - 1] - n


def nu_next(nu, i, n, m):
    """nu_{i+1} with 0-based i, using nu_{m+1} = nu_1 + n."""
    return nu[i + 1] if i + 1 < m else nu[0] + n


def box(mu, n, m):
    rngs = []
    for i in range(m):
        hi = nu_next(mu, i, n, m) - 1
        rngs.append(range(mu[i], hi + 1))
    return product(*rngs)


def LR(nu, lam, n, m):
    """list of (L_i, R_i), 0-based."""
    return [(max(nu[i], lam_prev(lam, i, n, m) + 1),
             min(nu_next(nu, i, n, m) - 1, lam[i]))
            for i in range(m)]


def conv_intervals(widths):
    """convolution of 1_[0,w-1] for w in widths, as a list starting at index 0.
    Returns [] if any width <= 0."""
    if any(w <= 0 for w in widths):
        return []
    cur = [1]
    for w in widths:
        new = [0] * (len(cur) + w - 1)
        # prefix-sum sliding window
        for j, c in enumerate(cur):
            if c:
                for k in range(w):
                    new[j + k] += c
        cur = new
    return cur


def f_nu(nu, lam, n, m):
    """(offset, coeffs) with f_nu(a)=coeffs[a-offset]; (0,[]) if empty."""
    lr = LR(nu, lam, n, m)
    lo = [L - nu[i] for i, (L, R) in enumerate(lr)]
    widths = [R - L + 1 for (L, R) in lr]
    if any(w <= 0 for w in widths):
        return (0, [])
    return (sum(lo), conv_intervals(widths))


def slice_sum(mu, lam, n, m, b):
    """G = a -> k(a,b) as (offset, coeffs) via the box-slice formula."""
    target = sum(mu) + b
    acc = {}
    for nu in box(mu, n, m):
        if sum(nu) != target:
            continue
        off, co = f_nu(nu, lam, n, m)
        for j, c in enumerate(co):
            acc[off + j] = acc.get(off + j, 0) + c
    if not acc:
        return (0, [])
    lo, hi = min(acc), max(acc)
    return (lo, [acc.get(i, 0) for i in range(lo, hi + 1)])


def slice_nus(mu, lam, n, m, b):
    target = sum(mu) + b
    return [nu for nu in box(mu, n, m) if sum(nu) == target]


# ---------- independent instrument: direct chain enumeration ----------

def k_table_chains(mu, lam, n, m):
    """{(a,b): #{mu<nu<kappa<lam hstrip chain, |nu/mu|=b, |kappa/nu|=a}}.
    Uses only is_shape/hstrip/size from cyl.py -- no L_i,R_i."""
    rngs = [range(mu[i], lam[i] + 1) for i in range(m)]
    S = [x for x in product(*rngs) if is_shape(x, n, m)]
    out = {}
    for nu in S:
        if not hstrip(mu, nu, n, m):
            continue
        b = size(mu, nu)
        for kap in S:
            if hstrip(nu, kap, n, m) and hstrip(kap, lam, n, m):
                a = size(nu, kap)
                out[(a, b)] = out.get((a, b), 0) + 1
    return out


# ---------- shape enumeration ----------

def shapes(n, m, mu1=0):
    """cylindric shapes with x_1 = mu1."""
    out = []
    for rest in product(*[range(mu1 + 1, mu1 + n) for _ in range(m - 1)]):
        x = (mu1,) + rest
        if is_shape(x, n, m):
            out.append(x)
    return out


def pairs(n, m, dmax):
    """(mu, lam) with mu_1 = 0, mu <= lam componentwise, lam cylindric, |lam/mu| <= dmax."""
    out = []
    for mu in shapes(n, m, 0):
        rngs = [range(mu[i], mu[i] + dmax + 1) for i in range(m)]
        for lam in product(*rngs):
            if sum(lam) - sum(mu) > dmax:
                continue
            if is_shape(lam, n, m):
                out.append((mu, lam))
    return out


def is_pf2(co):
    """PF2 on a coefficient list (support assumed trimmed to the list)."""
    if not co:
        return True
    for i in range(len(co)):
        left = co[i - 1] if i - 1 >= 0 else 0
        right = co[i + 1] if i + 1 < len(co) else 0
        if co[i] * co[i] < left * right:
            return False
    # interval support
    nz = [i for i, c in enumerate(co) if c != 0]
    if nz and (nz[-1] - nz[0] + 1 != len(nz)):
        return False
    return True
