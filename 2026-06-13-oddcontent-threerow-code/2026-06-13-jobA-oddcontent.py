#!/usr/bin/env python3
"""
Job A — the odd-content probe (decisive).

Armon-Swanson factor f^lambda(q,t) = f^lambda(q) * prod_{cell}(1 + t q^{c(cell)}).
At q=-1, (1+t q^c) vanishes (in t-adic order) iff c(cell) is ODD.  So the order of
the content product at q=-1 is OC(lambda) := #{cells : content odd}.

Proved order law:  ord_{q=-1} G_lambda = floor(|2-core(lambda)|/2|).

TEST: is ord_{q=-1} G_lambda a clean function of OC(lambda) or OC(core_2(lambda))?

Pure-Python (no SageMath in this container).
"""
from functools import lru_cache
from sympy import symbols, Poly, factor, QQ, Rational
import sympy

# ---------------------------------------------------------------------------
# partitions
# ---------------------------------------------------------------------------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest

# ---------------------------------------------------------------------------
# contents and odd-content count
# ---------------------------------------------------------------------------
def odd_content_count(lam):
    """#{cells (i,j): content j-i is odd}, i,j 0-indexed (row i, col j)."""
    oc = 0
    for i, part in enumerate(lam):
        for j in range(part):
            if (j - i) % 2 != 0:
                oc += 1
    return oc

# ---------------------------------------------------------------------------
# 2-core by repeated domino (2-rim-hook) removal
# ---------------------------------------------------------------------------
def two_core(lam):
    lam = list(lam)
    changed = True
    while changed:
        changed = False
        # try to remove a 2-rim-hook (border strip of size 2)
        rh = remove_a_rimhook(lam, 2)
        if rh is not None:
            lam = rh
            changed = True
    return tuple(lam)

def remove_a_rimhook(lam, size):
    """Remove one border strip of given size if possible; return new partition list or None.
    Uses beta-numbers: a rim hook of length h corresponds to b_i -> b_i - h being a 'gap'."""
    lam = [p for p in lam if p > 0]
    k = len(lam)
    if k == 0:
        return None
    # beta numbers b_i = lam_i + (k-1-i), strictly decreasing, i=0..k-1
    beta = [lam[i] + (k - 1 - i) for i in range(k)]
    betaset = set(beta)
    for idx in range(k):
        b = beta[idx]
        if b - size >= 0 and (b - size) not in betaset:
            newbeta = beta[:]
            newbeta[idx] = b - size
            newbeta.sort(reverse=True)
            # back to partition
            newlam = [newbeta[i] - (k - 1 - i) for i in range(k)]
            newlam = [p for p in newlam if p > 0]
            return newlam
    return None

# ---------------------------------------------------------------------------
# parity-twisted SYT descent statistic  s(T) = sum_{i in Des(T)} w_i,
#   w_i = 2i-1 if n-i odd else 0   (positions i=1..n-1)   [from memory]
# G_lambda(q) = sum_{T in SYT(lambda)} q^{s(T)}.   Cross-check: ord at q=-1.
# ---------------------------------------------------------------------------
def syt_list(lam):
    """Generate all SYT of shape lam as filling dicts; yield list of (row of value v)."""
    cells = [(i, j) for i, p in enumerate(lam) for j in range(p)]
    n = len(cells)
    # build by insertion of 1..n into addable cells
    results = []
    shape = [0] * (len(lam) + 1)  # current row lengths
    target = list(lam)
    rowofval = [0] * (n + 1)  # rowofval[v] = row index of value v (1..n)

    def rec(v):
        if v > n:
            results.append(tuple(rowofval[1:]))
            return
        for r in range(len(target)):
            # can place in row r at column shape[r] if shape[r] < target[r] and (r==0 or shape[r] < shape[r-1])
            if shape[r] < target[r] and (r == 0 or shape[r] < shape[r - 1]):
                shape[r] += 1
                rowofval[v] = r
                rec(v + 1)
                shape[r] -= 1
    rec(1)
    return results, n

def s_statistic(rowofval, n):
    # descent at i if value i+1 is in a strictly lower row than value i
    s = 0
    for i in range(1, n):
        if rowofval[i] > rowofval[i - 1]:  # rowofval indexed 0..n-1 for values 1..n
            # i is a descent (value i in lower row? value i+1 below value i)
            w = (2 * i - 1) if ((n - i) % 2 == 1) else 0
            s += w
    return s

def G_poly(lam):
    q = symbols('q')
    sl, n = syt_list(lam)
    poly = 0
    for rv in sl:
        poly += q ** s_statistic(rv, n)
    return Poly(poly, q), q

def ord_at_minus1_from_G(lam):
    poly, q = G_poly(lam)
    # multiplicity of root q=-1
    if poly.is_zero:
        return None
    mult = 0
    p = poly
    while p.eval(-1) == 0:
        p = p.diff(q)
        mult += 1
        if mult > 200:
            break
    # better: factor multiplicity via division
    return _root_mult(poly, q, -1)

def _root_mult(poly, q, root):
    mult = 0
    p = poly
    while True:
        if p.eval(root) == 0 and not p.is_zero:
            quo, rem = sympy.div(p.as_expr(), (q - root), q)
            p = Poly(quo, q)
            mult += 1
        else:
            break
        if mult > 500:
            break
    return mult

# ---------------------------------------------------------------------------
# Main probe
# ---------------------------------------------------------------------------
def floor_div(a, b):
    return a // b

def staircase_OC_closed(k):
    """OC(delta_k) for delta_k=(k,k-1,...,1)."""
    return odd_content_count(tuple(range(k, 0, -1)))

def run():
    print("=" * 78)
    print("STAIRCASE (2-core) ODD-CONTENT vs ORDER LAW  ord=floor(|delta_k|/2)")
    print("=" * 78)
    print(f"{'k':>3} {'|delta_k|':>9} {'floor(|.|/2)':>12} {'OC(delta_k)':>12} {'match?':>7}")
    sc_rows = []
    for k in range(0, 9):
        dk = tuple(range(k, 0, -1)) if k > 0 else ()
        size = k * (k + 1) // 2
        ordv = size // 2
        oc = odd_content_count(dk)
        sc_rows.append((k, size, ordv, oc, ordv == oc))
        print(f"{k:>3} {size:>9} {ordv:>12} {oc:>12} {str(ordv==oc):>7}")
    print()

    # Full census n<=14 using proved closed form; cross-check ord via SYT for n<=9
    print("=" * 78)
    print("FULL CENSUS: ord_{q=-1}G = floor(|2-core|/2)  vs  OC(lambda), OC(core_2)")
    print("=" * 78)
    mism_lam = 0   # ord vs OC(lambda)
    mism_core = 0  # ord vs OC(core_2)
    total = 0
    crosscheck_fail = 0
    examples = []
    for n in range(1, 15):
        for lam in partitions(n):
            core = two_core(lam)
            csize = sum(core)
            ord_closed = csize // 2
            oc_lam = odd_content_count(lam)
            oc_core = odd_content_count(core)
            total += 1
            if ord_closed != oc_lam:
                mism_lam += 1
            if ord_closed != oc_core:
                mism_core += 1
                if len(examples) < 12:
                    examples.append((lam, core, ord_closed, oc_core))
            # cross-check the closed form against actual G for small n
            if n <= 9:
                ord_syt = ord_at_minus1_from_G(lam)
                if ord_syt != ord_closed:
                    crosscheck_fail += 1
                    if crosscheck_fail <= 8:
                        print(f"  CROSSCHECK FAIL lam={lam}: SYT ord={ord_syt} closed={ord_closed}")
    print(f"\nTotal partitions n<=14: {total}")
    print(f"Cross-check (SYT G vs closed form, n<=9) failures: {crosscheck_fail}")
    print(f"ord  !=  OC(lambda)   mismatches: {mism_lam}")
    print(f"ord  !=  OC(core_2)   mismatches: {mism_core}")
    if examples:
        print("\nExamples ord != OC(core_2):")
        for lam, core, o, oc in examples:
            print(f"  lam={lam} core={core} ord={o} OC(core)={oc}")
    return sc_rows

if __name__ == "__main__":
    run()
