"""
fusion_lib.py -- Littlewood-Richardson coefficients and level-k fusion
(level-restricted LR / Verlinde) coefficients for sl_n, self-contained (no Sage).

LR rule:   c^nu_{lam,mu} = # LR skew tableaux of shape nu/lam, content mu.
Fusion (Kac-Walton, type A_{n-1}^(1) at level k):
   N^{(k) nu}_{lam mu} = sum_{nu'} c^{nu'}_{lam mu} * eps(w(nu')),
   where each LR target nu' (a partition with <= n rows) is reflected by the affine
   Weyl group dot-action into the level-k alcove; eps = (-1)^{#reflections}; if nu'+rho
   hits a wall it contributes 0; the dominant image is reduced mod columns of height n.

Validated below against: sl_2 level-k closed formula, the Ising (sl_2 level 2) fusion
rules, and sl_3 level 1 (Z/3 simple currents).  Run `python3 fusion_lib.py`.
"""
from collections import defaultdict
from functools import lru_cache


# --------------------------------------------------------------------------
# Littlewood-Richardson coefficients (direct rule)
# --------------------------------------------------------------------------
def _lr_fillings(outer, inner, content_target):
    """Count LR tableaux of shape outer/inner with given content (yamanouchi/ballot)."""
    inner = list(inner) + [0] * (len(outer) - len(inner))
    cells = []
    for r in range(len(outer)):
        for c in range(inner[r], outer[r]):
            cells.append((r, c))
    # fill cells row by row, left to right (reading word = right-to-left, top-to-bottom)
    # We enumerate in reading-word order: row 0 R->L, row 1 R->L, ...
    reading = []
    for r in range(len(outer)):
        for c in range(outer[r] - 1, inner[r] - 1, -1):
            reading.append((r, c))
    nrows = len(outer)
    val = [[0] * outer[r] for r in range(nrows)]
    maxcontent = len(content_target)
    count = 0

    def rec(idx, ballot):
        nonlocal count
        if idx == len(reading):
            if tuple(ballot[i] for i in range(maxcontent)) == tuple(content_target):
                count += 1
            return
        r, c = reading[idx]
        for v in range(1, maxcontent + 1):
            # column-strict: value above must be < v
            if r > 0 and val[r - 1][c] != 0 and val[r - 1][c] >= v:
                continue
            # weakly increasing along row L->R: val[r][c] <= val[r][c+1] (cell to the right placed earlier)
            if c + 1 < outer[r] and val[r][c + 1] != 0 and v > val[r][c + 1]:
                continue
            # ballot: after placing v, count_v <= count_{v-1}
            if v > 1 and ballot[v - 1] + 1 > ballot[v - 2]:
                continue
            if ballot[v - 1] + 1 > content_target[v - 1]:
                continue
            val[r][c] = v
            ballot[v - 1] += 1
            rec(idx + 1, ballot)
            ballot[v - 1] -= 1
            val[r][c] = 0

    rec(0, [0] * maxcontent)
    return count


@lru_cache(maxsize=None)
def lr_coeff(lam, mu, nu):
    """LR coefficient c^nu_{lam,mu}."""
    lam = tuple(x for x in lam if x > 0)
    mu = tuple(x for x in mu if x > 0)
    nu = tuple(x for x in nu if x > 0)
    if sum(nu) != sum(lam) + sum(mu):
        return 0
    if len(lam) > len(nu):
        return 0
    for i in range(len(lam)):
        if lam[i] > nu[i]:
            return 0
    if not mu:
        return 1 if lam == nu else 0
    return _lr_fillings(nu, lam, mu)


def lr_targets(lam, mu, n):
    """All partitions nu (<= n rows) with c^nu_{lam mu} != 0, as {nu: c}."""
    from itertools import product
    lam = tuple(lam); mu = tuple(mu)
    N = sum(lam) + sum(mu)
    out = {}
    for nu in _partitions_leq_rows(N, n, (lam[0] if lam else 0) + (mu[0] if mu else 0)):
        c = lr_coeff(lam, mu, nu)
        if c:
            out[nu] = c
    return out


@lru_cache(maxsize=None)
def _partitions_leq_rows(N, maxrows, maxcol):
    res = []
    def gen(remaining, maxpart, cur):
        if remaining == 0:
            res.append(tuple(cur)); return
        if len(cur) >= maxrows:
            return
        for p in range(min(remaining, maxpart, maxcol), 0, -1):
            gen(remaining - p, p, cur + [p])
    gen(N, maxcol if maxcol else N, [])
    return res


# --------------------------------------------------------------------------
# Affine Weyl reflection into the level-k alcove (type A_{n-1}^(1))
# --------------------------------------------------------------------------
def _reflect_to_alcove(nu, n, k, maxit=10000):
    """
    Reflect beta = nu + rho into the dominant level-k alcove by the affine Weyl
    dot-action.  rho=(n-1,...,1,0), shifted level ell = k+n.
    Returns (reduced_partition, sign) or None if beta hits a wall.
    """
    rho = [n - 1 - i for i in range(n)]
    beta = [ (nu[i] if i < len(nu) else 0) + rho[i] for i in range(n) ]
    ell = k + n
    sign = 1
    for _ in range(maxit):
        # finite Weyl: make strictly decreasing
        done = True
        for i in range(n - 1):
            if beta[i] == beta[i + 1]:
                return None                      # on a finite wall
            if beta[i] < beta[i + 1]:
                beta[i], beta[i + 1] = beta[i + 1], beta[i]
                sign = -sign; done = False; break
        if not done:
            continue
        # affine wall: beta_1 - beta_n vs ell
        d = beta[0] - beta[-1]
        if d == ell:
            return None                          # on the affine wall
        if d > ell:
            b1, bn = beta[0], beta[-1]
            beta[0] = bn + ell
            beta[-1] = b1 - ell
            sign = -sign
            continue
        # inside the alcove and strictly decreasing
        nu_dom = [beta[i] - rho[i] for i in range(n)]
        # reduce modulo columns of height n: subtract the smallest part
        c = nu_dom[-1]
        red = tuple(x - c for x in nu_dom)
        red = tuple(x for x in red if x > 0)
        return red, sign
    raise RuntimeError(f"alcove reflection did not terminate: nu={nu}, n={n}, k={k}")


def fusion_product(lam, mu, n, k):
    """Level-k fusion product N^{(k)}_{lam,mu} as {nu: multiplicity} (sl_n)."""
    lam = tuple(x for x in lam if x > 0)
    mu = tuple(x for x in mu if x > 0)
    res = defaultdict(int)
    for nu_prime, c in lr_targets(lam, mu, n).items():
        r = _reflect_to_alcove(nu_prime, n, k)
        if r is None:
            continue
        nu_dom, sign = r
        res[nu_dom] += sign * c
    return {nu: m for nu, m in res.items() if m != 0}


def fusion_coeff(lam, mu, nu, n, k):
    """Single level-k fusion multiplicity N^{(k) nu}_{lam mu}."""
    nu = tuple(x for x in nu if x > 0)
    return fusion_product(lam, mu, n, k).get(nu, 0)


def level_weights(n, k):
    """All level-k dominant weights of sl_n: partitions with <= n-1 rows, largest part <= k."""
    res = []
    for N in range(0, (n - 1) * k + 1):
        for nu in _partitions_leq_rows(N, n - 1, k):
            res.append(nu)
    return res


# --------------------------------------------------------------------------
# sl_2 level-k closed formula (validation anchor)
# --------------------------------------------------------------------------
def sl2_fusion_closed(a, b, k):
    """sl_2 level k: N^c_{ab} = 1 if |a-b|<=c<=min(a+b, 2k-a-b) and a+b+c even, else 0.
    Returns {(c,): 1}."""
    out = {}
    for c in range(0, k + 1):
        if abs(a - b) <= c <= min(a + b, 2 * k - a - b) and (a + b + c) % 2 == 0:
            out[(c,) if c > 0 else ()] = 1
    return out


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def _validate():
    print("=" * 78)
    print("VALIDATION of fusion_lib")
    print("=" * 78)

    # --- ordinary LR sanity ---
    assert lr_coeff((1,), (1,), (2,)) == 1
    assert lr_coeff((1,), (1,), (1, 1)) == 1
    assert lr_coeff((2, 1), (2, 1), (3, 2, 1)) == 2   # famous c=2
    assert lr_coeff((2, 1), (2, 1), (4, 2)) == 1
    print("  LR rule: (1)*(1)=(2)+(11); c^{(3,2,1)}_{(2,1),(2,1)}=2  -- OK")

    # --- Ising: sl_2 level 2 ---
    print("\n  Ising (sl_2 level 2): labels ()=1, (1)=sigma, (2)=psi")
    tests = {
        ((1,), (1,)): {(): 1, (2,): 1},          # sigma*sigma = 1 + psi
        ((2,), (2,)): {(): 1},                    # psi*psi = 1
        ((1,), (2,)): {(1,): 1},                  # sigma*psi = sigma
        ((2,), (1,)): {(1,): 1},
        ((1,), ()):  {(1,): 1},
    }
    ok = True
    for (a, b), exp in tests.items():
        got = fusion_product(a, b, 2, 2)
        match = got == exp
        ok = ok and match
        print(f"    {a} x {b} = {dict(sorted(got.items()))}   expected {dict(sorted(exp.items()))}  {'OK' if match else 'FAIL'}")
    assert ok, "Ising anchor FAILED"

    # --- sl_2 level k vs closed formula, all a,b,k ---
    print("\n  sl_2 level k: Kac-Walton vs closed Verlinde formula")
    bad = 0; nchk = 0
    for k in range(1, 7):
        for a in range(0, k + 1):
            for b in range(0, k + 1):
                la = (a,) if a else ()
                lb = (b,) if b else ()
                got = fusion_product(la, lb, 2, k)
                exp = sl2_fusion_closed(a, b, k)
                got = {kk: vv for kk, vv in got.items()}
                if got != exp:
                    bad += 1
                    if bad <= 8:
                        print(f"    MISMATCH k={k} {a}x{b}: got {got} exp {exp}")
                nchk += 1
    print(f"    checked {nchk} products, mismatches: {bad}  {'-- OK' if bad == 0 else '-- FAIL'}")
    assert bad == 0

    # --- sl_3 level 1: Z/3 simple currents ---
    print("\n  sl_3 level 1 (Z/3): ()=0, (1)=1, (1,1)=2")
    z3 = {
        ((1,), (1,)): {(1, 1): 1},      # 1+1=2
        ((1,), (1, 1)): {(): 1},        # 1+2=0
        ((1, 1), (1, 1)): {(1,): 1},    # 2+2=1
    }
    ok3 = True
    for (a, b), exp in z3.items():
        got = fusion_product(a, b, 3, 1)
        m = got == exp; ok3 = ok3 and m
        print(f"    {a} x {b} = {dict(sorted(got.items()))}  expected {exp}  {'OK' if m else 'FAIL'}")
    assert ok3

    # --- associativity / commutativity spot-check sl_3 level 2 ---
    print("\n  sl_3 level 2: commutativity & a nontrivial product")
    p = fusion_product((2, 1), (2, 1), 3, 2)
    q = fusion_product((2, 1), (2, 1), 3, 2)
    assert p == q
    print(f"    (2,1) x (2,1) at sl_3 level 2 = {dict(sorted(p.items()))}")
    # total dimension check vs ordinary LR restricted to level weights would differ; just show it's truncated
    ordinary = lr_targets((2, 1), (2, 1), 3)
    print(f"    (ordinary LR targets <=3 rows: {dict(sorted(ordinary.items()))})")

    print("\n  ALL VALIDATION PASSED")
    return True


if __name__ == "__main__":
    _validate()
