"""
Job B (2026-06-15 evening) — PROBE-FIRST: spin/cospin candidate for the 2-core
order law  ord_{q=-1} G_lambda = floor(|2-core(lambda)|/2).

Discipline (CODE.md): build the candidate, run THE (3,2,1) probe, decide. Do NOT
scaffold a dead route.

Candidate.  Spin generating function of standard domino tableaux of lambda
(Carre-Leclerc / Stanton-White / Schilling-Shimozono-White): a standard domino
tableau is a chain  core = mu_0 subset mu_1 subset ... subset mu_k = lambda  where
each mu_i / mu_{i-1} is a single domino (2 adjacent cells); spin(D) = #vertical
dominoes; S_lambda(q) = sum_D q^{spin(D)}.  Natural readings of S_lambda as a
candidate(lambda):
    deg  = max spin (top degree)
    ord  = order of vanishing of S_lambda(q) at q = -1  (the CSP/sieve reading)
    val  = S_lambda(-1)
We compare each to the KNOWN order law floor(|2-core|/2).
"""
from functools import lru_cache
from core_quotient import core_quotient

def _strip(p):
    return tuple(x for x in p if x > 0)

def add_dominoes(mu):
    """All partitions obtained from mu by adding ONE domino (H or V), with spin
    contribution (0 horizontal, 1 vertical)."""
    mu = list(mu)
    r = len(mu)
    out = []
    # work with rows 0..r (allow a new row)
    rows = mu + [0, 0]
    L = len(rows)
    # horizontal domino in row i: add 2 boxes at end of row i
    for i in range(L):
        new = rows[:]
        new[i] += 2
        if _is_partition(new):
            out.append((_strip(tuple(new)), 0))
    # vertical domino spanning rows i, i+1: add 1 box to each (same column)
    for i in range(L - 1):
        new = rows[:]
        # the new cell in row i is at column rows[i]+1; in row i+1 at rows[i+1]+1;
        # they must be the same column => rows[i] == rows[i+1]
        if rows[i] == rows[i + 1]:
            new[i] += 1
            new[i + 1] += 1
            if _is_partition(new):
                out.append((_strip(tuple(new)), 1))
    return out

def _is_partition(rows):
    prev = float('inf')
    for x in rows:
        if x < 0 or x > prev:
            return False
        prev = x
    return True

@lru_cache(maxsize=None)
def spin_gen(lam):
    """Return dict spin -> count of standard domino tableaux of shape lam
    (built up from its 2-core)."""
    lam = _strip(lam)
    core, _ = core_quotient(lam, 2)
    return _spin_from(core, lam)

@lru_cache(maxsize=None)
def _spin_from(start, target):
    start = _strip(start); target = _strip(target)
    if start == target:
        return {0: 1}
    if sum(start) > sum(target):
        return {}
    acc = {}
    for nxt, sp in add_dominoes(start):
        if sum(nxt) > sum(target):
            continue
        if not _contained(nxt, target):
            continue
        sub = _spin_from(nxt, target)
        for s, c in sub.items():
            acc[s + sp] = acc.get(s + sp, 0) + c
    return acc

def _contained(mu, lam):
    mu = list(mu); lam = list(lam)
    if len(mu) > len(lam):
        return False
    return all(mu[i] <= lam[i] for i in range(len(mu)))

# ---- candidate readings ----
def poly_from_dict(d):
    if not d:
        return [0]
    deg = max(d)
    return [d.get(i, 0) for i in range(deg + 1)]

def eval_at(coeffs, x):
    return sum(c * x**i for i, c in enumerate(coeffs))

def ord_at_minus1(coeffs):
    """multiplicity of (q+1) as a factor, i.e. order of vanishing at q=-1."""
    if all(c == 0 for c in coeffs):
        return float('inf')
    c = coeffs[:]
    k = 0
    # synthetic division by (q+1): repeatedly while value at -1 is 0
    while eval_at(c, -1) == 0 and len(c) > 1:
        # divide by (q+1)
        new = []
        rem = 0
        for coeff in reversed(c):
            rem = coeff - rem  # since root -1: synthetic division by (q-(-1))
            new.append(rem)
        new = list(reversed(new))
        # new[0] is remainder (should be 0)
        c = new[1:]
        k += 1
    return k

def candidate(lam):
    d = spin_gen(lam)
    coeffs = poly_from_dict(d)
    return dict(spin_gen=d, deg=max(d) if d else 0,
                val_at_m1=eval_at(coeffs, -1),
                ord_at_m1=ord_at_minus1(coeffs),
                ndom=(sum(_strip(lam)) - sum(core_quotient(lam, 2)[0])) // 2)

def orderlaw(lam):
    core, _ = core_quotient(lam, 2)
    return sum(core) // 2  # floor(|2-core|/2)

# =============================================================== THE PROBE
if __name__ == "__main__":
    print("=" * 72)
    print("THE (3,2,1) PROBE  — candidate(spin) vs order law floor(|2-core|/2)")
    print("=" * 72)
    witnesses = [
        ("delta_3 = (3,2,1)", (3, 2, 1)),
        ("delta_4 = (4,3,2,1)", (4, 3, 2, 1)),
        ("delta_5 = (5,4,3,2,1)", (5, 4, 3, 2, 1)),
        ("delta_6", (6, 5, 4, 3, 2, 1)),
        ("(2,2)  [bounded]", (2, 2)),
        ("(3,1)  [bounded]", (3, 1)),
        ("(4,4)  [bounded]", (4, 4)),
        ("(2,1)  [a 2-core]", (2, 1)),
    ]
    print(f"{'shape':22} {'orderlaw':>8} | {'cand.deg':>8} {'cand.ord@-1':>11} "
          f"{'cand.val@-1':>11} {'#dominoes':>9}  spin_gen")
    for name, lam in witnesses:
        ol = orderlaw(lam)
        c = candidate(lam)
        print(f"{name:22} {ol:>8} | {c['deg']:>8} {c['ord_at_m1']:>11} "
              f"{c['val_at_m1']:>11} {c['ndom']:>9}  {dict(sorted(c['spin_gen'].items()))}")

    print()
    print("=" * 72)
    print("DISCRIMINATION:  does the candidate READ THE 2-CORE?")
    print("=" * 72)
    print("(a) SAME 2-core, DIFFERENT 2-quotient  -> order law MUST agree, spin may differ")
    # build pairs with same 2-core
    from itertools import groupby
    pairs_same_core = [
        ((4, 2), (2, 2, 1, 1)),     # both 2-core () ?  check
        ((4,), (2, 2)),
        ((5, 1), (3, 3)),
        ((6, 2), (4, 4)),
    ]
    for l1, l2 in pairs_same_core:
        c1, _ = core_quotient(l1, 2); c2, _ = core_quotient(l2, 2)
        a, b = candidate(l1), candidate(l2)
        flag = "SAME-CORE" if c1 == c2 else "diff-core"
        print(f"  {str(l1):10} core={str(c1):8} ol={orderlaw(l1)} spin_deg={a['deg']} val@-1={a['val_at_m1']} | "
              f"{str(l2):12} core={str(c2):8} ol={orderlaw(l2)} spin_deg={b['deg']} val@-1={b['val_at_m1']}  [{flag}]")

    print()
    print("(b) DIFFERENT 2-core, SAME 2-quotient -> order law MUST differ; spin (on quotient) identical")
    # (3,2,1) has core (3,2,1), quotient empty.  Add the SAME quotient to different cores:
    # core () with no dominoes = (), core (1) ... Actually compare 2-cores of different size,
    # all with EMPTY quotient (no dominoes): every 2-core has empty quotient & spin_gen={0:1}.
    for lam in [(1,), (2, 1), (3, 2, 1), (4, 3, 2, 1), (5, 4, 3, 2, 1)]:
        core, quot = core_quotient(lam, 2)
        c = candidate(lam)
        print(f"  {str(lam):14} core={str(core):14} |core|/2={orderlaw(lam):>2}  "
              f"quotient={quot}  spin_gen={dict(sorted(c['spin_gen'].items()))}  cand.deg={c['deg']}")
