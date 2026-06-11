"""
t-core and t-quotient of a partition via the abacus (beta-set) construction.

Conventions
-----------
Take lambda = (l_1 >= l_2 >= ... >= l_ell), pad to length N (a multiple of t,
with N >= ell).  Beta-set  beta_i = l_i + (N - i),  i = 1..N  (distinct, >= 0).
Place a bead at each beta_i on an abacus with t runners (runner r holds the
beads with position == r mod t).

  * t-quotient: on runner r, the bead positions are  r, r+t, r+2t, ...; their
    "heights" q = (pos - r)//t form a beta-set whose partition is lambda^{(r)}.
  * t-core: slide every bead up its runner as far as it goes; read off the
    resulting partition.

The labelling of the runners (hence the cyclic order of the quotient) depends on
N mod t; we fix N divisible by t so the runners are r = 0..t-1 in the natural
order.  This matches the convention used in the rectangular-class character
formula chi^lambda(t^w) = eps * multinomial * prod f^{lambda^{(r)}}.

Verified against sympy / direct Murnaghan-Nakayama in __main__.
"""
from functools import lru_cache


def _beta_set(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + (N - 1 - i) for i in range(N)]   # i from 0: beta = l_i + (N-1-i)


def _partition_from_beta(beta):
    beta = sorted(beta, reverse=True)
    n = len(beta)
    lam = [beta[i] - (n - 1 - i) for i in range(n)]
    return tuple(x for x in lam if x > 0)


def t_quotient(lam, t, N=None):
    """Return (core, [lambda^{(0)},...,lambda^{(t-1)}]) for the given t."""
    lam = tuple(x for x in lam if x > 0)
    ell = len(lam)
    if N is None:
        N = ((ell // t) + 2) * t          # multiple of t, comfortably >= ell
    assert N % t == 0 and N >= ell
    beta = _beta_set(lam, N)
    # quotient
    quot = []
    for r in range(t):
        runner = sorted((b - r) // t for b in beta if b % t == r)
        # runner positions are a beta-set (heights) -> partition
        quot.append(_partition_from_beta(runner))
    # core: slide beads up each runner
    core_beta = []
    for r in range(t):
        heights = sorted((b - r) // t for b in beta if b % t == r)
        k = len(heights)
        # slid-up heights are 0,1,...,k-1 -> positions r, r+t, ..., r+(k-1)t
        core_beta += [r + h * t for h in range(k)]
    core = _partition_from_beta(core_beta)
    return core, quot


def f_syt(lam):
    """Number of standard Young tableaux of shape lam (hook length formula)."""
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    if n == 0:
        return 1
    import math
    prod = 1
    cols = _conj(lam)
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - 1 - j
            leg = cols[j] - 1 - i
            prod *= (arm + leg + 1)
    return math.factorial(n) // prod


def _conj(lam):
    if not lam:
        return ()
    m = max(lam)
    return tuple(sum(1 for x in lam if x > j) for j in range(m))


if __name__ == "__main__":
    from job1_tie_census import mn_character, partitions
    import math

    # --- sanity: 2-quotient of small shapes ---
    for lam in [(2, 2), (3, 1), (2, 1, 1), (4,), (1, 1, 1, 1), (3, 3), (4, 2)]:
        core, quot = t_quotient(lam, 2)
        print(f"lam={lam}: 2-core={core} 2-quot={quot} "
              f"(|q0|+|q1|={sum(map(sum,quot))})")

    # --- verify rectangular character formula chi^lambda(2^m) = eps*C(m;|q0|,|q1|)*f^q0*f^q1 ---
    print("\nRectangular-class formula chi^lambda(2^m) vs 2-quotient product:")
    ok = bad = empty_core_with_zero = 0
    for m in range(1, 9):
        for lam in partitions(2 * m):
            chi2m = mn_character(lam, tuple([2] * m))
            core, quot = t_quotient(lam, 2)
            if core != ():
                # nonempty 2-core => chi should vanish
                if chi2m != 0:
                    bad += 1
                continue
            a, b = sum(quot[0]), sum(quot[1])
            mult = math.comb(m, a)  # = C(m; a,b) since a+b=m
            prod = mult * f_syt(quot[0]) * f_syt(quot[1])
            if abs(chi2m) == prod:
                ok += 1
            else:
                bad += 1
                if bad <= 10:
                    print(f"  MISMATCH lam={lam}: chi={chi2m} pred=+-{prod} "
                          f"quot={quot}")
    print(f"  |chi^lambda(2^m)| == C(m;|q0|,|q1|) f^q0 f^q1 : ok={ok} bad={bad}")
