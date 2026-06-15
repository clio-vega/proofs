"""
d-core and d-quotient of a partition via the abacus (beta-set) method.

Conventions:
  lam is a tuple (descending, possibly with trailing zeros stripped).
  We use a beta-set of size N (a multiple of d, N >= len(lam)):
      beta_i = lam_i + (N - 1 - i),  i = 0..N-1  (lam padded with zeros),
  which is a set of N distinct nonnegative integers.

  d-core: on each of the d runners (residue r mod d), push all beads down to
  the bottom (smallest available positions on that runner); read back the core.

  d-quotient: on runner r, the bead positions p = (b - r)//d for b ≡ r (mod d)
  form a beta-set (of size = #beads on runner r) for the r-th quotient partition.

  The d-quotient depends on N mod d (a cyclic relabelling of runners); we FIX
  N ≡ 0 (mod d) so the runner labels are canonical (runner r ↔ residue r).

Verified below against hand computations:
  4-core((2,2)) = (2,2), 4-quotient = ((),(),(),())   [no 4-hooks]
  2-core((3,1)) = (), and the 2-quotient recovers |lam| = 2|core's worth|...
"""

def _strip(lam):
    return tuple(x for x in lam if x > 0)

def beta_set(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + (N - 1 - i) for i in range(N)]

def from_beta(beta):
    """Recover a partition from a beta-set (any size)."""
    beta = sorted(beta, reverse=True)
    N = len(beta)
    lam = [beta[i] - (N - 1 - i) for i in range(N)]
    return _strip(tuple(lam))

def core_quotient(lam, d):
    """Return (d_core, d_quotient) where d_quotient is a tuple of d partitions."""
    lam = _strip(lam)
    # N: multiple of d, >= len(lam), and big enough headroom.
    N = len(lam) + d
    while N % d != 0:
        N += 1
    beta = beta_set(lam, N)

    # ---- d-quotient ----
    quotient = []
    for r in range(d):
        positions = sorted(((b - r) // d for b in beta if b % d == r), reverse=True)
        k = len(positions)
        # these positions are a beta-set of size k
        part = tuple(positions[i] - (k - 1 - i) for i in range(k))
        quotient.append(_strip(part))

    # ---- d-core: push beads down each runner ----
    core_beta = []
    for r in range(d):
        cnt = sum(1 for b in beta if b % d == r)
        # lowest cnt positions on runner r: r, r+d, r+2d, ...
        core_beta.extend(r + d * t for t in range(cnt))
    core = from_beta(core_beta)
    return core, tuple(quotient)

def d_core(lam, d):
    return core_quotient(lam, d)[0]

def d_quotient(lam, d):
    return core_quotient(lam, d)[1]

# ---------------- tests ----------------
if __name__ == "__main__":
    def size(p): return sum(p)

    # invariant: |lam| = |core| + d * sum|quotient_r|
    from job1_tie_census import partitions
    bad = 0; tot = 0
    for n in range(0, 25):
        for lam in partitions(n) if n > 0 else [()]:
            for d in (2, 3, 4, 5):
                core, quot = core_quotient(lam, d)
                tot += 1
                lhs = size(lam)
                rhs = size(core) + d * sum(size(q) for q in quot)
                if lhs != rhs:
                    bad += 1
                    if bad <= 10:
                        print(f"FAIL size: lam={lam} d={d} core={core} quot={quot} "
                              f"{lhs} != {rhs}")
    print(f"size invariant |lam|=|core|+d*sum|quot|: {tot-bad}/{tot} OK")

    # hand checks
    assert core_quotient((2,2), 4)[0] == (2,2), core_quotient((2,2),4)
    assert all(q == () for q in core_quotient((2,2),4)[1]), core_quotient((2,2),4)
    # (4): 4-hook present, 4-core should be empty, quotient size 1
    c,q = core_quotient((4,), 4)
    assert c == (), (c,q)
    assert sum(size(x) for x in q) == 1, (c,q)
    # 2-core of (3,1) is empty (staircase 2-core are staircases (k,k-1,...,1))
    c,_ = core_quotient((3,1), 2)
    assert c == (), c
    # 2-core of (3,2,1) is itself (a staircase) ; |lam|=6 even so quotient size 0 -> core=lam
    c,q = core_quotient((3,2,1), 2)
    assert c == (3,2,1) and all(x==() for x in q), (c,q)
    print("hand checks PASSED")
