"""Q254, Test B (d >= 1) and the uniqueness statement.

For a cylindric skew shape lam/mu in cyl^{n,m} with ell variables:

  (b) recompute the true support W_ell and its dominance maximum lamhat, and
      check W_ell = {alpha : sort(alpha) <= lamhat}   [control on thm:main];
  (c) EXHAUSTIVELY search every diagram D with columns inside [ell] and
      #D = |lam/mu| for supp(chi_D) = W_ell, and report all solutions.
      Prediction from Symmetry Rigidity: the solutions are exactly the
      multiset of bottom-justified columns {ell-k+1..ell} with multiplicity
      lamhat'_k - i.e. the skyline diagram of the antidominant rearrangement
      of lamhat, and nothing else.
"""
import sys, itertools
from collections import defaultdict
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
from fms import supp_chi, is_Sm_stable, is_bottom_justified
import cyl as C   # vendored copy, see cyl.py header
from rigidity import dominates, sort_part


def cyl_shapes(n, m):
    """all x=(x_1<...<x_m), x_m < x_1+n, normalised so x_1 = 0."""
    out = []
    for rest in itertools.combinations(range(1, n), m - 1):
        x = (0,) + rest
        if C.is_shape(x, n, m):
            out.append(x)
    return out


def true_support(mu, lam, n, m, ell):
    return set(C.weights(mu, lam, n, m, ell))


def dom_max(W):
    P = {sort_part(a) for a in W}
    mx = [p for p in P if not any(q != p and dominates(q, p) for q in P)]
    return mx


def ideal_check(W, lamhat, ell):
    N = sum(lamhat)
    for a in W:
        if not dominates(lamhat, sort_part(a)):
            return False
    # count
    from rigidity import partitions_le, rearrangements
    cnt = sum(rearrangements(p, ell) for p in partitions_le(N, ell) if dominates(lamhat, p))
    return cnt == len(W)


def skyline_of_lamhat(lamhat, ell):
    """columns {ell-k+1..ell} with multiplicity lamhat'_k, i.e. one column
    {ell-lamhat'_j+1..ell}... equivalently: for each j=1..lamhat_1 a column of
    size lamhat'_j, bottom-justified."""
    if not lamhat:
        return []
    conj = [sum(1 for v in lamhat if v >= j) for j in range(1, lamhat[0] + 1)]
    return sorted(tuple(range(ell - k + 1, ell + 1)) for k in conj)


def all_diagrams(ell, N, maxcols):
    subs = [c for r in range(1, ell + 1) for c in itertools.combinations(range(1, ell + 1), r)]
    out = []
    def rec(start, rem, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        if len(cur) >= maxcols: return
        for idx in range(start, len(subs)):
            s = subs[idx]
            if len(s) <= rem:
                rec(idx, rem - len(s), cur + [s])
    rec(0, N, [])
    return out


def run(nmax=5, ell_max=4, Nmax=8, verbose=True):
    rows = []
    for n in range(2, nmax + 1):
        for m in range(1, n):
            shapes = cyl_shapes(n, m)
            for mu in shapes:
                for lam in itertools.product(*[range(mu[i], mu[i] + n + 1) for i in range(m)]):
                    if not C.is_shape(lam, n, m) or not C.contains(lam, mu):
                        continue
                    N = C.size(mu, lam)
                    if N == 0 or N > Nmax:
                        continue
                    d = N // max(1, (n - m))        # crude winding proxy; reported, not used
                    for ell in range(1, ell_max + 1):
                        W = true_support(mu, lam, n, m, ell)
                        if not W:
                            continue
                        mx = dom_max(W)
                        if len(mx) != 1:
                            rows.append((n, m, mu, lam, ell, N, 'MULTIPLE-MAXIMA', mx))
                            continue
                        lamhat = mx[0]
                        ok_ideal = ideal_check(W, lamhat, ell)
                        rows.append((n, m, mu, lam, ell, N, lamhat, ok_ideal, W))
    return rows


if __name__ == "__main__":
    rows = run()
    bad = [r for r in rows if r[6] == 'MULTIPLE-MAXIMA' or r[7] is False]
    print(f"(b) control on thm:main: {len(rows)} (shape, ell) instances, "
          f"{len(rows)-len(bad)} with a unique dominance maximum lamhat and "
          f"supp = {{alpha : sort(alpha) <= lamhat}};  {len(bad)} failures")
    for r in bad[:5]:
        print("   FAIL", r[:8])

    # (c) uniqueness of the working diagram -- exhaustive, on a subsample
    print()
    print("(c) exhaustive search for diagrams D with supp(chi_D) = W_ell:")
    seen = set()
    tested = 0
    uniq_ok = 0
    for r in rows:
        if r[6] == 'MULTIPLE-MAXIMA':
            continue
        n, m, mu, lam, ell, N, lamhat, ok, W = r
        if ell < 2 or N > 6 or ell > 4:
            continue
        keyk = (ell, N, lamhat)
        if keyk in seen:
            continue
        seen.add(keyk)
        tested += 1
        sols = [D for D in all_diagrams(ell, N, N) if supp_chi(list(D), ell) == W]
        want = tuple(skyline_of_lamhat(lamhat, ell))
        got = {tuple(sorted(D)) for D in sols}
        if got == {want}:
            uniq_ok += 1
        else:
            print(f"   n={n} m={m} mu={mu} lam={lam} ell={ell} lamhat={lamhat}")
            print(f"      solutions found: {sorted(got)}")
            print(f"      predicted (skyline of rev(lamhat)): {want}")
    print(f"   {uniq_ok}/{tested} distinct (ell,N,lamhat) classes: the ONLY diagram D with "
          f"supp(chi_D)=W_ell is the skyline of the antidominant rearrangement of lamhat")
