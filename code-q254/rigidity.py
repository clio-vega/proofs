"""Q254, Test 1: SYMMETRY RIGIDITY, exhaustive.

CLAIM to be falsified:  Let D be a nonempty diagram and m = max(union of its
columns).  Then supp(chi_D) is stable under the S_m-action permuting the first
m coordinates  IF AND ONLY IF  every nonempty column D_j equals
{m-|D_j|+1, ..., m}  (is "bottom-justified" in [m]).

If true, the only symmetric dual characters of flagged Weyl modules are the
straight-shape Schur polynomials s_lamhat(x_1..x_m), lamhat' = sorted column
sizes.

Exhaustive over all multisets of nonempty columns.  Column ORDER is irrelevant
(supp is a sumset), so multisets suffice and the enumeration is complete.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
from fms import supp_chi, active_rows, is_Sm_stable, is_bottom_justified, lam_hat_from_columns


def all_nonempty_subsets(m):
    out = []
    for r in range(1, m + 1):
        out.extend(itertools.combinations(range(1, m + 1), r))
    return out


def sort_part(a):
    return tuple(sorted((x for x in a if x > 0), reverse=True))


def dominates(lam, mu):
    if sum(lam) != sum(mu):
        return False
    s1 = s2 = 0
    for i in range(max(len(lam), len(mu))):
        s1 += lam[i] if i < len(lam) else 0
        s2 += mu[i] if i < len(mu) else 0
        if s2 > s1:
            return False
    return True


def partitions_le(N, maxlen):
    out = []
    def rec(rem, mx, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        if len(cur) == maxlen: return
        for p in range(min(rem, mx), 0, -1):
            rec(rem - p, p, cur + [p])
    rec(N, N, [])
    return out


def rearrangements(p, m):
    """#distinct arrangements of the multiset p padded with 0s to length m."""
    from math import factorial
    from collections import Counter
    c = Counter(list(p) + [0] * (m - len(p)))
    r = factorial(m)
    for v in c.values():
        r //= factorial(v)
    return r


def run(mmax, kmax, verbose=True):
    tot = sym = 0
    fail_fwd = []   # symmetric but some column not bottom-justified
    fail_bwd = []   # all bottom-justified but not symmetric
    fail_schur = [] # symmetric but supp != {alpha : sort(alpha) <= lamhat}
    for m in range(1, mmax + 1):
        subs = all_nonempty_subsets(m)
        for k in range(1, kmax + 1):
            for D in itertools.combinations_with_replacement(subs, k):
                if active_rows(D) != m:      # m must be the true top row
                    continue
                tot += 1
                n = m                        # all columns inside [m]; extra vars inert
                S = supp_chi(list(D), n)
                stable = is_Sm_stable(S, m)
                bj = all(is_bottom_justified(Dj, m) for Dj in D)
                if stable:
                    sym += 1
                    if not bj:
                        fail_fwd.append(D)
                    lamhat = lam_hat_from_columns(D)
                    N = sum(len(x) for x in D)
                    # (a) every point of S is dominated by lamhat
                    ok = all(dominates(lamhat, sort_part(a)) for a in S)
                    # (b) |S| = #{alpha in N^m : |alpha|=N, sort(alpha) <= lamhat}
                    cnt = 0
                    for p in partitions_le(N, m):
                        if dominates(lamhat, p):
                            cnt += rearrangements(p, m)
                    if not ok or len(S) != cnt:
                        fail_schur.append((D, lamhat, len(S), cnt))
                elif bj:
                    fail_bwd.append(D)
    if verbose:
        print(f"m<={mmax}, #columns<={kmax}:  {tot} diagrams, {sym} with S_m-stable support")
        print(f"  symmetric but NOT all columns bottom-justified : {len(fail_fwd)}")
        print(f"  all columns bottom-justified but NOT symmetric : {len(fail_bwd)}")
        print(f"  symmetric but supp != permutahedron of lamhat  : {len(fail_schur)}")
        for D in fail_fwd[:5]: print("   FWD FAIL", D)
        for D in fail_bwd[:5]: print("   BWD FAIL", D)
        for D in fail_schur[:5]: print("   SCHUR FAIL", D)
    return tot, sym, fail_fwd, fail_bwd, fail_schur


if __name__ == "__main__":
    for mmax, kmax in [(4, 5), (5, 4), (6, 3)]:
        run(mmax, kmax)
        print()
