"""
Verification for 2026-09-25 c2: SNP and Newton(s^c) = P_lambdahat,
proved WITHOUT citing Rado, from Lemma 4.1 (Robin Hood step) of
proofs/2026-09-20-c1-cylindric-M-convexity.tex.

Three checks, all exact (Fraction / integer arithmetic):

 V1  J = J'      sorted form vs subset form            (Lemma A of the writeup)
 V2  the CHAIN   Lemma 4.1 iterated from lambdahat down to sigma:
                 - every step has nu_a - nu_b >= 2   (so t = 1/(nu_a-nu_b) in (0,1/2])
                 - the statistic Phi = sum_r (N_r - S_r) strictly decreases
                 - the UNROLLED convex combination of PERMUTATIONS of lambdahat
                   reproduces alpha exactly, with coefficients >= 0 summing to 1
 V3  T1          conv(W) cap Z^l = W on the subset description

V2 is the point: it tests the route, not the conclusion.  A pure LP
feasibility check would confirm alpha in P_lambdahat while saying nothing
about whether MY argument produces it.
"""
from fractions import Fraction
from itertools import permutations
import sys
from collections import defaultdict

# ---------------------------------------------------------------- basics

def sort_desc(a):
    return tuple(sorted(a, reverse=True))

def psums(v):
    out, s = [], 0
    for x in v:
        s += x
        out.append(s)
    return out

def dominates(nu, sig):
    """nu |> sig (weakly): same size, partial sums of nu >= those of sig."""
    assert sum(nu) == sum(sig)
    N, S = psums(nu), psums(sig)
    return all(n >= s for n, s in zip(N, S))

def compositions(d, ell):
    """all alpha in N^ell with |alpha| = d"""
    if ell == 1:
        yield (d,); return
    for first in range(d + 1):
        for rest in compositions(d - first, ell - 1):
            yield (first,) + rest

def partitions_le(d, ell):
    """partitions of d with at most ell parts, padded to length ell"""
    seen = set()
    for a in compositions(d, ell):
        seen.add(sort_desc(a))
    return sorted(seen, reverse=True)

# ---------------------------------------------------------------- V1

def in_J_sorted(alpha, lamhat):
    """sort(alpha) <| lamhat  (dominance, which presupposes equal size)"""
    if sum(alpha) != sum(lamhat):
        return False
    return dominates(lamhat, sort_desc(alpha))

def in_J_subset(alpha, lamhat):
    """alpha([l]) = d  and  alpha(S) <= Lambda_{|S|} for all S subseteq [l].
       Tested over ALL subsets, not just the greedy ones -- no shortcut."""
    ell = len(alpha)
    d = sum(lamhat)
    if sum(alpha) != d:
        return False
    Lam = [0] + psums(sort_desc(lamhat))
    for mask in range(1 << ell):
        S = [i for i in range(ell) if mask >> i & 1]
        if sum(alpha[i] for i in S) > Lam[len(S)]:
            return False
    return True

def V1(dmax=9, ellmax=5):
    pairs = 0; pts = 0; bad = 0
    for ell in range(1, ellmax + 1):
        for d in range(0, dmax + 1):
            for lamhat in partitions_le(d, ell):
                pairs += 1
                for alpha in compositions(d, ell):
                    pts += 1
                    if in_J_sorted(alpha, lamhat) != in_J_subset(alpha, lamhat):
                        bad += 1
                        if bad <= 5:
                            print("  V1 MISMATCH", lamhat, alpha)
    print(f"V1  J = J' : {pairs} pairs (lamhat,ell), {pts} points tested, {bad} disagreements")
    return bad == 0

# ---------------------------------------------------------------- V2

def robin_hood(nu, sig):
    """Lemma 4.1 verbatim: returns (a, b) with a<b, nu_a >= nu_b + 2,
       tau = nu - e_a + e_b a partition, nu |> tau |>= sig.
       nu, sig padded to the same length ell.  Requires nu |> sig strictly."""
    ell = len(nu)
    N, S = psums(nu), psums(sig)
    # i minimal with nu_i > sig_i
    i = next(c for c in range(ell) if nu[c] > sig[c])
    # j minimal > i with nu_j < sig_j
    j = next(c for c in range(i + 1, ell) if nu[c] < sig[c])
    a = max(c for c in range(i, j) if nu[c] == nu[i])
    b = min(c for c in range(i + 1, j + 1) if nu[c] == nu[j])
    return a, b

def chain(lamhat, sigma):
    """Iterate Lemma 4.1 from lamhat down to sigma.
       Returns the list of (a,b,gap) steps.  Checks the two proof obligations:
       gap >= 2 at every step, and Phi strictly decreasing."""
    nu = tuple(lamhat)
    steps = []
    Phi = sum(n - s for n, s in zip(psums(nu), psums(sigma)))
    guard = 0
    while nu != sigma:
        assert dominates(nu, sigma), "invariant nu |>= sigma lost"
        a, b = robin_hood(nu, sigma)
        gap = nu[a] - nu[b]
        assert a < b, f"a<b failed: {a},{b}"
        assert gap >= 2, f"GAP < 2 at {nu} -> a={a},b={b},gap={gap}"   # obligation 1
        tau = list(nu); tau[a] -= 1; tau[b] += 1; tau = tuple(tau)
        assert all(tau[k] >= tau[k+1] for k in range(len(tau)-1)), f"tau not a partition: {tau}"
        assert all(x >= 0 for x in tau), f"tau not in N^l: {tau}"
        Phi2 = sum(n - s for n, s in zip(psums(tau), psums(sigma)))
        assert Phi2 == Phi - (b - a), f"Phi bookkeeping: {Phi}->{Phi2}, b-a={b-a}"
        assert Phi2 < Phi, "Phi did not strictly decrease"           # obligation 2
        Phi = Phi2
        steps.append((a, b, gap))
        nu = tau
        guard += 1
        assert guard < 10000
    return steps

def unroll(lamhat, sigma):
    """Build the EXPLICIT convex combination of permutations of lamhat equal to
       sigma, by unrolling the induction:  nu^{t+1} = (1-s) nu^t + s (a b) nu^t
       with s = 1/(nu^t_a - nu^t_b).  Represent each nu^t as a dict
       {permutation-image-of-lamhat : coefficient}.  This is obligation 3:
       the endpoints reached are permutations of lamhat, not merely points of J."""
    ell = len(lamhat)
    combo = {tuple(lamhat): Fraction(1)}
    nu = tuple(lamhat)
    while nu != sigma:
        a, b = robin_hood(nu, sigma)
        s = Fraction(1, nu[a] - nu[b])
        new = defaultdict(Fraction)
        for pt, c in combo.items():
            sw = list(pt); sw[a], sw[b] = sw[b], sw[a]
            new[pt] += (1 - s) * c
            new[tuple(sw)] += s * c
        combo = {k: v for k, v in new.items() if v != 0}
        tau = list(nu); tau[a] -= 1; tau[b] += 1
        nu = tuple(tau)
    return combo

def V2(dmax=9, ellmax=4):
    tested = 0; steps_tot = 0; bad = 0; maxsupport = 0
    for ell in range(1, ellmax + 1):
        for d in range(0, dmax + 1):
            for lamhat in partitions_le(d, ell):
                for sigma in partitions_le(d, ell):
                    if not dominates(lamhat, sigma):
                        continue
                    tested += 1
                    try:
                        st = chain(lamhat, sigma)
                    except AssertionError as e:
                        bad += 1
                        if bad <= 5: print("  V2 CHAIN FAIL", lamhat, sigma, e)
                        continue
                    steps_tot += len(st)
                    combo = unroll(lamhat, sigma)
                    maxsupport = max(maxsupport, len(combo))
                    # coefficients: nonnegative, sum to 1
                    if any(c < 0 for c in combo.values()) or sum(combo.values()) != 1:
                        bad += 1; print("  V2 NOT CONVEX", lamhat, sigma); continue
                    # every point in the support is a permutation of lamhat
                    if any(sort_desc(p) != sort_desc(lamhat) for p in combo):
                        bad += 1; print("  V2 NOT A PERMUTATION OF lamhat", lamhat, sigma); continue
                    # the combination reproduces sigma exactly
                    got = tuple(sum(c * p[k] for p, c in combo.items()) for k in range(ell))
                    if got != tuple(Fraction(x) for x in sigma):
                        bad += 1; print("  V2 WRONG POINT", lamhat, sigma, got); continue
    print(f"V2  chain+unroll : {tested} pairs (lamhat,sigma) with sigma <| lamhat, "
          f"{steps_tot} Robin Hood steps, max |support| = {maxsupport}, {bad} failures")
    return bad == 0

# ---------------------------------------------------------------- V3

def V3(dmax=8, ellmax=4):
    """SNP on the subset description: every lattice point of conv(J) is in J.
       conv(J) is contained in Q = {x >= 0, x([l]) = d, x(S) <= Lambda_{|S|}};
       so it suffices to check Q cap Z^l = J, which is what in_J_subset says.
       The non-trivial content: enumerate ALL lattice points of the bounding box
       satisfying Q and confirm each is in J (i.e. sort <| lamhat)."""
    bad = 0; pts = 0; pairs = 0
    for ell in range(1, ellmax + 1):
        for d in range(0, dmax + 1):
            for lamhat in partitions_le(d, ell):
                pairs += 1
                for alpha in compositions(d, ell):
                    if in_J_subset(alpha, lamhat):
                        pts += 1
                        if not in_J_sorted(alpha, lamhat):
                            bad += 1
    print(f"V3  Q cap Z^l = J : {pairs} pairs, {pts} lattice points of Q, {bad} outside J")
    return bad == 0

if __name__ == "__main__":
    ok = True
    ok &= V1(dmax=9, ellmax=5)
    ok &= V2(dmax=9, ellmax=4)
    ok &= V3(dmax=8, ellmax=4)
    print("ALL PASS" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)
