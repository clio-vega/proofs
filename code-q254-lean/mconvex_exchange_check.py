"""Differential check for the M-convex exchange axiom on
    J(lhat) = { alpha in N^ell : sort(alpha) <= lhat in dominance }.

Written from the DEFINITION, independently of the Lean file
TworowD4Kernel/MConvexExchange.lean.  Three things are checked.

  (EQ)  the two characterisations of J agree:
          sorted-dominance form   sort(alpha) |> dominated by lhat
          subset form             alpha(S) <= Lambda_{|S|} for all S, alpha([ell]) = d
        This equality is the first sentence of the proof of prop:perm-mconvex in
        proofs/2026-09-20-c1-cylindric-M-convexity.tex and is the ONE step the Lean
        file does not formalise (Lean works with the subset form throughout).

  (EX1) one-sided exchange: for all alpha,beta in J and i with beta_i < alpha_i,
        there is j with alpha_j < beta_j and alpha - e_i + e_j in J.
        This is exactly what the paper's proof establishes, and exactly what the
        Lean theorem `insupp_exchange_onesided` states.

  (EX2) two-sided (symmetric) exchange -- Murota's axiom for M-convex sets: the SAME j
        must additionally satisfy beta + e_i - e_j in J.  The paper asserts M-convexity
        (= EX2) but its proof only delivers EX1; see the session note.
"""
import itertools, sys
from collections import Counter


def partitions(d, maxlen, maxpart=None):
    """Partitions of d with at most maxlen parts, as tuples (weakly decreasing)."""
    if maxpart is None:
        maxpart = d
    if d == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for first in range(min(d, maxpart), 0, -1):
        for rest in partitions(d - first, maxlen - 1, first):
            yield (first,) + rest


def pad(nu, ell):
    return tuple(nu) + (0,) * (ell - len(nu))


def dominated(sigma, lhat, ell):
    """sigma <= lhat in dominance order (both padded to length ell, equal sums)."""
    s = t = 0
    for r in range(ell):
        s += sigma[r]
        t += lhat[r]
        if s > t:
            return False
    return True


def J_sorted(lhat, ell):
    """{alpha in N^ell : sum alpha = |lhat|, sort_desc(alpha) dominated by lhat}."""
    d = sum(lhat)
    L = pad(lhat, ell)
    out = set()
    for alpha in weak_compositions(d, ell):
        if dominated(tuple(sorted(alpha, reverse=True)), L, ell):
            out.add(alpha)
    return out


def J_subset(lhat, ell):
    """{alpha in N^ell : alpha([ell]) = d, alpha(S) <= Lambda_{|S|} for all S}."""
    d = sum(lhat)
    L = pad(lhat, ell)
    Lam = [0] * (ell + 1)
    for r in range(ell):
        Lam[r + 1] = Lam[r] + L[r]
    subsets = [S for k in range(ell + 1) for S in itertools.combinations(range(ell), k)]
    out = set()
    for alpha in weak_compositions(d, ell):
        if sum(alpha) != d:
            continue
        if all(sum(alpha[c] for c in S) <= Lam[len(S)] for S in subsets):
            out.add(alpha)
    return out


def weak_compositions(d, ell):
    if ell == 0:
        if d == 0:
            yield ()
        return
    for first in range(d + 1):
        for rest in weak_compositions(d - first, ell - 1):
            yield (first,) + rest


def step(alpha, i, j):
    a = list(alpha)
    a[i] -= 1
    a[j] += 1
    return tuple(a)


def main(dmax=9, ellmax=4):
    eq_bad = eq_tot = 0
    ex1_bad = ex1_tot = 0
    ex2_bad = ex2_tot = 0
    for ell in range(1, ellmax + 1):
        for d in range(0, dmax + 1):
            for lhat in partitions(d, ell):
                A = J_sorted(lhat, ell)
                B = J_subset(lhat, ell)
                eq_tot += 1
                if A != B:
                    eq_bad += 1
                    print("EQ FAIL", lhat, ell, sorted(A ^ B)[:4])
                Jset = A
                for alpha in Jset:
                    for beta in Jset:
                        for i in range(ell):
                            if not (beta[i] < alpha[i]):
                                continue
                            ex1_tot += 1
                            ex2_tot += 1
                            js = [j for j in range(ell) if alpha[j] < beta[j]]
                            ok1 = any(step(alpha, i, j) in Jset for j in js)
                            ok2 = any(step(alpha, i, j) in Jset
                                      and step(beta, j, i) in Jset for j in js)
                            if not ok1:
                                ex1_bad += 1
                                print("EX1 FAIL", lhat, ell, alpha, beta, i)
                            if not ok2:
                                ex2_bad += 1
                                print("EX2 FAIL", lhat, ell, alpha, beta, i)
    print(f"range: |lhat| <= {dmax}, ell <= {ellmax}")
    print(f"(EQ)  sorted-form vs subset-form : {eq_bad} disagreements / {eq_tot} pairs (lhat, ell)")
    print(f"(EX1) one-sided exchange         : {ex1_bad} failures / {ex1_tot} triples (alpha, beta, i)")
    print(f"(EX2) two-sided exchange         : {ex2_bad} failures / {ex2_tot} triples (alpha, beta, i)")


if __name__ == "__main__":
    main(*(int(x) for x in sys.argv[1:]))
