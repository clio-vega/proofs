#!/usr/bin/env python3
"""Differential check for TworowD4Kernel/GreedyChain.lean.

The Lean file transcribes def:shape, eq:hstrip and eq:greedy of
projects/proofs/2026-09-20-c1-cylindric-M-convexity.tex.  The type checker cannot
detect a wrong transcription, so this script re-implements the SAME definitions
independently, by brute-force enumeration of chains, and checks the four statements
the Lean file proves:

  A1  greedy_isCylindric        : each g^t is a cylindric shape
  A2  greedy_hstrip             : g^t  -hstrip->  g^{t+1}, and g^t subset lambda
  B   greedy_dominates          : x^t_i <= g^t_i for EVERY chain
  C   greedy_eq_lam_of_le       : g^s = lambda for all s >= l0
      nonzeroIncr_indep_of_len  : multiset of nonzero increments is l-independent

It also checks prop:noell's other half -- T_l nonempty iff l >= l0 -- which is NOT
formalised in Lean (it needs lem:greedy(2)).

A shape x : Z -> Z of type (n,m) is stored as the tuple (x_1,...,x_m) with
x_1 < ... < x_m < x_1 + n, extended by x_{i+m} = x_i + n.
"""

from itertools import combinations
import sys

# ---------------------------------------------------------------- the model


def get(v, n, m, i):
    """x_i for arbitrary i in Z, from the fundamental tuple v = (x_1,...,x_m)."""
    q, r = divmod(i - 1, m)
    return v[r] + q * n


def is_cylindric(v, n, m):
    """def:shape: strictly increasing with x_{i+m} = x_i + n.

    Periodicity is automatic from `get`; strict increase on all of Z reduces to
    x_1 < ... < x_m < x_1 + n."""
    return all(v[i] < v[i + 1] for i in range(m - 1)) and v[m - 1] < v[0] + n


def sub(p, q, n, m):
    """p subset q  <=>  p_i <= q_i for i = 1..m (periodicity does the rest)."""
    return all(p[i] <= q[i] for i in range(m))


def hstrip(p, q, n, m):
    """eq:hstrip: p_i <= q_i < p_{i+1} for all i; i = 1..m suffices."""
    return all(p[i] <= q[i] < get(p, n, m, i + 2) for i in range(m))


def greedy_step(lam, g, n, m):
    """eq:greedy: g^t_i = min(g^{t-1}_{i+1} - 1, lambda_i)."""
    return tuple(min(get(g, n, m, i + 2) - 1, lam[i]) for i in range(m))


def u(v):
    return sum(v)


# ---------------------------------------------------------------- enumeration


def shapes(n, m, lo, hi):
    """All cylindric shapes with x_1 in [lo, hi]."""
    out = []
    for x1 in range(lo, hi + 1):
        for rest in combinations(range(x1 + 1, x1 + n), m - 1):
            v = (x1,) + rest
            if is_cylindric(v, n, m):
                out.append(v)
    return out


def successors(p, lam, n, m):
    """All cylindric y with p -hstrip-> y and y subset lambda.

    Built from eq:hstrip DIRECTLY (a box p_i <= y_i <= p_{i+1}-1, intersected with
    y_i <= lam_i), then filtered for cylindricity -- lem:box is NOT assumed."""
    ranges = []
    for i in range(m):
        lo, hi = p[i], min(get(p, n, m, i + 2) - 1, lam[i])
        if lo > hi:
            return []
        ranges.append(range(lo, hi + 1))
    out = []

    def rec(i, acc):
        if i == m:
            y = tuple(acc)
            if is_cylindric(y, n, m) and hstrip(p, y, n, m):
                out.append(y)
            return
        for val in ranges[i]:
            rec(i + 1, acc + [val])

    rec(0, [])
    return out


def all_chains(mu, lam, ell, n, m):
    """Every chain mu = x^0 -hstrip-> ... -hstrip-> x^ell = lam."""
    chains = []

    def rec(path, t):
        if t == ell:
            if path[-1] == lam:
                chains.append(tuple(path))
            return
        for y in successors(path[-1], lam, n, m):
            rec(path + [y], t + 1)

    rec([mu], 0)
    return chains


# ---------------------------------------------------------------- the checks


def main():
    stats = dict(pairs=0, a1=0, a2=0, b_chains=0, b_coords=0, c_fix=0, c_indep=0,
                 nonempty=0, fail=0)
    MAXELL = 6          # enumerate chains up to this length
    EXTRA = 3           # check l-independence for l0 .. l0+EXTRA

    for n in range(2, 8):
        for m in range(1, min(4, n)):
            # mu normalised to mu_1 = 0 (everything is translation-covariant);
            # lam_1 ranges over 0..n so that |lam/mu| stays small enough to enumerate.
            for mu in shapes(n, m, 0, 0):
                for lam in shapes(n, m, 0, n):
                    if not sub(mu, lam, n, m):
                        continue
                    stats["pairs"] += 1

                    # --- build the greedy chain (eq:greedy), l-free
                    g = [mu]
                    for _ in range(MAXELL + EXTRA + 3):
                        g.append(greedy_step(lam, g[-1], n, m))

                    # --- A1: every g^t is a cylindric shape
                    for gt in g:
                        if not is_cylindric(gt, n, m):
                            print("A1 FAIL", n, m, mu, lam, gt); stats["fail"] += 1
                        else:
                            stats["a1"] += 1

                    # --- A2: g^t subset lam, and g^t -hstrip-> g^{t+1}
                    for t in range(len(g) - 1):
                        ok = sub(g[t], lam, n, m) and hstrip(g[t], g[t + 1], n, m)
                        if not ok:
                            print("A2 FAIL", n, m, mu, lam, t); stats["fail"] += 1
                        else:
                            stats["a2"] += 1

                    # --- l0 = min{t : g^t = lam}, if reached
                    # searched only up to MAXELL + 1, so that l0 + EXTRA + 1 stays
                    # inside g; l0 not found therefore means l0 > MAXELL, which is all
                    # the nonemptiness check below needs.
                    l0 = next((t for t, gt in enumerate(g[:MAXELL + 2]) if gt == lam),
                              None)

                    # --- C2: g^s = lam for all s >= l0
                    if l0 is not None:
                        for s in range(l0, len(g)):
                            if g[s] != lam:
                                print("C2 FAIL", n, m, mu, lam, s); stats["fail"] += 1
                            else:
                                stats["c_fix"] += 1

                    # --- B and prop:noell nonemptiness, by full enumeration
                    for ell in range(0, MAXELL + 1):
                        chains = all_chains(mu, lam, ell, n, m)

                        # prop:noell: T_l nonempty  <=>  l >= l0
                        want = (l0 is not None and ell >= l0)
                        if bool(chains) != want:
                            print("NONEMPTY FAIL", n, m, mu, lam, ell, len(chains), l0)
                            stats["fail"] += 1
                        else:
                            stats["nonempty"] += 1

                        # B: greedy dominates coordinatewise
                        for ch in chains:
                            stats["b_chains"] += 1
                            for t in range(ell + 1):
                                for i in range(m):
                                    if ch[t][i] > g[t][i]:
                                        print("B FAIL", n, m, mu, lam, ell, t, i)
                                        stats["fail"] += 1
                                    else:
                                        stats["b_coords"] += 1

                    # --- C3: multiset of nonzero increments is l-independent
                    if l0 is not None:
                        ref = None
                        for ell in range(l0, l0 + EXTRA + 1):
                            gam = [u(g[t + 1]) - u(g[t]) for t in range(ell)]
                            nz = tuple(sorted((x for x in gam if x != 0), reverse=True))
                            if ref is None:
                                ref = nz
                            elif nz != ref:
                                print("C3 FAIL", n, m, mu, lam, ell, nz, ref)
                                stats["fail"] += 1
                            else:
                                stats["c_indep"] += 1

    print()
    for k, v in stats.items():
        print(f"{k:12s} {v}")
    return 1 if stats["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
