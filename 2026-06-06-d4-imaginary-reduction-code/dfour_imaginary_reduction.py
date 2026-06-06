"""
d=4 imaginary-part reduction: does Im G_lambda(i) generalise to all shapes? (IR-4)
=================================================================================

Object:  G_lambda(i) = <s_lambda, psi^m>,   psi = h_2 + i e_2 = s_2 + i s_{11},  lambda |- 2m.

Central identity (to verify computationally):

    G_lambda(i) = sum_{k=0}^m  C(m,k) i^k  N_{lambda,k},
        N_{lambda,k} := <s_lambda, h_2^{m-k} e_2^k> = <s_lambda, s_2^{m-k} s_{11}^k>  >= 0.

Hence
    Re G_lambda = sum_{k even} (-1)^{k/2}     C(m,k) N_{lambda,k}
    Im G_lambda = sum_{k odd}  (-1)^{(k-1)/2} C(m,k) N_{lambda,k}      (*)

(*) is ALREADY a signed sum of nonnegative integers for EVERY lambda.  The two-row
trinomial bracket is the rank-2 instance N_{(2m-b,b),k} = T(m-k,b-k) - T(m-k,b-k-1).

This module:
  1. computes N_{lambda,k} by iterated Pieri in a fixed order (h_2's then e_2's),
  2. cross-checks the identity (*) against the direct/chain G_lambda(i),
  3. cross-checks N_{lambda,k} = vdist[k] / C(m,k) (commutativity of the chain model),
  4. tabulates the zero-set of Im G_lambda over all lambda |- 2m,
  5. tabulates the lower bound min_{lambda != (2,2)} |Im G_lambda|.

No SageMath required -- pure-Python Pieri rules (reused from the involution scratch core).
"""
import sys, os
from math import comb
from functools import lru_cache

sys.path.insert(0, os.path.dirname(__file__))  # local core.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                "..", "scratch", "2026-06-09-d4-involution"))
from core import (partitions, horizontal_strip_adds, vertical_strip_adds,
                  psi_power_schur, G_chains)

@lru_cache(maxsize=None)
def psi_pow(m):
    """psi^m in Schur basis (computed ONCE per m, not per shape)."""
    return psi_power_schur(m)

def G_direct_all(m):
    """{lambda: (Re, Im)} for all lambda |- 2m, from a single psi^m expansion."""
    return {lam: rim for lam, rim in psi_pow(m).items()}

# ---------- integer Schur algebra (no i) ----------

def pieri_h2(f):
    """Multiply an integer Schur-dict by h_2 = s_2 (horizontal 2-strips)."""
    out = {}
    for mu, c in f.items():
        for nu in horizontal_strip_adds(mu, 2):
            out[nu] = out.get(nu, 0) + c
    return out

def pieri_e2(f):
    """Multiply an integer Schur-dict by e_2 = s_{11} (vertical 2-strips)."""
    out = {}
    for mu, c in f.items():
        for nu in vertical_strip_adds(mu, 2):
            out[nu] = out.get(nu, 0) + c
    return out

@lru_cache(maxsize=None)
def Nk_table(m):
    """Return dict  k -> {lambda: N_{lambda,k}}  for all lambda |- 2m, k=0..m.

    N_{lambda,k} = <s_lambda, s_2^{m-k} s_{11}^k>, computed by applying e_2 k times
    then h_2 (m-k) times to the empty partition (fixed order; commutativity makes
    the order irrelevant to the final coefficient)."""
    table = {}
    for k in range(m + 1):
        f = {(): 1}
        for _ in range(k):
            f = pieri_e2(f)
        for _ in range(m - k):
            f = pieri_h2(f)
        table[k] = f  # {lambda: coeff}, all lambda |- 2m
    return table

def ReIm_from_Nk(lam, m):
    """(Re, Im) of G_lambda(i) from the signed-bracket decomposition (*)."""
    table = Nk_table(m)
    lam = tuple(lam)
    Re = Im = 0
    for k in range(m + 1):
        N = table[k].get(lam, 0)
        if N == 0:
            continue
        coeff = comb(m, k) * N
        r = k % 4
        if r == 0:
            Re += coeff
        elif r == 1:
            Im += coeff
        elif r == 2:
            Re -= coeff
        else:
            Im -= coeff
    return Re, Im

# ---------- verification suite ----------

def verify(mmax_chain=6, mmax_direct=11, verbose=True):
    """Cross-check (*) for all lambda |- 2m.

    - chain model (G_chains, exponential) up to m<=mmax_chain: full three-way check
      direct == chain == bracket, plus vdist[k] == C(m,k) N_{lam,k} and N >= 0;
    - direct (G_direct, symbolic psi^m) up to m<=mmax_direct: bracket == direct, N >= 0.
    """
    all_ok = True
    for m in range(1, mmax_direct + 1):
        n = 2 * m
        table = Nk_table(m)
        do_chain = m <= mmax_chain
        direct_all = G_direct_all(m)                     # one psi^m expansion
        for lam in partitions(n):
            ReIm_direct = direct_all.get(lam, (0, 0))    # symbolic psi^m
            ReIm_bracket = ReIm_from_Nk(lam, m)          # identity (*)
            if ReIm_direct != ReIm_bracket:
                all_ok = False
                print(f"  MISMATCH m={m} lam={lam}: direct={ReIm_direct} "
                      f"bracket={ReIm_bracket}")
            for k in range(m + 1):                       # N >= 0
                if table[k].get(lam, 0) < 0:
                    all_ok = False
                    print(f"  NEGATIVE N m={m} lam={lam} k={k}: {table[k][lam]}")
            if do_chain:
                (a, b), vdist = G_chains(lam)           # chain model
                if (a, b) != ReIm_direct:
                    all_ok = False
                    print(f"  CHAIN MISMATCH m={m} lam={lam}: chain={(a,b)} "
                          f"direct={ReIm_direct}")
                for k in range(m + 1):                   # vdist[k] = C(m,k) N
                    lhs = vdist.get(k, 0)
                    rhs = comb(m, k) * table[k].get(lam, 0)
                    if lhs != rhs:
                        all_ok = False
                        print(f"  vdist!=C*N m={m} lam={lam} k={k}: {lhs} vs {rhs}")
        if verbose:
            tag = "(+chain)" if do_chain else "(direct+bracket)"
            print(f"  m={m:2d} (n={n:2d}): {len(partitions(n)):4d} shapes {tag}",
                  flush=True)
    return all_ok


def two_row_bracket_check(mmax=10):
    """Verify N_{(2m-b,b),k} = T(m-k, b-k) - T(m-k, b-k-1), the two-row formula,
    where T(p, j) = [x^j] (1+x+x^2)^p is the central trinomial coefficient."""
    @lru_cache(maxsize=None)
    def trinomial(p, j):
        # [x^j] (1+x+x^2)^p ; j may be out of range -> 0
        if j < 0 or j > 2 * p:
            return 0
        # T(p,j) = sum_t C(p,t) C(p-t, j-2t)
        s = 0
        for t in range(0, p + 1):
            s += comb(p, t) * comb(p - t, j - 2 * t) if 0 <= j - 2 * t <= p - t else 0
        return s
    ok = True
    for m in range(1, mmax + 1):
        table = Nk_table(m)
        for b in range(0, m + 1):           # two-row (2m-b, b), need 2m-b >= b
            lam = (2 * m - b, b) if b > 0 else (2 * m,)
            for k in range(m + 1):
                N = table[k].get(lam, 0)
                formula = trinomial(m - k, b - k) - trinomial(m - k, b - k - 1)
                if N != formula:
                    ok = False
                    print(f"  TWO-ROW MISMATCH m={m} b={b} k={k}: N={N} formula={formula}")
    return ok


def zero_set(mmax=12):
    """All lambda |- 2m (m<=mmax) with Im G_lambda = 0.  Returns dict m -> list."""
    result = {}
    for m in range(1, mmax + 1):
        zeros = []
        for lam in partitions(2 * m):
            Re, Im = ReIm_from_Nk(lam, m)
            if Im == 0:
                zeros.append((lam, Re))
        result[m] = zeros
    return result


def im_lower_bound(mmax=12):
    """min_{lambda |- 2m, lambda != (2,2)} |Im G_lambda|, with a minimiser."""
    rows = []
    for m in range(1, mmax + 1):
        best = None
        for lam in partitions(2 * m):
            if lam == (2, 2):
                continue
            Re, Im = ReIm_from_Nk(lam, m)
            a = abs(Im)
            if best is None or a < best[0]:
                best = (a, lam)
        rows.append((m, best[0], best[1]))
    return rows


if __name__ == "__main__":
    print("=" * 70)
    print("TASK 1: verify identity  Im G = sum_{k odd} (-1)^((k-1)/2) C(m,k) N_{lam,k}")
    print("        (three-way cross-check; N >= 0; vdist[k] = C(m,k) N_{lam,k})")
    print("=" * 70)
    ok = verify(mmax_chain=6, mmax_direct=11)
    print(f"\n  ALL CHECKS PASS (chain m<=6, direct+bracket m<=11): {ok}\n", flush=True)

    print("=" * 70)
    print("        two-row trinomial-bracket specialisation check")
    print("=" * 70)
    ok2 = two_row_bracket_check(mmax=10)
    print(f"  two-row N = T(m-k,b-k) - T(m-k,b-k-1)  PASS (m<=10): {ok2}\n")

    print("=" * 70)
    print("TASK 2: ZERO-SET of Im G_lambda over all lambda |- 2m, m <= 12")
    print("=" * 70)
    zs = zero_set(mmax=12)
    total_zeros = 0
    for m, zeros in zs.items():
        total_zeros += len(zeros)
        tag = "  <-- (2,2) only" if zeros == [((2, 2), 0)] else ""
        if zeros:
            print(f"  m={m:2d} (n={2*m:2d}): {len(zeros)} zero(s): "
                  f"{[(lam, f'Re={Re}') for lam, Re in zeros]}{tag}")
        else:
            print(f"  m={m:2d} (n={2*m:2d}): NO Im-zeros")
    print(f"\n  total Im-zeros across m<=12: {total_zeros}")

    print("\n" + "=" * 70)
    print("TASK 3: lower bound  min_{lambda != (2,2)} |Im G_lambda|  vs m")
    print("=" * 70)
    print(f"  {'m':>3} {'min|Im|':>8}  minimiser")
    for m, mn, lam in im_lower_bound(mmax=12):
        print(f"  {m:>3} {mn:>8}  {lam}")
