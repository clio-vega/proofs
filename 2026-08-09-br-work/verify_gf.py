"""
Verify the generating-function identity:
    Frob(M_{n,r}^{(k)}) = sum_{j=0}^k C(k,j) (-1)^{k-j} F_{r(j+1)-1}
where
    F_m = [deg n] E-bar^m = sum_mu eps(mu) m^{l(mu)} / z_mu * p_mu

and hence
    span_{r>=1, 0<=k<=n} { Frob(M_{n,r}^{(k)}) }  =  span { F_m : m >= 1 }
                                                   =  span { P_ell : 1 <= ell <= n }
where P_ell = sum_{l(mu)=ell} eps(mu)/z_mu * p_mu.

We verify by:
  (a) Computing Frob(M_{6,3}^{(k)}) via probe.py for k = 0..6.
  (b) Computing sum_{j=0}^k C(k,j)(-1)^{k-j} F_{r(j+1)-1} in power-sum basis via F_m formula.
  (c) Verifying (a) = (b).
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, '/home/clio/projects/probes/2026-08-09-bhattacharya-rhoades-wreath-superspace')
from probe import frobenius_op_tilde_slice, sign_twist
from symmfunc import partitions_of, chi, z_lambda, p_to_s, s_to_p, print_s
from sympy import S, Matrix, Rational, Integer, expand, binomial
from collections import defaultdict


def epsilon(mu):
    """Sign of any permutation of cycle-type mu: (-1)^(n - l(mu))"""
    return (-1) ** (sum(mu) - len(mu))


def F_m_powersum(n, m):
    """Return F_m = [deg n] Ebar^m in power-sum basis as dict {mu -> coeff}."""
    out = {}
    for mu in partitions_of(n):
        c = Rational(epsilon(mu) * m ** len(mu), z_lambda(mu))
        if c != 0:
            out[mu] = c
    return out


def F_m_schur(n, m):
    """Return F_m in Schur basis."""
    p = F_m_powersum(n, m)
    return p_to_s(p, n)


def combine_p_dicts(dicts_with_coeffs):
    """Given list of (coeff, dict), sum them."""
    out = defaultdict(lambda: S.Zero)
    for c, d in dicts_with_coeffs:
        for k, v in d.items():
            out[k] += c * v
    return {k: expand(v) for k, v in out.items() if expand(v) != 0}


def Frob_slice_from_formula(n, r, k):
    """Frob(M_{n,r}^{(k)}) via the binomial formula."""
    terms = []
    for j in range(k + 1):
        coef = Integer((-1) ** (k - j)) * binomial(k, j)
        terms.append((coef, F_m_powersum(n, r * (j + 1) - 1)))
    return combine_p_dicts(terms)


def dicts_equal(d1, d2):
    keys = set(d1.keys()) | set(d2.keys())
    for k in keys:
        if expand(d1.get(k, 0) - d2.get(k, 0)) != 0:
            return False, k
    return True, None


def main():
    n = 6
    r = 3
    print(f"=== VERIFY: n = {n}, r = {r} ===")
    print()
    print(f"Expected: Frob(M_{{6,3}}^{{(k)}}) via probe.py  ==  formula via F_m")
    print()

    for k in range(n + 1):
        # (a) via probe
        raw = frobenius_op_tilde_slice(n, r, k)
        via_probe_s = sign_twist(raw)  # in Schur basis
        via_probe_p = s_to_p(via_probe_s, n)

        # (b) via formula
        via_formula_p = Frob_slice_from_formula(n, r, k)

        ok, bad_key = dicts_equal(via_probe_p, via_formula_p)
        print(f"  k = {k}: {'MATCH' if ok else f'MISMATCH at {bad_key}'}")
        if not ok:
            print(f"    probe: {via_probe_p}")
            print(f"    formula: {via_formula_p}")
            print(f"    diff at {bad_key}: probe={via_probe_p.get(bad_key,0)}, formula={via_formula_p.get(bad_key,0)}")


if __name__ == "__main__":
    main()
