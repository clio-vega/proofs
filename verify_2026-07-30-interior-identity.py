"""
Verify the interior identity M^{(c)}_lambda = v_lambda(t) * P_lambda for
n in {3, 4}, c in {3, 4, 5}, |lambda| up to a reasonable bound.

APPROACH.  We build M^{(c)}_mu layer by layer using vDEZ Cor 4.2 Pieri,
  DROPPING non-dominant lambda'+e_j targets (which vanish in the polynomial
  interpretation; the reflection contributions in the fusion quotient
  correspond to non-dominant weights being identified with signed dominant
  images, but for our purposes -- computing polynomial reps -- we drop).

For each interior mu, we compare M_mu against v_mu * P_mu and record whether
they are equal as polynomials in Q(t)[x_1, ..., x_n].

For each boundary mu, we record M_mu's P-basis expansion (which is generally
NOT scalar * P_mu; it mixes multiple P's).
"""

import sys, os
import sympy as sp
from sympy import expand, cancel, together, factor, simplify, Poly, Rational, symbols
from itertools import permutations, combinations
from functools import lru_cache

# Add path so we can import existing code
sys.path.insert(0, '/home/clio/projects/scratch/2026-07-24-pieri-affine-CR')
from hall_littlewood_engine import hall_littlewood_P, t_sym, v_lambda, partitions
from demazure_engine import make_vars

t = t_sym

# ---------------------------------------------------------------------------
# Helpers: pairing, positive roots of A_{n-1}, e_hat, V_{lambda, nu}(t)
# ---------------------------------------------------------------------------

def pair(x, y):
    return sum(a*b for a, b in zip(x, y))

def positive_roots(n):
    """R_0^+ of type A_{n-1}: {e_i - e_j : i < j}."""
    R = []
    for i in range(n):
        for j in range(i+1, n):
            v = [0]*n
            v[i] = 1
            v[j] = -1
            R.append(tuple(v))
    return R

def theta_root(n):
    """Highest root theta = e_1 - e_n."""
    v = [0]*n
    v[0] = 1
    v[-1] = -1
    return tuple(v)

def e_hat_t(alpha, R_plus, n):
    """e_hat(alpha) = prod_{beta in R^+} t^{<alpha, beta>/2}.
    Simply-laced A: t^{sum <alpha, beta>/2}.
    For alpha = e_i - e_j, this equals t^{j-i}."""
    total = sp.Rational(0)
    for beta in R_plus:
        total += sp.Rational(pair(alpha, beta), 2)
    return t**total

def V_coeff(lam, nu, c, n):
    """vDEZ Eq (2.11) simplified for simply-laced A_{n-1}, single parameter t."""
    R_plus = positive_roots(n)
    theta_v = theta_root(n)
    result = sp.Integer(1)
    for beta in R_plus:
        p_lam = pair(lam, beta)
        p_nu = pair(nu, beta)
        if p_lam == 0 and p_nu > 0:
            ehat = e_hat_t(beta, R_plus, n)
            result *= (1 - t * ehat) / (1 - ehat)
    # h_hat_t = t * e_hat(theta) = t * t^{n-1} = t^n
    h_hat = t * e_hat_t(theta_v, R_plus, n)
    for beta in R_plus:
        p_lam = pair(lam, beta)
        p_nu = pair(nu, beta)
        if p_lam == c and p_nu < 0:
            ehat_neg = e_hat_t(tuple(-x for x in beta), R_plus, n)
            result *= (1 - t * h_hat * ehat_neg) / (1 - h_hat * ehat_neg)
    return sp.simplify(result)


# ---------------------------------------------------------------------------
# Dominant partitions in P_c^+
# ---------------------------------------------------------------------------

def is_dominant(lam):
    return all(lam[i] >= lam[i+1] for i in range(len(lam)-1)) and lam[-1] >= 0

def in_Pc(lam, c, n):
    if not is_dominant(lam):
        return False
    theta_v = theta_root(n)
    if pair(lam, theta_v) > c:
        return False
    return True

def is_interior(lam, c, n):
    """Interior means <lam, beta> < c for all beta in R_0^+ (equivalently lam_1 - lam_n < c)."""
    if not is_dominant(lam):
        return False
    theta_v = theta_root(n)
    return pair(lam, theta_v) < c

def enumerate_alcove(n, c, max_size):
    """All (lam1, ..., lam_n) dominant, in P_c^+, |lam| <= max_size."""
    parts = []
    for size in range(0, max_size + 1):
        for p in partitions(size, max_len=n):
            p3 = tuple(list(p) + [0]*(n - len(p)))
            if in_Pc(p3, c, n):
                parts.append(p3)
    return sorted(parts, key=lambda p: (sum(p), tuple(-a for a in p)))


# ---------------------------------------------------------------------------
# Ordinary HL Pieri (from Macdonald III (5.7) specialised to horizontal 1-strip)
# We compute it via LINEAR ALGEBRA rather than hard-coding the psi formula.
# ---------------------------------------------------------------------------

def ordinary_HL_pieri(lam, n, X, P_dict):
    """Given lam (partition, tuple of length n with trailing zeros), compute
       m_omega1 * P_lam and decompose as sum over horizontal-1-strip mu:
           m_omega1 * P_lam = sum_mu psi_{mu/lam}(t) * P_mu.
       Returns dict {mu: psi_{mu/lam}}. Uses P_dict for known P's."""
    m_omega1 = sum(X)
    lhs = expand(m_omega1 * P_dict[lam])
    size = sum(lam) + 1
    mu_list = [tuple(list(p) + [0]*(n - len(p))) for p in partitions(size, max_len=n)]
    # ensure we have P_mu for each
    for mu in mu_list:
        if mu not in P_dict:
            P_dict[mu] = expand(hall_littlewood_P(mu, n, X))
    # Solve linear system: lhs = sum_mu psi_mu * P_mu
    def coeff_at(f, part):
        p_ = list(part) + [0]*(n - len(part))
        return Poly(f, X).coeff_monomial(tuple(p_))
    monom_basis = mu_list
    M = sp.Matrix(len(monom_basis), len(mu_list),
                  lambda i, j: coeff_at(P_dict[mu_list[j]], monom_basis[i]))
    b = sp.Matrix(len(monom_basis), 1,
                  lambda i, _: coeff_at(lhs, monom_basis[i]))
    sol = M.solve(b)
    return {mu_list[j]: sp.expand(sol[j]) for j in range(len(mu_list))}


# ---------------------------------------------------------------------------
# Build M layer by layer using drop-non-dominant vDEZ Pieri
# ---------------------------------------------------------------------------

def build_M_iterated(n, c, max_size, verbose=False):
    """Build M^{(c)}_mu using ansatz for interior mu and Pieri from interior
       precursor for boundary mu. See build_M_via_ansatz."""
    X = make_vars(n)
    all_lam = enumerate_alcove(n, c, max_size)
    P_dict = {}
    v_dict = {}
    for lam in all_lam:
        P_dict[lam] = expand(hall_littlewood_P(lam, n, X))
        v_dict[lam] = sp.expand(v_lambda(list(lam), n=n))
    # Also precompute P for out-of-alcove partitions of size <= max_size (they
    # arise when we compute m_omega1 * M_boundary via ordinary Pieri)
    for size in range(0, max_size + 1):
        for p in partitions(size, max_len=n):
            p3 = tuple(list(p) + [0]*(n - len(p)))
            if p3 not in P_dict:
                P_dict[p3] = expand(hall_littlewood_P(p3, n, X))
                v_dict[p3] = sp.expand(v_lambda(list(p3), n=n))
    return build_M_via_ansatz(n, c, max_size, X, P_dict, v_dict, verbose=verbose)


def build_M_via_ansatz(n, c, max_size, X, P_dict, v_dict, verbose=False):
    """Build M's using: interior mu -> ansatz M = v*P; boundary mu -> Pieri from
       an interior precursor with unique boundary target.

       Then VERIFY that vDEZ Pieri holds (as polynomial identity, drop-non-dom)
       from all interior precursors."""
    m_omega1 = sum(X)
    W0_omega1 = []
    for k in range(n):
        v = [0]*n
        v[k] = 1
        W0_omega1.append(tuple(v))

    all_lam = enumerate_alcove(n, c, max_size)
    M_dict = {}
    status_dict = {}

    for lam in all_lam:
        interior = is_interior(lam, c, n)
        status_dict[lam] = 'interior' if interior else 'boundary'
        if interior:
            M_dict[lam] = expand(v_dict[lam] * P_dict[lam])

    # For boundary mu: compute from Pieri equation from ANY interior precursor.
    # We'll use the first one we find.
    for lam in all_lam:
        if lam in M_dict:
            continue
        # Find an interior precursor
        found = False
        for i in range(n):
            lam_p = list(lam)
            lam_p[i] -= 1
            lam_p = tuple(lam_p)
            if not is_dominant(lam_p):
                continue
            if not in_Pc(lam_p, c, n):
                continue
            if not is_interior(lam_p, c, n):
                continue
            # Use this precursor
            LHS_poly = expand(m_omega1 * M_dict[lam_p])
            # RHS = sum_j V_{lam_p, e_j} * M_{lam_p+e_j} over dominant, in P_c^+
            V_target = None
            RHS_others = sp.Integer(0)
            for j in range(n):
                nu = tuple(1 if k == j else 0 for k in range(n))
                mu_j = tuple(lam_p[k] + nu[k] for k in range(n))
                if not is_dominant(mu_j):
                    continue
                if not in_Pc(mu_j, c, n):
                    continue
                V = sp.simplify(V_coeff(lam_p, nu, c, n))
                if mu_j == lam:
                    V_target = V
                else:
                    # Must be known (interior); use its M value
                    RHS_others += V * M_dict[mu_j]
            if V_target is None:
                continue
            # Solve: LHS = V_target * M_lam + RHS_others
            M_dict[lam] = sp.expand((LHS_poly - RHS_others) / V_target)
            found = True
            break
        if not found:
            print(f"  WARNING: no interior precursor found for {lam}")

    return M_dict, status_dict, P_dict, v_dict


def verify_interior_identity(n, c, max_size, verbose=True):
    """Build M's, verify interior mu satisfies M_mu = v_mu * P_mu."""
    X = make_vars(n)
    M_dict, status_dict, P_dict, v_dict = build_M_iterated(n, c, max_size, verbose=verbose)

    results = []
    for lam, M in M_dict.items():
        is_int = (status_dict[lam] == 'interior')
        expected = expand(v_dict[lam] * P_dict[lam])
        diff = sp.expand(M - expected)
        matches = sp.simplify(diff) == 0
        results.append({
            'lam': lam,
            'size': sum(lam),
            'interior': is_int,
            'matches_vP': matches,
            'M_poly': M,
        })

    return results, M_dict, P_dict, v_dict


def verify_vdez_pieri_polynomial(n, c, max_size, M_dict, P_dict, v_dict, verbose=True):
    """Check vDEZ Cor 4.2 (drop-non-dom) as POLYNOMIAL identity for each
       interior precursor lam' with |lam'| <= max_size - 1."""
    X = make_vars(n)
    m_omega1 = sum(X)
    W0_omega1 = []
    for k in range(n):
        v = [0]*n
        v[k] = 1
        W0_omega1.append(tuple(v))

    passes = []
    fails = []
    for lam in M_dict:
        if sum(lam) >= max_size:
            continue
        if not is_interior(lam, c, n):
            continue
        LHS = expand(m_omega1 * M_dict[lam])
        RHS = sp.Integer(0)
        for nu in W0_omega1:
            mu = tuple(lam[k] + nu[k] for k in range(n))
            if not is_dominant(mu):
                continue
            if not in_Pc(mu, c, n):
                continue
            if mu not in M_dict:
                continue
            V = sp.simplify(V_coeff(lam, nu, c, n))
            RHS += V * M_dict[mu]
        RHS = expand(RHS)
        diff = expand(LHS - RHS)
        ok = sp.simplify(diff) == 0
        if ok:
            passes.append(lam)
        else:
            fails.append((lam, diff))

    return passes, fails


def report_results(n, c, max_size):
    print(f"\n{'='*70}")
    print(f"n = {n}, c = {c}, max_size = {max_size}")
    print(f"{'='*70}\n")

    results, M_dict, P_dict, v_dict = verify_interior_identity(n, c, max_size, verbose=False)

    print(f"{'lam':<20}{'size':<6}{'status':<12}{'M = v*P?'}")
    print("-" * 60)
    int_pass, int_fail, bdry_pass, bdry_fail = 0, 0, 0, 0
    for r in results:
        status = 'interior' if r['interior'] else 'boundary'
        matches = 'YES' if r['matches_vP'] else 'no'
        print(f"{str(r['lam']):<20}{r['size']:<6}{status:<12}{matches}")
        if r['interior']:
            if r['matches_vP']:
                int_pass += 1
            else:
                int_fail += 1
        else:
            if r['matches_vP']:
                bdry_pass += 1
            else:
                bdry_fail += 1

    print()
    print(f"Interior:   {int_pass} match v*P, {int_fail} do NOT match v*P")
    print(f"Boundary:   {bdry_pass} match v*P, {bdry_fail} do NOT match v*P (expected non-match)")

    # Verify Pieri consistency from interior precursors
    passes, fails = verify_vdez_pieri_polynomial(n, c, max_size, M_dict, P_dict, v_dict, verbose=False)
    print(f"\nvDEZ Pieri (drop-non-dom) polynomial identity from interior precursors:")
    print(f"  {len(passes)} pass, {len(fails)} fail")
    if fails:
        print(f"  Fails:")
        for lam, diff in fails[:5]:
            print(f"    {lam}: diff = {sp.factor(diff)}")

    return results, M_dict, passes, fails


if __name__ == "__main__":
    # Test 1: n=3, c=3, |lam|<=4 (matches the 07-21 ter anchor plus extension)
    r1, M1, p1, f1 = report_results(3, 3, 4)

    # Test 2: n=3, c=4, |lam|<=5
    r2, M2, p2, f2 = report_results(3, 4, 5)

    # Test 3: n=3, c=3, |lam|<=5 (push a bit further)
    r3, M3, p3, f3 = report_results(3, 3, 5)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY: interior identity M^{(c)}_lam = v_lam * P_lam")
    print("="*70)
    for name, r in [("n=3,c=3,|lam|<=4", r1),
                    ("n=3,c=4,|lam|<=5", r2),
                    ("n=3,c=3,|lam|<=5", r3)]:
        int_pass = sum(1 for x in r if x['interior'] and x['matches_vP'])
        int_total = sum(1 for x in r if x['interior'])
        print(f"  {name}: {int_pass}/{int_total} interior lam match v_lam*P_lam")
