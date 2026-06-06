"""
Young's q-seminormal form for H_q(S_n) with T_i^2 = (q-1)T_i + q.

Following Ram "Seminormal representations of Hecke and Iwahori-Hecke algebras"
(J. Alg 1997).

For lambda |- n, basis {v_T : T in SYT(lambda)} of V^lambda.  For each SYT T and
each simple reflection s_i (1 <= i <= n-1):

  Let d = c(T^{-1}(i+1)) - c(T^{-1}(i)), content difference.

  (i)   If i, i+1 same row (d = 1):  T_i v_T = q v_T.
  (ii)  If i, i+1 same column (d = -1):  T_i v_T = -v_T.
  (iii) Otherwise:  T_i v_T = alpha(d) v_T + gamma(d) v_{s_i T},
        where alpha(d) = (q-1)/(1-q^d),  gamma(d) = 1 + alpha(d) = (q-q^d)/(1-q^d).

Verified:  trace = q-1, det = -q for 2x2 block, so eigenvalues q, -1.

Usage:  build matrices T_i for each irrep, verify braid relations, then use.
"""

import sympy as sp
from sympy import Rational, simplify, expand, factor, symbols, together
from itertools import permutations

q = symbols('q')


# --------- Partitions and SYT -------------------------------------------------

def partitions(n):
    """All partitions of n, as decreasing tuples."""
    if n == 0:
        return [()]
    result = []

    def gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(current))
            return
        for k in range(min(remaining, max_part), 0, -1):
            gen(remaining - k, k, current + [k])

    gen(n, n, [])
    return result


def shape_cells(lam):
    """Cells of shape lam: list of (row, col), 0-indexed, row-major."""
    cells = []
    for r, l in enumerate(lam):
        for c in range(l):
            cells.append((r, c))
    return cells


def syt_of_shape(lam):
    """Return list of SYT of shape lam, each as a dict (r,c) -> entry in {1,...,n}
    (equivalently, as a tuple of tuples representing the filling)."""
    cells = shape_cells(lam)
    n = len(cells)
    result = []

    # Build SYT by placing 1, 2, ..., n; at each step only cells where row>=1 entry
    # and col>=1 entry are already placed.
    def place(i, filled, pos_of):
        if i == n + 1:
            # Store as tuple of (entry-ordered) positions:  pos_of[k] = cell of k
            result.append(dict(pos_of))
            return
        # i can go in any cell (r,c) currently "addable": cell is empty, and (r-1,c) filled or r=0, and (r,c-1) filled or c=0
        for (r, c) in cells:
            if (r, c) in filled:
                continue
            if r > 0 and (r-1, c) not in filled:
                continue
            if c > 0 and (r, c-1) not in filled:
                continue
            filled.add((r, c))
            pos_of[i] = (r, c)
            place(i + 1, filled, pos_of)
            filled.remove((r, c))
            del pos_of[i]

    place(1, set(), {})
    return result


def content(cell):
    r, c = cell
    return c - r


def axial_dist(T, i):
    return content(T[i+1]) - content(T[i])


def swap_values(T, i):
    """Return s_i T: swap values i and i+1.  Return None if not standard."""
    new_T = dict(T)
    new_T[i], new_T[i+1] = T[i+1], T[i]
    # Check standard:
    # After swapping entry i and i+1: we need the new filling to still be standard
    # i.e., rows increasing left-to-right, cols increasing top-to-bottom.
    # Easiest: just check that (r_i, c_i) and (r_{i+1}, c_{i+1}) are not adjacent in same row/col.
    r1, c1 = T[i]
    r2, c2 = T[i+1]
    if r1 == r2:  # same row (must be adjacent then, since filling is SYT)
        return None
    if c1 == c2:  # same column
        return None
    return new_T


def same_row(T, i):
    return T[i][0] == T[i+1][0]


def same_col(T, i):
    return T[i][1] == T[i+1][1]


# --------- Build T_i matrices for irrep V^lambda -----------------------------

def build_Ti_on_Vlambda(lam, i):
    """Return T_i as a matrix on V^lambda (as sympy Matrix over Q(q))."""
    syts = syt_of_shape(lam)
    d = len(syts)
    # Order syts canonically (by some key).
    syts_keyed = sorted(syts, key=lambda T: tuple(T[k] for k in sorted(T)))
    index = {}
    for idx, T in enumerate(syts_keyed):
        # Hashable key: sorted tuple of (entry, pos)
        key = tuple((k, T[k]) for k in sorted(T))
        index[key] = idx

    def key_of(T):
        return tuple((k, T[k]) for k in sorted(T))

    M = sp.zeros(d, d)
    for idx, T in enumerate(syts_keyed):
        if same_row(T, i):
            M[idx, idx] = q
        elif same_col(T, i):
            M[idx, idx] = sp.Integer(-1)
        else:
            T2 = swap_values(T, i)
            idx2 = index[key_of(T2)]
            d_T = axial_dist(T, i)
            alpha = (q - 1) / (1 - q**d_T)
            gamma = 1 + alpha  # = (q - q^d)/(1 - q^d)
            # T_i v_T = alpha v_T + gamma v_{s_i T}
            M[idx, idx] += alpha
            M[idx2, idx] += gamma
            # Note: this fills COLUMN idx.  When we later fill column idx2 (coming
            # from T2), it will fill in T_i v_{T2}:  alpha' v_{T2} + gamma' v_T.
            # So both columns get filled correctly when we iterate over all T.
    # Simplify
    for a in range(d):
        for b in range(d):
            M[a, b] = together(M[a, b])
    return M, syts_keyed, index


def build_Ti_matrices(lam, n):
    """Build T_1, ..., T_{n-1} on V^lambda."""
    mats = {}
    syts = None
    index = None
    for i in range(1, n):
        M, syts, index = build_Ti_on_Vlambda(lam, i)
        mats[i] = M
    return mats, syts, index


def verify_relations(mats, n):
    """Verify T_i^2 = (q-1) T_i + q, braid T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1},
    commutativity for |i-j| >= 2."""
    ok = True
    for i in range(1, n):
        Ti = mats[i]
        lhs = expand(Ti * Ti)
        rhs = expand((q-1) * Ti + q * sp.eye(Ti.shape[0]))
        diff = lhs - rhs
        diff.simplify()
        if not all(simplify(diff[a, b]) == 0 for a in range(diff.shape[0]) for b in range(diff.shape[1])):
            print(f"  FAIL: T_{i}^2 != (q-1) T_{i} + q")
            ok = False
    for i in range(1, n-1):
        Ti = mats[i]
        Tip1 = mats[i+1]
        lhs = expand(Ti * Tip1 * Ti)
        rhs = expand(Tip1 * Ti * Tip1)
        diff = lhs - rhs
        if not all(simplify(diff[a, b]) == 0 for a in range(diff.shape[0]) for b in range(diff.shape[1])):
            print(f"  FAIL: braid T_{i} T_{i+1} T_{i} != T_{i+1} T_{i} T_{i+1}")
            ok = False
    for i in range(1, n):
        for j in range(i+2, n):
            Ti = mats[i]
            Tj = mats[j]
            lhs = expand(Ti * Tj)
            rhs = expand(Tj * Ti)
            diff = lhs - rhs
            if not all(simplify(diff[a, b]) == 0 for a in range(diff.shape[0]) for b in range(diff.shape[1])):
                print(f"  FAIL: commutativity T_{i} T_{j} != T_{j} T_{i}")
                ok = False
    return ok


# --------- Build Pi_q on V^lambda --------------------------------------------

def staircase_word(n):
    word = []
    for k in range(1, n):
        for j in range(k, 0, -1):
            word.append(j)
    return word


def build_Pi_on_Vlambda(lam, n):
    mats, syts, index = build_Ti_matrices(lam, n)
    d = len(syts)
    I = sp.eye(d)
    Pi = I
    for i in staircase_word(n):
        Fi = I + mats[i]
        Pi = Pi * Fi
    # Simplify entries
    for a in range(d):
        for b in range(d):
            Pi[a, b] = together(expand(Pi[a, b]))
    return Pi, mats, syts, index


def build_PH_on_Vlambda(lam, n, odd_indices):
    mats, syts, index = build_Ti_matrices(lam, n)
    d = len(syts)
    I = sp.eye(d)
    PH = I
    for i in odd_indices:
        Fi = I + mats[i]
        PH = PH * Fi
    for a in range(d):
        for b in range(d):
            PH[a, b] = together(expand(PH[a, b]))
    return PH, mats, syts, index


# --------- Main: test for n=4, n=5 --------------------------------------------

def test_lambda(lam, n):
    print(f"\n--- lambda = {lam}, n = {n} ---")
    mats, syts, index = build_Ti_matrices(lam, n)
    d = len(syts)
    print(f"  dim V^lambda = {d}")
    if d == 0:
        return
    ok = verify_relations(mats, n)
    print(f"  Hecke relations: {'OK' if ok else 'FAIL'}")
    if not ok:
        return

    # Build Pi_q
    I = sp.eye(d)
    Pi = I
    for i in staircase_word(n):
        Fi = I + mats[i]
        Pi = Pi * Fi
    Pi = sp.simplify(Pi)
    # Trace
    tr_pi = sp.simplify(sum(Pi[a, a] for a in range(d)))
    tr_pi = factor(tr_pi)
    print(f"  tr Pi_q  = {tr_pi}")

    # Rank via nullspace
    try:
        null = Pi.nullspace()
        rank = d - len(null)
    except Exception as e:
        print(f"  rank: error {e}")
        rank = None
    print(f"  rank Pi_q = {rank}  (V^H dim expected)")

    # If rank > 0, extract nonzero eigenvalue via trace / rank
    if rank and rank > 0:
        e_lam = sp.simplify(tr_pi / rank)
        e_lam = factor(e_lam)
        print(f"  e_lambda = tr/rank = {e_lam}")

        # Check divisibility by (1+q)^{floor(n/2)}
        target = (1 + q)**(n // 2)
        quot = sp.cancel(e_lam / target)
        quot_simplified = sp.factor(quot)
        # Check if quotient is a polynomial (in Q[q])
        # i.e., no 1/(1+q) factors
        poly_check = sp.Poly(sp.expand(e_lam), q)
        # Check (1+q)^{floor(n/2)} divides e_lambda
        try:
            quot_poly = sp.div(poly_check, sp.Poly(sp.expand(target), q), q)
            # sp.div returns (quotient, remainder)
            rem = quot_poly[1]
            q_poly = quot_poly[0]
            is_div = (rem.as_expr() == 0)
        except Exception as ex:
            print(f"  Divisibility check error: {ex}")
            is_div = None
        if is_div is True:
            print(f"  (1+q)^{n//2} DIVIDES e_lambda: YES")
            print(f"    e_lambda / (1+q)^{n//2} = {factor(q_poly.as_expr())}")
        elif is_div is False:
            print(f"  (1+q)^{n//2} DIVIDES e_lambda: NO  (remainder = {rem.as_expr()})")
        else:
            print(f"  (1+q)^{n//2} DIVIDES: unknown")

    return Pi, mats, syts, index


if __name__ == "__main__":
    for n in [3, 4]:
        print(f"\n{'='*70}\n  n = {n}\n{'='*70}")
        for lam in partitions(n):
            test_lambda(lam, n)
