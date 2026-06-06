"""
Extract exact symbolic eigenvalues at n=5 for each irrep V_lambda using
Young's seminormal form restricted to each irrep (d x d matrix, manageable).

We'll compute Pi_q on each irrep V_lambda as a d x d symbolic matrix, then
factor the characteristic polynomial to extract eigenvalues as polynomials in q.

This gives EXACT eigenvalues (polynomials in q), from which we can check
divisibility by (1+q)^2 directly.
"""

import sympy as sp
from sympy import symbols, expand, factor, simplify, Poly, together, cancel
from itertools import permutations

q = symbols('q')


# --- Partitions and SYT -------------------------------------------------------

def partitions(n):
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
    cells = []
    for r, l in enumerate(lam):
        for c in range(l):
            cells.append((r, c))
    return cells


def syt_of_shape(lam):
    cells = shape_cells(lam)
    n = len(cells)
    result = []

    def place(i, filled, pos_of):
        if i == n + 1:
            result.append(dict(pos_of))
            return
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
    # Content of cell of i minus content of cell of i+1 (Ram/Hoefsmit convention
    # that makes braid relations close with b=1, c=a(d)a(-d)+q).
    return content(T[i]) - content(T[i+1])


def swap_values(T, i):
    new_T = dict(T)
    new_T[i], new_T[i+1] = T[i+1], T[i]
    r1, c1 = T[i]; r2, c2 = T[i+1]
    if r1 == r2 or c1 == c2:
        return None
    return new_T


def same_row(T, i):
    return T[i][0] == T[i+1][0]


def same_col(T, i):
    return T[i][1] == T[i+1][1]


# --- Young seminormal form using Hoefsmit/Murphy formula -------------------
# For H_q with T_i^2 = (q-1)T_i + q:
#   If i, i+1 same row: T_i v_T = q v_T.
#   If same col: T_i v_T = -v_T.
#   Else: let d = axial_dist(T, i); let
#       T_i v_T = ((q-1)/(1 - q^d)) v_T + ? v_{s_i T}
#   We'll tune the off-diagonal so braid relation holds.
#
# Standard Hoefsmit formula (Ram 1997): T_i v_T = a(d) v_T + b(d) v_{s_i T}
# where a(d) = (q-1)/(1-q^d), b(d) = sqrt(a(d) a(-d) + q) or similar.
# To avoid square roots, use Murphy's basis variant: use v_T such that
#   T_i v_T = a(d) v_T + v_{s_i T}
#   T_i v_{s_i T} = a(-d) v_{s_i T} + (1 - a(d) a(-d)) v_T  ???
# Actually the CORRECT non-symmetric form (Ram):
#   T_i v_T = a(d) v_T + (1 + a(d)) v_{s_i T}   for T "before" s_i T
#   T_i v_{s_i T} = a(-d) v_{s_i T} + (1 + a(-d)) v_T / (something) ...
# Let me use a slightly different approach: define T_i acting by
#   T_i v_T = a(d) v_T + v_{s_i T}   where a(d) = (q-1)/(1 - q^d)
#   T_i v_{s_i T} = ? v_T + ? v_{s_i T}
# must satisfy T_i^2 = (q-1) T_i + q.
# The 2x2 block in basis (v_T, v_{s_i T}) should have det = -q, tr = q-1.
# M = [[a(d), b], [c, a(-d)]] with a(d) + a(-d) = tr = q-1 and a(d)*a(-d) - bc = -q.
# a(d) + a(-d) = (q-1)/(1-q^d) + (q-1)/(1-q^{-d}) = (q-1)[(1-q^{-d}) + (1-q^d)] / [(1-q^d)(1-q^{-d})]
#              = (q-1)[2 - q^d - q^{-d}] / [2 - q^d - q^{-d}] = q-1 ✓
# a(d)*a(-d) = (q-1)^2 / [(1-q^d)(1-q^{-d})] = (q-1)^2 / (2 - q^d - q^{-d}).
# Want: bc = a(d) a(-d) + q.
# Pick b = 1 (as in Ram), then c = a(d)a(-d) + q.
# So:
#   T_i v_T = a(d) v_T + v_{s_i T}
#   T_i v_{s_i T} = (a(d) a(-d) + q) v_T + a(-d) v_{s_i T}

def build_Ti_seminormal(lam, i, syts_keyed, index):
    """Build T_i on V_lambda using Hoefsmit seminormal form (non-symmetric).

    For the pair (T, s_i T) with T.index < s_i T.index:
      T_i v_T = a(d) v_T + v_{s_i T}
      T_i v_{s_i T} = (a(d) a(-d) + q) v_T + a(-d) v_{s_i T}
    """
    d = len(syts_keyed)
    M = sp.zeros(d, d)
    seen = set()
    for idx, T in enumerate(syts_keyed):
        if same_row(T, i):
            M[idx, idx] = q
        elif same_col(T, i):
            M[idx, idx] = sp.Integer(-1)
        else:
            T2 = swap_values(T, i)
            key2 = tuple((k, T2[k]) for k in sorted(T2))
            idx2 = index[key2]
            if (idx, idx2) in seen or (idx2, idx) in seen:
                continue
            seen.add((idx, idx2))
            # canonical order: smaller idx first
            i1, i2 = min(idx, idx2), max(idx, idx2)
            # T_small -> T_large under s_i ? Either way, determine T_small
            if i1 == idx:
                T_small, T_big = T, T2
            else:
                T_small, T_big = T2, T
            d_small = axial_dist(T_small, i)   # axial dist of T_small
            # For T_big = s_i T_small, axial_dist(T_big, i) = -d_small
            a_d = (q - 1) / (1 - q**d_small)
            a_mind = (q - 1) / (1 - q**(-d_small))
            # Column i1 (v_{T_small}): T_i v_{T_small} = a(d) v_{T_small} + v_{T_big}
            M[i1, i1] = a_d
            M[i2, i1] = sp.Integer(1)
            # Column i2 (v_{T_big}): T_i v_{T_big} = (a(d) a(-d) + q) v_{T_small} + a(-d) v_{T_big}
            M[i2, i2] = a_mind
            M[i1, i2] = a_d * a_mind + q
    return sp.Matrix([[together(M[a, b]) for b in range(d)] for a in range(d)])


def verify_hecke_relations(lam, n, mats):
    """Check T_i^2 = (q-1)T_i + q and braid relations."""
    for i in range(1, n):
        Ti = mats[i]
        lhs = expand(Ti * Ti)
        rhs = expand((q-1)*Ti + q*sp.eye(Ti.shape[0]))
        diff = lhs - rhs
        for a in range(diff.shape[0]):
            for b in range(diff.shape[1]):
                if simplify(diff[a, b]) != 0:
                    return False, f"T_{i}^2 failed at ({a},{b}): {simplify(diff[a,b])}"
    for i in range(1, n-1):
        lhs = expand(mats[i] * mats[i+1] * mats[i])
        rhs = expand(mats[i+1] * mats[i] * mats[i+1])
        for a in range(lhs.shape[0]):
            for b in range(lhs.shape[1]):
                if simplify(lhs[a, b] - rhs[a, b]) != 0:
                    return False, f"braid T_{i}T_{i+1}T_{i}: ({a},{b})"
    for i in range(1, n):
        for j in range(i+2, n):
            lhs = expand(mats[i] * mats[j])
            rhs = expand(mats[j] * mats[i])
            for a in range(lhs.shape[0]):
                for b in range(lhs.shape[1]):
                    if simplify(lhs[a, b] - rhs[a, b]) != 0:
                        return False, f"commute T_{i}T_{j}"
    return True, "all good"


def staircase_word(n):
    word = []
    for k in range(1, n):
        for j in range(k, 0, -1):
            word.append(j)
    return word


def compute_Pi_on_irrep(lam, n):
    syts = syt_of_shape(lam)
    syts_keyed = sorted(syts, key=lambda T: tuple(T[k] for k in sorted(T)))
    index = {tuple((k, T[k]) for k in sorted(T)): idx
             for idx, T in enumerate(syts_keyed)}
    d = len(syts_keyed)
    if d == 0:
        return None, 0

    mats = {}
    for i in range(1, n):
        mats[i] = build_Ti_seminormal(lam, i, syts_keyed, index)

    ok, msg = verify_hecke_relations(lam, n, mats)
    if not ok:
        return None, d

    I = sp.eye(d)
    Pi = I
    for i in staircase_word(n):
        Fi = I + mats[i]
        Pi = Pi * Fi
        Pi = sp.Matrix([[together(Pi[a, b]) for b in range(d)] for a in range(d)])
    return Pi, d


def analyze_irrep(lam, n):
    print(f"\n--- lambda = {lam}, n = {n} ---")
    Pi, d = compute_Pi_on_irrep(lam, n)
    if Pi is None:
        print(f"  SKIPPED (relations failed or d=0)")
        return
    print(f"  dim V = {d}")

    # rank
    nullspace = Pi.nullspace()
    rank = d - len(nullspace)
    print(f"  rank Pi_q = {rank}  (= dim V^H)")
    if rank == 0:
        return

    # trace
    tr = sum(Pi[a, a] for a in range(d))
    tr = together(cancel(tr))
    tr = factor(tr)
    print(f"  tr Pi_q = {tr}")

    # Check divisibility of trace by (1+q)^{n//2}
    s = n // 2
    target = (1+q)**s
    try:
        tr_poly = Poly(expand(tr), q)
        div_poly = Poly(expand(target), q)
        quot, rem = sp.div(tr_poly, div_poly, q)
        if rem.as_expr() == 0:
            print(f"  (1+q)^{s} DIVIDES tr: YES, tr/(1+q)^{s} = {factor(quot.as_expr())}")
        else:
            print(f"  (1+q)^{s} DIVIDES tr: NO (remainder = {rem.as_expr()})")
    except Exception as ex:
        print(f"  div check error: {ex}")

    # Characteristic polynomial to extract eigenvalues
    if d <= 6:  # only for small dim
        x = symbols('x')
        try:
            cp = Pi.charpoly(x).as_expr()
            cp = factor(cp)
            print(f"  char poly = {cp}")
        except Exception as ex:
            print(f"  charpoly error: {ex}")


if __name__ == "__main__":
    for n in [3, 4, 5]:
        print(f"\n{'='*70}")
        print(f"  n = {n}")
        print(f"{'='*70}")
        for lam in partitions(n):
            analyze_irrep(lam, n)
