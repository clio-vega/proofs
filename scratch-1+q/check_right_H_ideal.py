"""
Test the key structural claim:

  Pi_q \in A * P_H  in H_q(S_n), where A = H_q(S_n), P_H = prod over i odd of (1+T_i).

Equivalently:  Pi_q * T_i = q * Pi_q  for every odd i in {1, 3, 5, ..., n-1 or n-2}.

If true for all such i, then Pi_q = X * P_H for some X, and the theorem
  (1+q)^{floor(n/2)} divides e_lambda(q)
follows immediately.

We verify for n = 3, 4, 5 by building Pi_q in the T-basis regular representation
and checking the relation.

Conventions match build_factors.py:
  T_i^2 = (q-1) T_i + q.
  Left multiplication on T-basis:
    T_i T_w = T_{s_i w}          if  ell(s_i w) > ell(w)
            = q T_{s_i w} + (q-1) T_w    otherwise.
  (where s_i acts on one-line notation by swapping VALUES i and i+1).
"""

import sympy as sp
from sympy import symbols, expand, zeros, Poly
from itertools import permutations

q = symbols('q')


def num_inversions(w):
    inv = 0
    for a in range(len(w)):
        for b in range(a+1, len(w)):
            if w[a] > w[b]:
                inv += 1
    return inv


def left_mult_Ti_on_Tw(i, w):
    """T_i * T_w, returning [(coeff, new_w), ...]."""
    pos_i = w.index(i)
    pos_ip1 = w.index(i + 1)
    new_w = list(w)
    new_w[pos_i] = i + 1
    new_w[pos_ip1] = i
    new_w = tuple(new_w)
    if pos_i < pos_ip1:
        return [(sp.Integer(1), new_w)]
    else:
        return [(q, new_w), (q - 1, w)]


def right_mult_Tw_by_Ti(w, i):
    """T_w * T_i, returning [(coeff, new_w), ...].
    Right mult by T_i:
        T_w * T_i = T_{w s_i}     if ell(w s_i) > ell(w)
                  = q T_{w s_i} + (q-1) T_w  otherwise.
    Here s_i is the transposition swapping POSITIONS i and i+1 in the one-line
    notation (since w s_i means apply s_i first, then w, so in one-line we swap
    positions).
    """
    new_w = list(w)
    new_w[i-1], new_w[i] = new_w[i], new_w[i-1]
    new_w = tuple(new_w)
    # ell(w s_i) > ell(w) iff w(i) < w(i+1) in one-line (ascending at position i)
    if w[i-1] < w[i]:
        return [(sp.Integer(1), new_w)]
    else:
        return [(q, new_w), (q - 1, w)]


def build_Ti_left(i, basis, perm_to_idx):
    N = len(basis)
    M = sp.zeros(N, N)
    for j, v in enumerate(basis):
        for coeff, u in left_mult_Ti_on_Tw(i, v):
            M[perm_to_idx[u], j] += coeff
    for a in range(N):
        for b in range(N):
            M[a, b] = expand(M[a, b])
    return M


def build_Fi_left(i, basis, perm_to_idx):
    """(1 + T_i) as left-mult matrix."""
    M = build_Ti_left(i, basis, perm_to_idx)
    return sp.eye(len(basis)) + M


def build_Ti_right(i, basis, perm_to_idx):
    """T_i acting by RIGHT multiplication.  If M = M_{T_i}^right, then
    (a * T_i) as algebra element has LEFT-mult matrix M_a * M_{T_i}^left.
    But here we want:  if element x has LEFT-mult matrix M_x, then the element
    x * T_i has LEFT-mult matrix M_x * M_{T_i}^left.  So we just use M_{T_i}^left.
    No separate construction needed.
    """
    return build_Ti_left(i, basis, perm_to_idx)


def build_Pi(staircase, basis, perm_to_idx):
    """Build Pi_q = prod_k (1 + T_{staircase[k]})  as left-mult matrix."""
    N = len(basis)
    M = sp.eye(N)
    for i in staircase:
        M = M * build_Fi_left(i, basis, perm_to_idx)
    for a in range(N):
        for b in range(N):
            M[a, b] = expand(M[a, b])
    return M


def staircase_word(n):
    """Return the staircase reduced word for w_0 in S_n:
       s_1 (s_2 s_1) (s_3 s_2 s_1) ... (s_{n-1} ... s_1)."""
    word = []
    for k in range(1, n):
        for j in range(k, 0, -1):
            word.append(j)
    return word


def test_n(n):
    print(f"\n{'=' * 70}")
    print(f"Testing n = {n}")
    print(f"{'=' * 70}")

    N = 1
    for k in range(1, n+1):
        N *= k
    basis = sorted(permutations(range(1, n+1)),
                   key=lambda w: (num_inversions(w), w))
    perm_to_idx = {w: i for i, w in enumerate(basis)}

    stair = staircase_word(n)
    print(f"Staircase word: {stair}  (length {len(stair)} = C({n},2) = {n*(n-1)//2})")

    print("Building Pi_q (this may take a while for n=5)...")
    Pi = build_Pi(stair, basis, perm_to_idx)
    print(f"Pi_q built as {N}x{N} matrix over Z[q].")

    # Odd indices s_1, s_3, ... generating H
    H_odd = list(range(1, n, 2))
    print(f"H = <{['s_' + str(i) for i in H_odd]}> = (S_2)^{len(H_odd)}")
    print(f"Poincare polynomial of H: (1+q)^{len(H_odd)}")
    print(f"floor(n/2) = {n // 2}, matches len(H_odd) = {len(H_odd)}: "
          f"{'YES' if n // 2 == len(H_odd) else 'NO'}")

    for i in H_odd:
        print(f"\n  Checking Pi_q * T_{i} = q * Pi_q ...")
        Ti = build_Ti_left(i, basis, perm_to_idx)
        LHS = Pi * Ti
        RHS = q * Pi
        # Check equality
        diff = LHS - RHS
        # Reduce each entry mod q-expansion
        nonzero = 0
        for a in range(N):
            for b in range(N):
                d = expand(diff[a, b])
                if d != 0:
                    nonzero += 1
        if nonzero == 0:
            print(f"    >>> HOLDS: Pi_q * T_{i} = q * Pi_q  (all entries match).")
        else:
            print(f"    !!! FAILS: {nonzero} nonzero entries in diff.")
            # Print first few
            count = 0
            for a in range(N):
                for b in range(N):
                    d = expand(diff[a, b])
                    if d != 0 and count < 5:
                        print(f"        diff[{a},{b}] = {d}")
                        count += 1

    # Also: does Pi_q * T_even = ??
    print(f"\n  Spot-check: Pi_q * T_2 (even index) -- expected NOT equal q*Pi_q")
    if n >= 3:
        T2 = build_Ti_left(2, basis, perm_to_idx)
        LHS = Pi * T2
        RHS = q * Pi
        diff = LHS - RHS
        nonzero = 0
        for a in range(N):
            for b in range(N):
                if expand(diff[a, b]) != 0:
                    nonzero += 1
        print(f"    nonzero entries in Pi_q*T_2 - q*Pi_q: {nonzero}  (should be > 0 for sanity)")

    return Pi, basis, perm_to_idx


if __name__ == "__main__":
    # Start with small n
    for n in [3, 4]:
        test_n(n)
