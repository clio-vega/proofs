"""
Verify the algebra identity:
    Pi_q = A * P_H * B
where for n = 4 or n = 5:
    P_H = (1+T_1)(1+T_3)        (Poincare = (1+q)^2 = (1+q)^{floor(n/2)})
    A = (1+T_1)(1+T_2)
    B = (1+T_2)(1+T_1) for n=4
    B = (1+T_2)(1+T_1)(1+T_4)(1+T_3)(1+T_2)(1+T_1) for n=5

And verify divisibility of e_lambda by (1+q)^{floor(n/2)} for each irrep with
V_lambda^H != 0.

Strategy:  use the LEFT-regular representation (N x N matrices, N = n!),
compute char poly of Pi_q, factor it, extract nonzero eigenvalues.
For n=5 (N=120), use specialization + interpolation.
"""

import sympy as sp
from sympy import symbols, expand, factor, Poly, zeros, eye, Rational
from itertools import permutations
import time

q = symbols('q')


def num_inversions(w):
    inv = 0
    for a in range(len(w)):
        for b in range(a+1, len(w)):
            if w[a] > w[b]:
                inv += 1
    return inv


def left_mult_Ti_on_Tw(i, w):
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


def build_setup(n):
    basis = sorted(permutations(range(1, n+1)),
                   key=lambda w: (num_inversions(w), w))
    perm_to_idx = {w: i for i, w in enumerate(basis)}
    return basis, perm_to_idx


def build_prod(indices, basis, perm_to_idx):
    """Build prod_{k}(1 + T_{indices[k]}) as left-mult matrix."""
    N = len(basis)
    M = sp.eye(N)
    for i in indices:
        Ti = build_Ti_left(i, basis, perm_to_idx)
        Fi = sp.eye(N) + Ti
        M = M * Fi
    for a in range(N):
        for b in range(N):
            M[a, b] = expand(M[a, b])
    return M


def staircase_word(n):
    word = []
    for k in range(1, n):
        for j in range(k, 0, -1):
            word.append(j)
    return word


def verify_factorization(n):
    print(f"\n{'='*70}")
    print(f"Verifying Pi_q = A * P_H * B for n = {n}")
    print(f"{'='*70}")

    basis, perm_to_idx = build_setup(n)
    N = len(basis)
    print(f"  N = n! = {N}")

    # Compute Pi_q directly from staircase
    print("  Building Pi_q from staircase word...")
    t0 = time.time()
    stair = staircase_word(n)
    print(f"    staircase = {stair}")
    Pi = build_prod(stair, basis, perm_to_idx)
    print(f"    Pi_q built. time = {time.time() - t0:.1f}s")

    # Compute A, P_H, B separately
    # For n=4: stair = [1,2,1,3,2,1].  A = positions 1,2 = [1,2]; P_H = pos 3,4 = [1,3]; B = pos 5,6 = [2,1].
    # For n=5: stair = [1,2,1,3,2,1,4,3,2,1]. A = [1,2]; P_H = [1,3]; B = [2,1,4,3,2,1].
    A_word = stair[:2]
    PH_word = stair[2:4]
    B_word = stair[4:]
    print(f"    A word: {A_word}")
    print(f"    P_H word: {PH_word}")
    print(f"    B word: {B_word}")

    # Confirm P_H_word = [1, 3]
    assert sorted(PH_word) == [1, 3], f"P_H word mismatch: {PH_word}"
    print(f"    P_H = (1+T_1)(1+T_3)  confirmed from staircase.")

    A = build_prod(A_word, basis, perm_to_idx)
    PH = build_prod(PH_word, basis, perm_to_idx)
    B = build_prod(B_word, basis, perm_to_idx)

    # Verify Pi_q = A * P_H * B
    APHB = A * PH * B
    for a in range(N):
        for b in range(N):
            APHB[a, b] = expand(APHB[a, b])

    diff = Pi - APHB
    max_diff = 0
    for a in range(N):
        for b in range(N):
            d = expand(diff[a, b])
            if d != 0:
                max_diff += 1
    if max_diff == 0:
        print(f"  >>> Pi_q = A * P_H * B  ALGEBRA IDENTITY CONFIRMED.")
    else:
        print(f"  !!! FAILED: {max_diff} nonzero diff entries.")
        return None

    return Pi, A, PH, B, basis, perm_to_idx


def check_divisibility_n4():
    """For n=4, factor charpoly of Pi_q symbolically and check divisibility."""
    n = 4
    Pi, A, PH, B, basis, perm_to_idx = verify_factorization(n)
    N = len(basis)

    # Characteristic polynomial of Pi_q (24x24)
    print("\n  Computing characteristic polynomial of Pi_q (may take a minute)...")
    x = symbols('x')
    t0 = time.time()
    charpoly = Pi.charpoly(x).as_expr()
    charpoly = factor(charpoly)
    print(f"    time = {time.time()-t0:.1f}s")
    print(f"    charpoly = {charpoly}")

    # Check divisibility of each nonzero factor by (1+q)^{floor(n/2)}
    target = (1 + q)**(n // 2)
    print(f"\n  Target divisor: (1+q)^{n//2} = {target}")
    # Parse factors
    if hasattr(charpoly, 'args'):
        # list of (factor, mult)
        for arg in charpoly.args:
            print(f"    factor: {arg}")


def check_divisibility_n5_numeric():
    """For n=5, specialize q at several values and compute eigenvalues."""
    n = 5
    result = verify_factorization(n)
    if result is None:
        return
    Pi, A, PH, B, basis, perm_to_idx = result
    N = len(basis)

    # Expected multiplicities:
    #   (5): mult 1 (dim=1, V^H=1),        eigenvalue e_(5)
    #   (4,1): mult 8 (dim=4, V^H=2),      e_(4,1)
    #   (3,2): mult 10 (dim=5, V^H=2),     e_(3,2)
    #   (3,1,1): mult 6 (dim=6, V^H=1),    e_(3,1,1)
    #   (2,2,1): mult 5 (dim=5, V^H=1),    e_(2,2,1)
    # Plus 90 zero eigenvalues.
    expected_mults = [1, 8, 10, 6, 5]  # sum = 30, plus 90 zeros = 120
    print(f"\n  Expected multiplicity pattern: {expected_mults} + 90 zeros = 120")

    # Specialize q at integer values, compute eigenvalues
    q_values = [2, 3, 5, 7, 11]
    eigvals_data = {}  # q -> sorted list of nonzero eigenvalues (by magnitude, then as tuples)

    import numpy as np
    for q_val in q_values:
        Pi_num = Pi.subs(q, q_val)
        # Convert to numpy
        Pi_np = np.array(Pi_num, dtype=object)
        # Get int matrix
        Pi_int = np.zeros((N, N), dtype=int)
        for a in range(N):
            for b in range(N):
                Pi_int[a, b] = int(Pi_num[a, b])
        # Eigenvalues via numpy (floats, just for inspection)
        eigs = np.linalg.eigvals(Pi_int.astype(float))
        eigs_nonzero = [e for e in eigs if abs(e) > 1e-6]
        eigs_sorted = sorted([(round(e.real), round(e.imag)) for e in eigs_nonzero])
        from collections import Counter
        counter = Counter(eigs_sorted)
        print(f"\n  q = {q_val}:")
        print(f"    rank of Pi_q at q={q_val}: {sum(1 for e in eigs if abs(e) > 1e-6)}")
        for val, ct in sorted(counter.items(), key=lambda x: -x[1]):
            print(f"      eig = {val}, mult = {ct}")
        eigvals_data[q_val] = counter


if __name__ == "__main__":
    print("*" * 70)
    print("Step 1: verify Pi_q = A * P_H * B for n=4")
    print("*" * 70)
    check_divisibility_n4()
    print()
    print("*" * 70)
    print("Step 2: verify Pi_q = A * P_H * B for n=5, compute eigenvalues numerically")
    print("*" * 70)
    check_divisibility_n5_numeric()
