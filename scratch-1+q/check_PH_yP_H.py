"""
Test whether P_H * T_w * P_H is divisible by (1+q)^{floor(n/2)} in H_q(S_n)
for all w in S_n.  If yes, this closes the proof.

For n=4, s = 2.  For n=5, s = 2.  (Same s.)

We test: in the regular rep, the matrix of P_H T_w P_H (as left-mult operator)
has every entry divisible by (1+q)^s.
"""

import sympy as sp
from sympy import symbols, expand, Poly, zeros, eye
from itertools import permutations

q = symbols('q')


def num_inversions(w):
    return sum(1 for a in range(len(w)) for b in range(a+1, len(w)) if w[a] > w[b])


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


def build_Tw_left(w_reduced, basis, perm_to_idx):
    """Build T_w as left-mult matrix, given w as a list of simple reflection indices."""
    N = len(basis)
    M = sp.eye(N)
    for i in w_reduced:
        M = M * build_Ti_left(i, basis, perm_to_idx)
    for a in range(N):
        for b in range(N):
            M[a, b] = expand(M[a, b])
    return M


def check_divisibility_by_1pq_s(M, s):
    """Check each entry of M is divisible by (1+q)^s."""
    divisor = (1 + q)**s
    bad_entries = 0
    for a in range(M.shape[0]):
        for b in range(M.shape[1]):
            entry = expand(M[a, b])
            if entry == 0:
                continue
            # Polynomial division
            try:
                ent_poly = Poly(entry, q)
                div_poly = Poly(expand(divisor), q)
                quot, rem = sp.div(ent_poly, div_poly, q)
                if rem.as_expr() != 0:
                    bad_entries += 1
            except Exception as ex:
                bad_entries += 1
    return bad_entries


def test_n(n):
    print(f"\n{'='*70}")
    print(f"Testing P_H y P_H divisibility by (1+q)^{n//2} for n = {n}")
    print(f"{'='*70}")

    basis = sorted(permutations(range(1, n+1)),
                   key=lambda w: (num_inversions(w), w))
    perm_to_idx = {w: i for i, w in enumerate(basis)}
    N = len(basis)
    s = n // 2

    # Build P_H = prod_{i odd}(1 + T_i)
    odd_gens = list(range(1, n, 2))
    print(f"  H-generators (odd): {odd_gens}")
    PH = sp.eye(N)
    for i in odd_gens:
        PH = PH * (sp.eye(N) + build_Ti_left(i, basis, perm_to_idx))
    for a in range(N):
        for b in range(N):
            PH[a, b] = expand(PH[a, b])

    # Test each T_w
    all_divisible = True
    checked = 0
    for w in basis:
        # Build T_w from reduced word (compute inversions count to get reduced word)
        # Easiest: use the Coxeter descent procedure
        reduced = []
        w_cur = list(w)
        while True:
            # find leftmost descent
            found = False
            for i in range(len(w_cur) - 1):
                if w_cur[i] > w_cur[i+1]:
                    # this is an inversion: apply s_i on the left (swap values)
                    # we want to build T_w bottom-up, so collect reflections
                    # Actually let's just use: w = s_i * w' where ell(w') < ell(w)
                    # If w has inversion at (i, i+1) in one-line notation (w_cur[i] > w_cur[i+1]),
                    # then w' = s_{w_cur[i+1]} w swaps values to reduce...
                    # Easier: build w directly from sorted.
                    pass
            break

        # Simpler: w = pi_{w_inv_count} sorting process. Use sympy permutation?
        # Alternative: use position-based left multiplication to compute T_w.
        # Since basis sorts by inversions, and left_mult_Ti swaps values i,i+1,
        # we can compute T_w by repeatedly applying s_i to identity.
        # But actually the left-mult matrix M_{T_w} where T_w is defined by the
        # canonical reduced expression...
        # Let's just compute M_{T_w} as: start with T_e = I, then apply left_mult_Ti
        # to build T_w via its reduced expression.
        # Need reduced expression for w.
        reduced = []
        cur = list(w)
        # bubble sort: find inversion, apply s to remove it; record s_i used
        while True:
            found = False
            for i in range(len(cur) - 1):
                # if cur has values a, b at positions i, i+1 with a > b, applying
                # left s_j where j = b (value) swaps values j and j+1 = b and b+1 = a?
                # Wait, we need a = b+1 for simple reflection.
                if cur[i] > cur[i+1] and cur[i] == cur[i+1] + 1:
                    val = cur[i+1]
                    # apply s_{val} on left: swap values val and val+1
                    cur[i], cur[i+1] = cur[i+1], cur[i]
                    reduced.append(val)
                    found = True
                    break
            if not found:
                if cur == sorted(cur):
                    break
                # no adjacent simple inversion; do any swap that reduces by other means
                # Actually for Sn this should always work for adjacent inversions
                # But could have cur = [3, 1, 2]: swap positions 0,1 -> val 1, get [1,3,2]; swap pos 1,2 val 2 -> [1,2,3]
                # so it does always work with adjacent inversions of consecutive values
                # If stuck, print:
                print(f"    stuck at {cur}")
                break
        # reduced gives reduced word (as list of simple reflection indices) to get w
        # from sorted (identity). So T_w = T_{s_{reduced[k]}} ... T_{s_{reduced[0]}}?
        # Hmm, check direction.
        # left_mult_Ti_on_Tw(i, w) computes T_i * T_w. Starting from T_e=T_{(1,2,..,n)}
        # and applying left mults: T_{s_{i_k}} ... T_{s_{i_1}} * T_e = T_{s_{i_k} ... s_{i_1}}.
        # So if we want T_w with w = s_{i_1} s_{i_2} ... s_{i_k} (reading left-to-right
        # in reduced word), we apply mults in reverse order i_k, i_{k-1}, ..., i_1.
        # But reduced built above gives the reflections applied to cur (starting from w)
        # to reach identity. If reduced = [r_1, r_2, ..., r_k], that means
        # s_{r_k} ... s_{r_2} s_{r_1} * w = identity, so w = s_{r_1} s_{r_2} ... s_{r_k}.
        # To build T_w, apply left mult by s_{r_1}, then by s_{r_2}, ...
        # Wait actually that's NOT standard. Standard: T_{s_{r_1}} * T_w' where w' = s_{r_1} w,
        # and ell(w') < ell(w). So T_w can be computed as:
        #   T_w = T_{r_1} T_{r_2} ... T_{r_k} (product of generators in reduced order)
        # Hmm, actually if w = s_{r_1} s_{r_2} ... s_{r_k}, then T_w = T_{r_1} T_{r_2} ... T_{r_k}.
        # And as left-mult matrix on T-basis: (M_{T_{r_1}}) (M_{T_{r_2}}) ... (M_{T_{r_k}}).
        # My reduced word is such that s_{r_k} ... s_{r_1} * w = id, so w = s_{r_1} s_{r_2} ... s_{r_k}.

        # Wait actually let me re-examine. Initially cur = w. We apply s_val (by swapping
        # values val, val+1 at positions). That corresponds to left-multiplying w by
        # s_val: cur_new = s_val * cur = s_val * w, which has smaller length.
        # So after k steps, cur = s_{r_k} s_{r_{k-1}} ... s_{r_1} w = id.
        # Hence w = s_{r_1}^{-1} ... s_{r_k}^{-1} = s_{r_1} ... s_{r_k} (since s^{-1} = s).
        # So w = s_{r_1} s_{r_2} ... s_{r_k} (reduced word).
        # T_w = T_{r_1} T_{r_2} ... T_{r_k}.
        # As left-mult matrix: product of M_{r_1}, M_{r_2}, ..., M_{r_k} in that order.

        M_Tw = sp.eye(N)
        for r in reduced:
            M_Tw = M_Tw * build_Ti_left(r, basis, perm_to_idx)
        for a in range(N):
            for b in range(N):
                M_Tw[a, b] = expand(M_Tw[a, b])

        # Compute P_H * T_w * P_H
        Prod = PH * M_Tw * PH
        for a in range(N):
            for b in range(N):
                Prod[a, b] = expand(Prod[a, b])

        bad = check_divisibility_by_1pq_s(Prod, s)
        checked += 1
        status = "OK" if bad == 0 else f"FAIL ({bad} entries not divisible)"
        print(f"    w = {w}: P_H * T_w * P_H  divisibility by (1+q)^{s}: {status}")
        if bad > 0:
            all_divisible = False

    print(f"\n  Total checked: {checked} elements.")
    print(f"  All divisible by (1+q)^{s}: {'YES' if all_divisible else 'NO'}")


if __name__ == "__main__":
    test_n(4)
