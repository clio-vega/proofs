"""
Corrected version: test divisibility of P_H T_w P_H by (1+q)^s for n=4.

Also test the CORRECT algebraic statement that might hold:
  P_H * H_q * P_H lies in (1+q)^s * H_q  ... which we now expect to be FALSE.

If it's false, we need a refined statement.
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


def reduced_word(w):
    """Compute a reduced word for w: list of simple reflection indices
    s.t. w = s_{r_1} s_{r_2} ... s_{r_k} (in Coxeter notation)."""
    n = len(w)
    cur = list(w)
    reduced = []
    while cur != sorted(cur):
        found = False
        for j in range(1, n):
            pos_j = cur.index(j)
            pos_jp1 = cur.index(j+1)
            if pos_jp1 < pos_j:  # val j+1 before val j; apply s_j to reduce
                cur[pos_j] = j+1
                cur[pos_jp1] = j
                reduced.append(j)
                found = True
                break
        if not found:
            raise RuntimeError(f"stuck at {cur}")
    # reduced applied to w reaches id; so w = s_{reduced[0]} s_{reduced[1]} ... s_{reduced[-1]}
    # Actually: s_{r_k} ... s_{r_1} w = id => w = s_{r_1} s_{r_2} ... s_{r_k}
    # So the reduced word for w is [r_1, r_2, ..., r_k].
    # But reduced list is built in order r_1, r_2, ..., r_k. Wait:
    # iter 1: apply s_{r_1} on left of cur=w to get s_{r_1} w. cur becomes s_{r_1} w.
    # iter 2: apply s_{r_2} on left. cur becomes s_{r_2} s_{r_1} w.
    # ... iter k: cur = s_{r_k} ... s_{r_1} w = id.
    # So w = (s_{r_k} ... s_{r_1})^{-1} = s_{r_1} ... s_{r_k}.
    # Hence reduced word for w (reading left to right) = [r_1, r_2, ..., r_k].
    # This is the standard form.
    return reduced


def build_Tw_from_reduced(reduced, basis, perm_to_idx):
    """Build T_w = T_{r_1} T_{r_2} ... T_{r_k} as left-mult matrix."""
    N = len(basis)
    M = sp.eye(N)
    for r in reduced:
        M = M * build_Ti_left(r, basis, perm_to_idx)
    for a in range(N):
        for b in range(N):
            M[a, b] = expand(M[a, b])
    return M


def div_info(poly_expr, s):
    """Return max power of (1+q) dividing poly_expr."""
    expr = expand(poly_expr)
    if expr == 0:
        return float('inf')
    for k in range(s, -1, -1):
        try:
            ent_poly = Poly(expr, q)
            div_poly = Poly(expand((1+q)**k), q)
            quot, rem = sp.div(ent_poly, div_poly, q)
            if rem.as_expr() == 0:
                return k
        except Exception:
            pass
    return 0


def min_divisibility(M, s):
    """Return the MINIMUM over entries of (max k: (1+q)^k divides entry)."""
    mn = float('inf')
    for a in range(M.shape[0]):
        for b in range(M.shape[1]):
            entry = expand(M[a, b])
            if entry == 0:
                continue
            k = div_info(entry, s)
            if k < mn:
                mn = k
    return mn


def test_n(n):
    print(f"\n{'='*70}")
    print(f"Testing divisibility of P_H T_w P_H by (1+q)^s (s={n//2}) for n={n}")
    print(f"{'='*70}")

    basis = sorted(permutations(range(1, n+1)),
                   key=lambda w: (num_inversions(w), w))
    perm_to_idx = {w: i for i, w in enumerate(basis)}
    N = len(basis)
    s = n // 2

    odd_gens = list(range(1, n, 2))
    print(f"  Odd generators (H): {odd_gens}")

    PH = sp.eye(N)
    for i in odd_gens:
        PH = PH * (sp.eye(N) + build_Ti_left(i, basis, perm_to_idx))
    for a in range(N):
        for b in range(N):
            PH[a, b] = expand(PH[a, b])

    # Check each w
    min_divs = {}
    for w in basis:
        red = reduced_word(w)
        Tw = build_Tw_from_reduced(red, basis, perm_to_idx)
        Prod = PH * Tw * PH
        for a in range(N):
            for b in range(N):
                Prod[a, b] = expand(Prod[a, b])
        md = min_divisibility(Prod, s)
        min_divs[w] = md
        status = "DIV" if md >= s else f"not div ((1+q)^{md} only)"
        print(f"    w={w} (red={red}):  min div = (1+q)^{md}  -- {status}")

    # Summary
    overall = min(min_divs.values())
    print(f"\n  Overall: min_{{w, entries}} div = (1+q)^{overall}")
    print(f"  So (1+q)^{overall} | P_H T_w P_H for all w, but NOT (1+q)^{overall+1}.")

    # Print double-coset representatives
    # H = <s_1, s_3>. Double cosets H\S_n/H represented by minimal-length elts.
    print(f"\n  Double-coset minimal reps and their divisibility:")
    # Find minimal-length rep per double coset.
    # For each w, coset HwH = {h1 * w * h2: h1, h2 in H}.
    if n == 4:
        H_elems = [(1,2,3,4), (2,1,3,4), (1,2,4,3), (2,1,4,3)]  # H = <s_1, s_3>
    elif n == 5:
        H_elems = [(1,2,3,4,5), (2,1,3,4,5), (1,2,4,3,5), (2,1,4,3,5),
                   (1,2,3,5,4)]  # <s_1, s_3>... wait n=5 odd gens are 1, 3; <s_1, s_3> in S_5
        # Actually H = <s_1, s_3> in S_5 = S_2 x S_2 x S_1 on {1,2},{3,4},{5}, size 4
        H_elems = [(1,2,3,4,5), (2,1,3,4,5), (1,2,4,3,5), (2,1,4,3,5)]

    def mult_perm(a, b):
        """Return a*b as perm: (a*b)(i) = a(b(i))."""
        return tuple(a[b[i]-1] for i in range(len(a)))

    coset_reps = {}
    for w in basis:
        coset = set()
        for h1 in H_elems:
            for h2 in H_elems:
                coset.add(mult_perm(mult_perm(h1, w), h2))
        # Find min-length rep
        min_rep = min(coset, key=lambda u: (num_inversions(u), u))
        coset_reps.setdefault(min_rep, []).append(w)

    for rep, members in sorted(coset_reps.items(), key=lambda x: num_inversions(x[0])):
        print(f"    double-coset rep {rep} (length {num_inversions(rep)}): "
              f"{len(members)} elements; div = (1+q)^{min_divs[rep]}")


if __name__ == "__main__":
    test_n(4)
