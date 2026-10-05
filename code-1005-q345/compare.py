"""The two-sided comparison.  Side A = Samuel (chains.T_tensor, root basis).
Side B = LS (ls_side.I_table, increasing chains).  Disjoint code paths."""
from collections import Counter
from itertools import permutations
from schub import length
from chains import T_tensor
from ls_side import I_table, increasing_chains
from math import factorial


def multinom(alpha):
    m = sum(alpha); r = factorial(m)
    for a in alpha: r //= factorial(a)
    return r


def run(n, verbose_examples=2):
    perms = list(permutations(range(1, n+1)))
    # tallies
    p1_bad = p1_t = 0
    p2_bad = p2_t = 0; p2_bad_flat = p2_t_flat = 0; p2_bad_big = p2_t_big = 0
    p3_bad = p3_t = 0
    strict_ex = []
    maxpart_hist = Counter(); len_hist = Counter(); ntypes = Counter()
    I_gt_N = 0

    for u in perms:
        for w in perms:
            m = length(w) - length(u)
            if m < 1:
                continue
            T = T_tensor(u, w, n)
            if not T:
                continue
            I = I_table(u, w, n)
            # --- P1: N(p) depends only on content ---
            bycontent = {}
            for p, c in T.items():
                key = tuple(sorted(p))
                bycontent.setdefault(key, set()).add(c)
            # every word of a given content must appear with the same count;
            # words absent from T have count 0, so check the full fibre
            from itertools import product as iproduct
            seen = {}
            for p, c in T.items():
                seen[p] = c
            for key in bycontent:
                alpha = [0]*(n-1)
                for k in key: alpha[k-1] += 1
                vals = set()
                for perm_word in set(permutations(key)):
                    vals.add(seen.get(perm_word, 0))
                p1_t += 1
                if len(vals) != 1:
                    p1_bad += 1
                    if p1_bad == 1:
                        print("   P1 MISMATCH", u, w, key, vals)

            # --- P2 / P3 ---
            types = set(I) | {tuple(Counter(p)[k+1] for k in range(n-1)) for p in T}
            ntypes[m] += len(types)
            for alpha in types:
                psorted = tuple(sum([[k+1]*alpha[k] for k in range(n-1)], []))
                Nval = T.get(psorted, 0)
                Ival = I.get(alpha, 0)
                mx = max(alpha) if alpha else 0
                maxpart_hist[mx] += 1
                p2_t += 1
                if mx <= 1:
                    p2_t_flat += 1
                else:
                    p2_t_big += 1
                if Ival != Nval:
                    p2_bad += 1
                    if mx <= 1: p2_bad_flat += 1
                    else: p2_bad_big += 1
                    if len(strict_ex) < verbose_examples:
                        strict_ex.append((u, w, alpha, Ival, Nval))
                if Ival > Nval:
                    I_gt_N += 1
                # P3: brief's guess
                p3_t += 1
                if Ival != multinom(alpha)*Nval:
                    p3_bad += 1
            for ch in increasing_chains(u, w, n):
                len_hist[len(ch)] += 1

    print(f"n={n}")
    print(f"  P1  N(p) depends only on content : {p1_t-p1_bad}/{p1_t} content-classes constant, {p1_bad} bad")
    print(f"  P2  I_alpha == N(p_alpha)        : {p2_t-p2_bad}/{p2_t} agree, {p2_bad} differ")
    print(f"        restricted to max part <=1 : {p2_t_flat-p2_bad_flat}/{p2_t_flat} agree, {p2_bad_flat} differ")
    print(f"        restricted to max part >=2 : {p2_t_big-p2_bad_big}/{p2_t_big} agree, {p2_bad_big} differ")
    print(f"  P2' I_alpha > N(p_alpha) ever?   : {I_gt_N} violations of I <= N")
    print(f"  P3  I_alpha == multinom*N        : {p3_t-p3_bad}/{p3_t} agree, {p3_bad} differ")
    print(f"  non-vacuity: max-part histogram over tested alphas = {dict(sorted(maxpart_hist.items()))}")
    print(f"  non-vacuity: increasing-chain LENGTH distribution  = {dict(sorted(len_hist.items()))}")
    print(f"  distinct types tested, by degree m                 = {dict(sorted(ntypes.items()))}")
    for ex in strict_ex:
        print(f"  example of I != N: u={ex[0]} w={ex[1]} alpha={ex[2]}  I={ex[3]}  N={ex[4]}")


if __name__ == "__main__":
    import sys
    for n in (3, 4):
        run(n)
