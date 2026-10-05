"""The load-bearing identity of Theorem 4, verified against a DIFFERENT mechanism.

   S_u * h_alpha  =  sum_w I_alpha(u,w) S_w          (LS side)
   S_u * Y^beta   =  sum_w N_{w/u}(p_beta) S_w       (Samuel side, my (**) of 10-04)

LHS computed with the divided-difference Schubert table (no chains anywhere).
RHS computed with the increasing-chain enumerator / the root-basis tensor.
This makes Theorem 4 independent of Sottile's [So96], which I have NOT read.
"""
from collections import Counter
from itertools import permutations, combinations_with_replacement
from schub import length, schubert_table, pmul, extract_coeff
from chains import T_tensor
from ls_side import I_table
from transition import comps_le_delta, word

def hpoly(a, k, N):
    """h_a(x_1..x_k) as a dict exponent->coeff in N variables."""
    out = Counter()
    for c in combinations_with_replacement(range(k), a):
        e = [0]*N
        for i in c: e[i] += 1
        out[tuple(e)] += 1
    return dict(out)

def run(n):
    N = n+1
    tab = schubert_table(n, N)
    perms = list(permutations(range(1, n+1)))
    byl = {}
    for v in perms: byl.setdefault(length(v), []).append(v)
    maxm = n*(n-1)//2
    badI = totI = badN = totN = 0
    for m in range(1, maxm+1):
        # --- LS side ---
        for al in comps_le_delta(n, m):
            H = {tuple([0]*N): 1}
            for k in range(1, n):
                if al[k-1]: H = pmul(H, hpoly(al[k-1], k, N))
            for u in perms:
                if length(u)+m > maxm: continue
                P = pmul(tab[u], H)
                I = {}
                for w in byl[length(u)+m]:
                    I[w] = I_table(u, w, n).get(al, 0)
                for w in byl[length(u)+m]:
                    totI += 1
                    if extract_coeff(P, w, N) != I[w]:
                        badI += 1
                        if badI == 1: print("   PIERI/LS FAILS", u, w, al)
        # --- Samuel side ---
        for be in comps_le_delta(n, m):
            Y = {tuple([0]*N): 1}
            for k in range(1, n):
                Yk = {tuple([1 if i == j else 0 for i in range(N)]): 1 for j in range(k)}
                for _ in range(be[k-1]): Y = pmul(Y, Yk)
            for u in perms:
                if length(u)+m > maxm: continue
                P = pmul(tab[u], Y)
                for w in byl[length(u)+m]:
                    totN += 1
                    if extract_coeff(P, w, N) != T_tensor(u, w, n).get(word(be), 0):
                        badN += 1
                        if badN == 1: print("   MONK/SAMUEL FAILS", u, w, be)
    print(f"n={n}:  S_u h_alpha = sum_w I_alpha(u,w) S_w : {totI-badI}/{totI}, {badI} bad")
    print(f"n={n}:  S_u Y^beta  = sum_w N_(w/u)(beta) S_w: {totN-badN}/{totN}, {badN} bad")

if __name__ == "__main__":
    import sys
    for n in [int(a) for a in sys.argv[1:]] or [3,4]:
        run(n)
