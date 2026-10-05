"""The transition between the two coordinate systems.

Degree-m part of H^*(Fl_n):
   Samuel coordinates:  N_{w/u}(p) = <S_u Y_p, S_w>,  Y_p = prod_s (x_1+..+x_{p_s})
   LS coordinates:      I_alpha(u,w) = <S_u h_alpha, S_w>,
                        h_alpha = prod_k h_{alpha_k}(x_1..x_k)

Both families expanded at u = e give square/rectangular matrices
   I_mat[alpha][v] = I_alpha(e,v)          (square: Lehmer code bijection)
   N_mat[beta][v]  = N_v(p_beta)
The transition eta = N_mat . I_mat^{-1} is then used UNCHANGED for every (u,w).
"""
from fractions import Fraction
from collections import Counter
from itertools import permutations, combinations_with_replacement
from schub import length
from chains import T_tensor
from ls_side import I_table


def comps_le_delta(n, m):
    """compositions alpha with n-1 parts, |alpha|=m, alpha_k <= n-k."""
    out = []
    def rec(k, rem, acc):
        if k == n:
            if rem == 0: out.append(tuple(acc))
            return
        for a in range(0, min(n-k, rem)+1):
            rec(k+1, rem-a, acc+[a])
    rec(1, m, [])
    return out


def word(alpha):
    return tuple(sum([[k+1]*alpha[k] for k in range(len(alpha))], []))


def inv_matrix(M):
    """exact inverse of a square matrix of ints (list of lists) over Q."""
    N = len(M)
    A = [[Fraction(M[i][j]) for j in range(N)] + [Fraction(int(i == j)) for j in range(N)]
         for i in range(N)]
    for c in range(N):
        piv = next((r for r in range(c, N) if A[r][c] != 0), None)
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [x/pv for x in A[c]]
        for r in range(N):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [A[r][j]-f*A[c][j] for j in range(2*N)]
    return [row[N:] for row in A]


def run(n):
    perms = list(permutations(range(1, n+1)))
    e = tuple(range(1, n+1))
    maxm = n*(n-1)//2
    byl = {}
    for v in perms:
        byl.setdefault(length(v), []).append(v)

    print(f"=== n = {n} ===")
    for m in range(1, maxm+1):
        A = comps_le_delta(n, m)
        V = byl[m]
        assert len(A) == len(V), (m, len(A), len(V))     # Lehmer-code bijection
        # --- I_mat at u = e ---
        Ie = {v: I_table(e, v, n) for v in V}
        I_mat = [[Ie[v].get(al, 0) for v in V] for al in A]
        Iinv = inv_matrix(I_mat)
        if Iinv is None:
            print(f"  m={m}: I_mat SINGULAR -- h_alpha is NOT a basis")
            continue
        # determinant via integrality of the inverse:
        unimod = all(x.denominator == 1 for row in Iinv for x in row)
        # --- N_mat at u = e ---
        BETAS = []
        for c in combinations_with_replacement(range(1, n), m):
            b = [0]*(n-1)
            for k in c: b[k-1] += 1
            BETAS.append(tuple(b))
        Te = {v: T_tensor(e, v, n) for v in V}
        N_mat = [[Te[v].get(word(be), 0) for v in V] for be in BETAS]
        # --- eta = N_mat . Iinv ---
        eta = [[sum(Fraction(N_mat[i][k])*Iinv[k][j] for k in range(len(V)))
                for j in range(len(A))] for i in range(len(BETAS))]
        eta_int = all(x.denominator == 1 for row in eta for x in row)

        # --- the test: eta is universal, i.e. works for every (u,w) ---
        bad = tot = 0; nontrivial = 0
        for u in perms:
            lu = length(u)
            if lu+m > maxm: continue
            for w in byl.get(lu+m, []):
                T = T_tensor(u, w, n)
                if not T and length(w)-length(u) != 0:
                    pass
                I = I_table(u, w, n)
                if not T and not I: continue
                for i, be in enumerate(BETAS):
                    lhs = T.get(word(be), 0)
                    rhs = sum(eta[i][j]*I.get(A[j], 0) for j in range(len(A)))
                    tot += 1
                    if lhs != rhs:
                        bad += 1
                        if bad == 1:
                            print(f"    UNIVERSALITY FAILS u={u} w={w} beta={be} lhs={lhs} rhs={rhs}")
                    if lhs != 0 and sum(1 for j in range(len(A)) if eta[i][j] != 0 and I.get(A[j],0) != 0) > 1:
                        nontrivial += 1
        print(f"  m={m}: |A|=|V|={len(A)}, |beta|={len(BETAS)}, "
              f"I_mat unimodular over Z: {unimod}, eta integral: {eta_int}; "
              f"universality {tot-bad}/{tot} ({nontrivial} with >=2 terms on the right), {bad} bad")


if __name__ == "__main__":
    for n in (3, 4):
        run(n)
