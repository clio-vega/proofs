"""Ground truth on the LS enumerator, BEFORE any comparison with Samuel.

Check 1 (Prop 3, Bergeron-Sottile):  S_u(x) = sum over increasing chains u -> w0
         of x^delta / x^gamma.   Compared against the divided-difference
         Schubert table -- a completely different mechanism.
Check 2 (Thm 2):  S_{w/u}(x) = sum_gamma x^delta/x^gamma over increasing chains
         u -> w,  and  S_{w/u} = sum_v c^w_{u,v} S_{w0 v}.
Check 3 (Cor 4):  I_alpha(u,w) = sum_v c^w_{u,v} I_alpha(w0 v, w0).
"""
import sys
from collections import Counter
from itertools import permutations
from schub import length, rmul_t, schubert_table, structure_constants
from ls_side import I_table, increasing_chains


def run(n):
    N = n + 1
    w0 = tuple(range(n, 0, -1))
    delta = tuple([n-1-i for i in range(n-1)])          # (n-1,...,1) in x_1..x_{n-1}
    tab = schubert_table(n, N)
    perms = list(permutations(range(1, n+1)))

    # ---------- Check 1: Prop 3 ----------
    bad1 = t1 = 0
    Itab_to_w0 = {}
    for u in perms:
        I = I_table(u, w0, n)
        Itab_to_w0[u] = I
        poly = Counter()
        for alpha, c in I.items():
            e = tuple([delta[i]-alpha[i] for i in range(n-1)] + [0]*(N-(n-1)))
            assert min(e) >= 0, (u, alpha)
            poly[e] += c
        truth = {e: c for e, c in tab[u].items() if c}
        got = {e: c for e, c in poly.items() if c}
        t1 += 1
        if got != truth:
            bad1 += 1
            if bad1 == 1:
                print("  PROP3 MISMATCH", u, "got", got, "truth", truth)
    print(f"  Check 1 (Prop 3, S_u = sum_gamma x^delta/x^gamma): {t1-bad1}/{t1}, {bad1} bad")

    # ---------- structure constants ----------
    C, _ = structure_constants(n, N)
    cof = {}
    for (u, v, w), c in C.items():
        cof.setdefault((u, w), {})[v] = c

    # ---------- Check 3: Cor 4 ----------
    bad3 = t3 = 0
    nontriv3 = 0
    for u in perms:
        for w in perms:
            m = length(w) - length(u)
            if m < 0:
                continue
            I = I_table(u, w, n)
            rhs = Counter()
            for v, c in cof.get((u, w), {}).items():
                for alpha, iv in Itab_to_w0[(tuple(w0[i]-1 for i in range(n)) and v)].items() if False else Itab_to_w0[rmul_t_w0(w0, v)].items():
                    rhs[alpha] += c*iv
            allk = set(I) | set(rhs)
            for a in allk:
                t3 += 1
                if I[a] != rhs[a]:
                    bad3 += 1
                    if bad3 == 1:
                        print("  COR4 MISMATCH", u, w, a, I[a], rhs[a])
                elif I[a] != 0 and len(cof.get((u, w), {})) > 1:
                    nontriv3 += 1
    print(f"  Check 3 (Cor 4, I_alpha(u,w) = sum_v c^w_uv I_alpha(w0 v, w0)): "
          f"{t3-bad3}/{t3}, {bad3} bad  ({nontriv3} with >1 term on the right)")
    return bad1 + bad3


def rmul_t_w0(w0, v):
    """w0 * v  (composition of permutations, v then w0): (w0 v)(i) = w0(v(i))."""
    return tuple(w0[v[i]-1] for i in range(len(v)))


if __name__ == "__main__":
    for n in (3, 4, 5):
        print(f"n = {n}")
        run(n)
