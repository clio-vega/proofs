"""PLANTED NEGATIVE CONTROLS on the LS enumerator.

Claim the instrument makes: "I enumerate exactly the lex-STRICTLY-increasing
chains in the labelled Bruhat order, with b = u(i)."
Three ways to break it, one per ingredient.  Each must turn Prop 3 and Cor 4 RED.
"""
from collections import Counter
from itertools import permutations
from schub import length, rmul_t, schubert_table, structure_constants
import ls_side

GOOD_edges = ls_side.labelled_edges

def make_broken(mode):
    def edges(u, n):
        out = []; lu = length(u)
        for i in range(1, n+1):
            for j in range(i+1, n+1):
                w = rmul_t(u, i, j)
                if length(w) != lu+1: continue
                if mode == "b_wrong":
                    b = u[j-1]                 # the OTHER convention, b=u(j)
                else:
                    b = u[i-1]
                rng = range(i, j) if mode != "k_wrong" else range(i, j+1)
                for k in rng:
                    if k > n-1: continue
                    out.append((k, b, w))
        return out
    return edges

def check(n):
    N = n+1
    w0 = tuple(range(n, 0, -1))
    delta = tuple(n-1-i for i in range(n-1))
    tab = schubert_table(n, N)
    perms = list(permutations(range(1, n+1)))
    bad = 0
    for u in perms:
        poly = Counter()
        for alpha, c in ls_side.I_table(u, w0, n).items():
            e = tuple([delta[i]-alpha[i] for i in range(n-1)] + [0]*(N-(n-1)))
            if min(e) < 0: return "NEGATIVE EXPONENT (fires)"
            poly[e] += c
        if {e: c for e, c in poly.items() if c} != {e: c for e, c in tab[u].items() if c}:
            bad += 1
    return f"{bad} Prop-3 mismatches"

n = 4
print("baseline (correct enumerator)      :", check(n))

for mode, desc in [("weak", "lex condition weakened to NON-strict (allow repeated labels)"),
                   ("b_wrong", "second coordinate taken as b = u(j) instead of u(i)"),
                   ("k_wrong", "label range widened to i <= k <= j")]:
    if mode == "weak":
        orig = ls_side.increasing_chains
        def weak(u, w, n, _o=orig):
            target = length(w)
            def rec(v, last):
                if v == w: yield []; return
                if length(v) >= target: return
                for (k, b, v2) in ls_side.labelled_edges(v, n):
                    if (k, b) < last: continue          # NON-strict: bug
                    for rest in rec(v2, (k, b)): yield [(k, b)]+rest
            yield from rec(u, (0, 0))
        ls_side.increasing_chains = weak
        print(f"control '{desc}':", check(n))
        ls_side.increasing_chains = orig
    else:
        ls_side.labelled_edges = make_broken(mode)
        print(f"control '{desc}':", check(n))
        ls_side.labelled_edges = GOOD_edges
print("baseline again (restored)          :", check(n))
