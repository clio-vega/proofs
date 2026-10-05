"""The ELEMENTARY direction:  LS from Samuel, with no basis theorem.

In Z[x_1,...,x_{n-1}] put Y_p = x_1 + ... + x_p (so x_p = Y_p - Y_{p-1}, Y_0=0).
Then h_alpha = prod_k h_{alpha_k}(x_1..x_k) is an integral polynomial in the Y's:
        h_alpha = sum_beta chi_alpha(beta) Y^beta,     chi_alpha(beta) in Z,
an identity of POLYNOMIALS (no quotient, no cohomology, no triangularity).
Applying the functional Lambda_{w/u}(f) = <S_u f, S_w> gives
        I_alpha(u,w) = sum_beta chi_alpha(beta) N_{w/u}(p_beta).
chi depends on neither u nor w.  This script computes chi and tests that.
"""
from collections import Counter
from itertools import permutations
from schub import length
from chains import T_tensor
from ls_side import I_table
from transition import comps_le_delta, word

# polynomials in Y_1..Y_{n-1}: dict  exponent-tuple -> int
def ymul(f, g):
    h = Counter()
    for e1, c1 in f.items():
        for e2, c2 in g.items():
            h[tuple(a+b for a, b in zip(e1, e2))] += c1*c2
    return {e: c for e, c in h.items() if c}

def yadd(f, g):
    h = Counter(f); h.update(g)
    return {e: c for e, c in h.items() if c}

def yscale(f, s):
    return {e: c*s for e, c in f.items() if c*s}

def Yvar(p, n):
    """the generator Y_p as a monomial."""
    e = [0]*(n-1); e[p-1] = 1
    return {tuple(e): 1}

def xvar(p, n):
    """x_p = Y_p - Y_{p-1}."""
    f = Yvar(p, n)
    return f if p == 1 else yadd(f, yscale(Yvar(p-1, n), -1))

def h_in_Y(a, k, n):
    """h_a(x_1,...,x_k) as a polynomial in Y_1..Y_{n-1}."""
    # recursion h_a(x_[k]) = sum_{r=0}^{a} x_k^r h_{a-r}(x_[k-1]),  h_a(x_[0]) = [a==0]
    cache = {}
    def H(a, k):
        if a == 0: return {tuple([0]*(n-1)): 1}
        if k == 0: return {}
        if (a, k) in cache: return cache[(a, k)]
        xk = xvar(k, n); out = {}; pw = {tuple([0]*(n-1)): 1}
        for r in range(a+1):
            out = yadd(out, ymul(pw, H(a-r, k-1)))
            pw = ymul(pw, xk)
        cache[(a, k)] = out
        return out
    return H(a, k)

def chi(alpha, n):
    """h_alpha = prod_k h_{alpha_k}(x_1..x_k), expanded in Y-monomials."""
    out = {tuple([0]*(n-1)): 1}
    for k in range(1, n):
        if alpha[k-1]:
            out = ymul(out, h_in_Y(alpha[k-1], k, n))
    return out

def run(n):
    perms = list(permutations(range(1, n+1)))
    maxm = n*(n-1)//2
    tot = bad = 0; nonunit = 0; terms = Counter()
    CH = {}
    for m in range(1, maxm+1):
        for al in comps_le_delta(n, m):
            CH[al] = chi(al, n)
            terms[len(CH[al])] += 1
            if len(CH[al]) > 1: nonunit += 1
    for u in perms:
        for w in perms:
            m = length(w)-length(u)
            if m < 1: continue
            T = T_tensor(u, w, n); I = I_table(u, w, n)
            if not T and not I: continue
            for al in comps_le_delta(n, m):
                lhs = I.get(al, 0)
                rhs = sum(c*T.get(word(be), 0) for be, c in CH[al].items())
                tot += 1
                if lhs != rhs:
                    bad += 1
                    if bad == 1: print("   CHI FAILS", u, w, al, lhs, rhs)
    print(f"n={n}: I_alpha(u,w) = sum_beta chi_alpha(beta) N_{{w/u}}(p_beta): "
          f"{tot-bad}/{tot}, {bad} bad;  {nonunit} of the chi_alpha have >1 term; "
          f"term-count histogram {dict(sorted(terms.items()))}")

if __name__ == "__main__":
    for n in (3, 4, 5):
        run(n)
