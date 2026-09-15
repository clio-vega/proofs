"""
The conjectured family: for each d | gcd(e,f), transport Murnaghan-Nakayama
through the d-quotient.

  MNd(e,d)(u) = (-1)^{#{ i in {d,2d,...,e-d} : u_i = 1 }}      u in {0,1}^{e-1}

i.e. the sign of the number of jumped beads lying on the SAME d-runner as b.
d=1 gives the honest MN sign (-1)^{height}; d=e gives the constant weight 1.
"""
import sympy as sp
from engine import all_words, word_to_comp

def MNd(e, d):
    """dict word -> +-1"""
    assert e % d == 0
    return {u: sp.Integer(-1)**sum(u[i-1] for i in range(d, e, d)) for u in all_words(e-1)}

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]


def gamma_weight(e, d, gam):
    """W(u) = (-1)^{|u|} * prod_{i=1..e-1} gam[i mod d]^{u_i}.
    gam: dict on Z/d with gam[0] = 1.  d must divide e."""
    assert e % d == 0 and sp.simplify(gam[0]-1) == 0
    out = {}
    for u in all_words(e-1):
        val = sp.Integer(-1)**sum(u)
        for i, ui in enumerate(u, start=1):
            if ui: val *= gam[i % d]
        out[u] = sp.simplify(val)
    return out
