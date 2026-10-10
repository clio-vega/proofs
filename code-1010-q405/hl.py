"""Independent engine: Hall-Littlewood P_lambda via Macdonald's symmetrisation formula
   P_la(x;t) = (1/v_la(t)) * sum_{w in S_n} w( x^la prod_{i<j} (x_i - t x_j)/(x_i - x_j) )
   Macdonald, Symmetric Functions and Hall Polynomials, III (1.2)/(2.2).
   v_la(t) = prod_{i>=0} [m_i(la)]_t!  with m_0 counted over the n-|support| zero parts.
   Cross-checked below against P_la(x;0)=s_la and P_la(x;1)=m_la."""
import sympy as sp, itertools
from sympy import Symbol, factor, cancel, simplify
t = sp.Symbol('t')

def qfac(m, t):
    r = sp.Integer(1)
    for i in range(1, m+1):
        r *= (1 - t**i)/(1 - t)
    return sp.simplify(r)

def v_la(la, n, t):
    la = list(la) + [0]*(n-len(la))
    from collections import Counter
    c = Counter(la)
    r = sp.Integer(1)
    for part, mult in c.items():
        r *= qfac(mult, t)
    return sp.simplify(r)

def P(la, xs, t=t):
    n = len(xs)
    lam = list(la) + [0]*(n-len(la))
    assert len(lam) == n
    tot = sp.Integer(0)
    for w in itertools.permutations(range(n)):
        y = [xs[w[i]] for i in range(n)]
        term = sp.prod([y[i]**lam[i] for i in range(n)])
        for i in range(n):
            for j in range(i+1, n):
                term *= (y[i] - t*y[j])/(y[i] - y[j])
        tot += term
    return sp.cancel(sp.together(sp.simplify(tot/v_la(lam, n, t))))
