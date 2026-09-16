"""Phase 0a -- the dictionary check, with its normalisation, in ONE convention.

Claim (derived, not fitted):  with q = -t and eps = 1+t,
   P_{(e)}(X;q) = q_e(X;q)/(1-q),   sum_r q_r(X;q) u^r = prod_i (1-q x_i u)/(1-x_i u)
                = exp( sum_n (1-q^n) p_n u^n / n ),
hence
   P_{(e)}(X;-t) = sum_{rho |- e} eps^{ell(rho)-1} / z_rho * prod_i [rho_i]_{-t} * p_rho,
   [n]_{-t} = 1 + (-t) + ... + (-t)^{n-1}.

TEST: M[P_{(e)}(X;-t)] s_lam  ==  g_e s_lam  (broken border strips, (1+t)^{m-1} t^{ht}).
No free scalar anywhere: this is a check, not a fit.
"""
from engine import *
import sympy as sp

def z(rho):
    from collections import Counter
    c = Counter(rho); r = 1
    for i, m in c.items(): r *= (i**m) * sp.factorial(m)
    return sp.Integer(r)

def bracket_n(n):
    return sum((-t)**i for i in range(n))

def P_row(e):
    """P_{(e)}(X;-t) in the p-basis: dict rho -> coefficient."""
    out = {}
    for rho in parts_of(e):
        c = (1 + t)**(len(rho) - 1) / z(rho)
        for r in rho: c *= bracket_n(r)
        out[rho] = sp.expand(sp.simplify(c))
    return out

def mult_p_rho(state, rho):
    for a in rho: state = op_p(state, a)
    return state

def M_P(state, e):
    out = {}
    for rho, c in P_row(e).items():
        for mu, v in mult_p_rho(state, rho).items():
            out[mu] = out.get(mu, 0) + c * v
    return clean(out)

if __name__ == '__main__':
    for e in range(1, 6):
        print('P_(%d)(X;-t) =' % e,
              ' + '.join('(%s) p_%s' % (sp.factor(c), ''.join(map(str, r)))
                         for r, c in sorted(P_row(e).items()) if c != 0))
    print()
    bad = 0; tot = 0
    for e in range(1, 6):
        for n in range(0, 7):
            for lam in parts_of(n):
                A = M_P({lam: sp.Integer(1)}, e)
                B = op_g({lam: sp.Integer(1)}, e)
                tot += 1
                d = sub(A, B)
                if d:
                    bad += 1
                    if bad < 4: print('MISMATCH', e, lam, d)
    print('M[P_(e)(X;-t)] == g_e :  %d/%d agree' % (tot - bad, tot))
