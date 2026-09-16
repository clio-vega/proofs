"""Phase 0b/0c.

(0b)  N_e - Psi_{e,1} = M[ d/dt P_{(e)}(X;-t) |_{t=-1} ]
                      = M[ (1/2) sum_{a=1}^{e-1} p_a p_{e-a}  -  ((e-1)/2) p_e ].
      Derived from  [n]_{-t} = n - eps*C(n,2) + O(eps^2),  eps = 1+t:
        P_(e) = sum_rho eps^{l(rho)-1}/m_rho! * (1 - eps*(e-l(rho))/2 + O(eps^2)) p_rho.
      This is the exact statement that "component number is only a FILTRATION":
      the eps-linear part of g_e is NOT the 2-component operator alone; it is
      the 2-component operator MINUS the height-weighted connected operator.
      The mismatch Psi_{e,1} is precisely the merging term.

(0c)  N_e |empty> = 0 ?    (must be checked, not assumed)
"""
from fast import *
import itertools

def M_correction(state, e):
    """M[ (1/2) sum_{a<e} p_a p_{e-a} - ((e-1)/2) p_e ] with integer bookkeeping x2."""
    out = {}
    for a in range(1, e):
        for mu, v in op_p(op_p(state, a), e - a).items():
            out[mu] = out.get(mu, 0) + v          # this is 2 * (1/2) sum
    for mu, v in op_p(state, e).items():
        out[mu] = out.get(mu, 0) - (e - 1) * v
    return {k: v for k, v in out.items() if v}

def sub(a, b):
    o = dict(a)
    for k, v in b.items(): o[k] = o.get(k, 0) - v
    return {k: v for k, v in o.items() if v}

def parts(n, cap=None):
    if cap is None: cap = n
    if n == 0: return [()]
    return [(p,) + r for p in range(min(n, cap), 0, -1) for r in parts(n - p, p)]

if __name__ == '__main__':
    bad = 0; tot = 0
    for e in range(2, 7):
        for n in range(0, 9):
            for lam in parts(n):
                st = {lam: 1}
                lhs = sub(op_N(st, e), op_Psi(st, e))
                rhs = {k: v for k, v in
                       ((k, v) for k, v in M_correction(st, e).items())}
                # lhs must equal rhs/2
                d = sub({k: 2 * v for k, v in lhs.items()}, rhs)
                tot += 1
                if d:
                    bad += 1
                    if bad < 4: print('MISMATCH e=%d lam=%s' % (e, lam), d)
    print('(0b)  N_e - Psi_{e,1} = M[dP/dt|_{-1}] :  %d/%d agree' % (tot - bad, tot))

    print()
    print('(0c)  N_e |empty> :', op_N({(): 1}, 3), op_N({(): 1}, 4), op_N({(): 1}, 5))
    print('      N_e |(1)>   :', op_N({(1,): 1}, 3), '   N_e |(2)> :', op_N({(2,): 1}, 3))
