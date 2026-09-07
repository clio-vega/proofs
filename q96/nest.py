"""Nested commutators ad(R_{m e}(-1)) applied to R_e(t)."""
import sympy as sp
from engine import R_abacus, t

def Rvec(vec, e, tt):
    out = {}
    for lam, c in vec.items():
        for mu, d in R_abacus(lam, e, tt).items():
            out[mu] = out.get(mu, 0) + c*d
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def ad(f_list, e, tt):
    """returns callable lam -> dict computing ad(R_{f_k}(-1))...ad(R_{f_1}(-1)) applied to R_e(tt),
    i.e. [[...[R_e(t), R_{f1}(-1)], R_{f2}(-1)],...]."""
    m1 = sp.Integer(-1)
    def X(lam):
        v = R_abacus(lam, e, tt)
        return {k: sp.expand(val) for k, val in v.items() if sp.expand(val) != 0}
    cur = X
    for f in f_list:
        prev = cur
        def make(prev, f):
            def g(lam):
                a = Rvec(prev(lam), f, m1)                 # R_f(-1) . prev
                b = prev_on(prev, R_abacus(lam, f, m1))    # prev . R_f(-1)
                out = dict(b)
                for k, v in a.items():
                    out[k] = sp.expand(out.get(k, 0) - v)
                return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}
            return g
        cur = make(prev, f)
    return cur

def prev_on(prev, vec):
    out = {}
    for lam, c in vec.items():
        for mu, d in prev(lam).items():
            out[mu] = out.get(mu, 0) + c*d
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}
