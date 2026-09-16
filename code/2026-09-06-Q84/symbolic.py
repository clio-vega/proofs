"""symbolic.py -- Corollary B, checked directly in t.

ord(N_e) = infinity forces ord_{Q(t)}(R_e(t)) = infinity, because
ad^d commutes with d/dt (the M_{p_a} are t-free) and R'_e(-1) = -Psi_{e,1}
= -(N_e - M[...]).  Here we watch it happen: compute the hook coefficient
    c_{e,d}(t) := < s_{(2,1^{d+e-2})} , ad_{M_{p_1}}^d ( R_e(t) ) s_empty >
symbolically and check
    (i)  c_{e,d}(-1) = 0                  (R_e(-1) = M_{p_e} has order 0)
    (ii) c'_{e,d}(-1) = -(-1)^{d+e+1}     (the Theorem A witness, with sign)
    (iii) c_{e,d}(t) is not identically 0 -> ord R_e(t) = infinity.
Also E_e(t) = g_e - R_e: same coefficient with the opposite eps-linear part.
"""
import sympy as sp
from math import comb
from engine import t, parts_of, op_R, op_p, op_g, clean

def ad_d(op, e, d, lam=()):
    out = {}
    for r in range(d + 1):
        c = (-1)**r * comb(d, r)
        st = {lam: sp.Integer(1)}
        for _ in range(d - r): st = op_p(st, 1)
        st = op(st, e)
        for _ in range(r): st = op_p(st, 1)
        for k, v in st.items(): out[k] = out.get(k, 0) + c * v
    return clean(out)

print(' e  d   c_{e,d}(t) = <s_(2,1^{d+e-2}), ad^d R_e(t) . 1>          c(-1)  c\'(-1)  pred')
for e in (2, 3, 4):
    for d in (1, 2, 3):
        mu = (2,) + (1,) * (d + e - 2)
        c = sp.expand(ad_d(op_R, e, d).get(mu, sp.Integer(0)))
        c0 = sp.simplify(c.subs(t, -1)) if c != 0 else 0
        c1 = sp.simplify(sp.diff(c, t).subs(t, -1)) if c != 0 else 0
        print('%2d %2d   %-52s %5s  %5s  %5d' %
              (e, d, sp.factor(c), c0, c1, -(-1)**(d + e + 1)))
