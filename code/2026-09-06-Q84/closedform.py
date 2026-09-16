"""closedform.py -- the hook coefficient of ad^d R_e(t), at integer t, larger range.

Conjecture from the symbolic run (e<=4, d<=3):
    < s_{(2,1^{d+e-2})} , ad_{M_p1}^d ( R_e(t) ) . 1 >  =  (-1)^d t^{e-2} (1+t).
Independent of d except for the sign.  Tested here at t = -3..4 for e=2..6,
d=1..7 using an integer-arithmetic engine (no sympy).
"""
from math import comb
from fast import op_p, ribbons

def op_R_int(state, e, tv):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, e):
            out[mu] = out.get(mu, 0) + w * tv**h
    return {k: v for k, v in out.items() if v}

def ad_d(e, d, tv):
    out = {}
    for r in range(d + 1):
        c = (-1)**r * comb(d, r)
        st = {(): 1}
        for _ in range(d - r): st = op_p(st, 1)
        st = op_R_int(st, e, tv)
        for _ in range(r): st = op_p(st, 1)
        for k, v in st.items(): out[k] = out.get(k, 0) + c * v
    return out

bad = 0; tot = 0
for e in range(2, 7):
    for d in range(1, 8):
        if d + e > 13: continue
        mu = (2,) + (1,) * (d + e - 2)
        for tv in (-3, -2, -1, 0, 1, 2, 4):
            got = ad_d(e, d, tv).get(mu, 0)
            pred = (-1)**d * tv**(e - 2) * (1 + tv)
            tot += 1
            if got != pred:
                bad += 1
                if bad < 5: print('MISMATCH e=%d d=%d t=%d got %d pred %d' % (e, d, tv, got, pred))
print('hook coeff of ad^d R_e(t) == (-1)^d t^{e-2}(1+t) : %d/%d' % (tot - bad, tot))
