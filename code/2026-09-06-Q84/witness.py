"""witness.py -- the hook witness.

CLAIM (proved by hand, tested here).  Let ad(X) := [X, M_{p_1}].  For every
e >= 2 and every d >= 1,

    < s_{(2,1^{d+e-2})} ,  ad^d(N_e) . s_empty >  =  (-1)^{d+e+1}   (never 0).

Mechanism: inside the arm-1 hook mu = (2,1^m) the ONLY two-component broken
border strips kappa/nu are  nu = (1^{j'}) -> kappa = (2,1^{j'+e-2}) with j' >= 1
(j' = 0 makes the cell (1,1) join the leg, giving one component).  Each such
has f^nu = 1, f^{mu/kappa} = 1, ht = e-2.  So the coefficient is
    (-1)^e * sum_{r=0}^{d-1} (-1)^r C(d,r)  =  (-1)^e * ( -(-1)^d )  =  (-1)^{d+e+1}.
The alternating sum is truncated at r = d-1 exactly because j' = 0 is excluded:
the missing r = d term is what makes the sum nonzero. That single excluded
chain is the whole theorem.
"""
import itertools, sys
from fast import op_p, op_N

def nested_p1(e, d, lam=()):
    out = {}
    for r in range(d + 1):
        c = (-1)**r
        from math import comb
        c *= comb(d, r)
        st = {lam: 1}
        for _ in range(d - r): st = op_p(st, 1)
        st = op_N(st, e)
        for _ in range(r): st = op_p(st, 1)
        for k, v in st.items(): out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}

if __name__ == '__main__':
    print(' e   d   mu                       coeff   predicted (-1)^(d+e+1)')
    ok = True
    for e in range(2, 7):
        for d in range(1, 9):
            if d + e > 14: continue
            mu = (2,) + (1,) * (d + e - 2)
            got = nested_p1(e, d).get(mu, 0)
            pred = (-1)**(d + e + 1)
            flag = 'OK' if got == pred else '*** MISMATCH ***'
            ok &= (got == pred)
            print('%2d %3d   %-22s %5d   %5d  %s' % (e, d, str(mu), got, pred, flag))
        sys.stdout.flush()
    print('\nALL AGREE' if ok else '\nFAILED')
