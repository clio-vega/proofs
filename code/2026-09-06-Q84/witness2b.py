"""witness2b.py -- Theorem 2 closed form, tested at integer t over a wide range.

    < s_{(e,1^d)} , ad_{M_p1}^d ( R_e(t) ) . 1 >
        = (-1)^d [ sum_{b=0}^{min(d,e-1)} C(d,b) t^b  -  C(d-1,e-1) t^{e-1} ].

Mechanism (proved by hand):  inside mu = (e,1^d) a connected e-ribbon kappa/nu is
either (i) nu = empty, kappa = (e-b,1^b)   [it must touch the corner (1,1)], or
(ii) nu = (a,1^{b'}), kappa = (a,1^{b'+e}) [a pure leg extension].
Case (ii) sums by Vandermonde to the CONSTANT C(d-1,e-1) t^{e-1}, independent
of r -- so the d-th finite difference annihilates it, EXCEPT that Vandermonde
fails at s = d-r = 0, i.e. exactly at r = d.  Case (i) lives only at r = d.
So F(r) = const + spike at r=d, and ad^d keeps only the spike.
"""
from math import comb
from fast import op_p, ribbons

def op_R_int(state, e, tv):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, e):
            out[mu] = out.get(mu, 0) + w * tv**h
    return {k: v for k, v in out.items() if v}

def ad_d(e, d, tv, mu):
    out = 0
    for r in range(d + 1):
        st = {(): 1}
        for _ in range(d - r): st = op_p(st, 1)
        st = op_R_int(st, e, tv)
        for _ in range(r): st = op_p(st, 1)
        out += (-1)**r * comb(d, r) * st.get(mu, 0)
    return out

def pred(e, d, tv):
    s = sum(comb(d, b) * tv**b for b in range(0, min(d, e - 1) + 1))
    return (-1)**d * (s - comb(d - 1, e - 1) * tv**(e - 1))

bad = tot = 0
for e in range(2, 7):
    for d in range(1, 8):
        if d + e > 13: continue
        mu = (e,) + (1,) * d
        for tv in (-3, -2, -1, 0, 1, 2, 3, 5):
            g, p = ad_d(e, d, tv, mu), pred(e, d, tv)
            tot += 1
            if g != p:
                bad += 1
                if bad < 5: print('MISMATCH e=%d d=%d t=%d got %s pred %s' % (e, d, tv, g, p))
print('Theorem 2 closed form: %d/%d agree' % (tot - bad, tot))
print('value at t=0 is (-1)^d :',
      all(ad_d(e, d, 0, (e,) + (1,) * d) == (-1)**d
          for e in range(2, 7) for d in range(1, 7) if d + e <= 12))
