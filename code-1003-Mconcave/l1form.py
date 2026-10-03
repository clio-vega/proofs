"""TEST the reformulation:  Lambda(nu) = (n-m)/2 - (1/2) sum_i |nu_i - (lam_{i-1}+1)|.

Derivation.  With a_i = lam_{i-1}+1 and the identity (t)_+ = ((t)+|t|)/2,
    sum_i (nu_i - a_i)_+ = (1/2)(S_nu - sum a_i) + (1/2) sum_i |nu_i - a_i|,
and sum_i a_i = |lam| - n + m (cyclic: sum_i lam_{i-1} = |lam| - n).  Plug into
    Lambda = (1/2)S_nu - (1/2)|lam| + n - m - sum_i (nu_i - a_i)_+
and S_nu cancels:   Lambda = (n-m)/2 - (1/2) sum_i |nu_i - a_i|.
CALIBRATION first: tent.Lambda_formula is itself verified against len(co) in tent.py.
"""
import gen, tent
from collections import Counter

def a_vec(lam, n, m):
    return [tent.lam_at(lam, i-1, n, m) + 1 for i in range(m)]

def Lambda_l1(nu, lam, n, m):
    a = a_vec(lam, n, m)
    return (n - m)/2 - sum(abs(nu[i] - a[i]) for i in range(m))/2

st = Counter(); wit = []
# calibration: sum_i a_i == |lam| - n + m
for m in range(2, 7):
    for n in range(m, m+5):
        for (mu, lam) in gen.pairs(n, m, 7):
            if sum(a_vec(lam,n,m)) != sum(lam) - n + m:
                st['A_BAD'] += 1
            else:
                st['A_ok'] += 1
            for nu in gen.box(mu, n, m):
                st['nu'] += 1
                v1 = tent.Lambda_formula(nu, lam, n, m)
                v2 = Lambda_l1(nu, lam, n, m)
                if abs(v1 - v2) < 1e-9: st['L_ok'] += 1
                else:
                    st['L_BAD'] += 1
                    if len(wit) < 4: wit.append((m,n,mu,lam,nu,v1,v2))
                # also: against the TRUE half-width (len(co)-1)/2, when nonzero
                off, co = gen.f_nu(nu, lam, n, m)
                if co:
                    if abs(v2 - (len(co)-1)/2) < 1e-9: st['L_vs_true_ok'] += 1
                    else:
                        st['L_vs_true_BAD'] += 1
                        if len(wit) < 8: wit.append(('TRUE',m,n,mu,lam,nu,v2,(len(co)-1)/2))
print(dict(st))
for w in wit: print('  ', w)
