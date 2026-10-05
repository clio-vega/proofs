"""Independent confirmation of Lenart-Sottile's Lemma 5, used in Theorem 5.4.

Derivation (S5.2 of the paper): h_alpha = sum_v I_alpha(e,v) S_v by Lemma 4.2 at
u=e, plus Poincare duality <S_a S_b, S_w0> = delta_{b, w0 a} with a = w0 v, gives
psi_alpha(S_{w0 v}) = I_alpha(e,v).  LS Lemma 5 says psi_alpha(f) is the
coefficient of x^{delta-alpha} in the normal form of f.  So the prediction is

        I_alpha(e,v)  =  [ x^{delta-alpha} ] S_{w0 v},

the left side from the increasing-chain enumerator, the right from the
divided-difference Schubert table.  Two different mechanisms."""
from itertools import permutations
from schub import length, schubert_table
from ls_side import I_table
from transition import comps_le_delta

for n in (3, 4, 5):
    N = n+1
    tab = schubert_table(n, N)
    w0 = tuple(range(n, 0, -1))
    delta = [n-1-i for i in range(n-1)]
    e = tuple(range(1, n+1))
    bad = tot = 0
    for v in permutations(range(1, n+1)):
        m = length(v)
        w0v = tuple(w0[v[i]-1] for i in range(n))
        Iv = I_table(e, v, n)
        for al in comps_le_delta(n, m):
            ex = tuple([delta[i]-al[i] for i in range(n-1)] + [0]*(N-(n-1)))
            tot += 1
            if Iv.get(al, 0) != tab[w0v].get(ex, 0):
                bad += 1
                if bad == 1: print("   LEM5 MISMATCH", v, al)
    print(f"n={n}: LS Lemma 5, I_alpha(e,v) = [x^(delta-alpha)] S_(w0 v): {tot-bad}/{tot}, {bad} bad")
