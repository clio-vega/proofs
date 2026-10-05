"""Proposition 6.1 (the discriminator) and the block-factorial refutation.
  LS single block:     I_{a e_k}(u,w) = <S_u h_a(x_1..x_k), S_w>   -- Pieri
  Samuel single block: N_{w/u}(k^a)   = <S_u Y_k^a, S_w>           -- iterated Monk
Pieri is multiplicity-free, Monk is not; hence no scalar factor relates them."""
from collections import Counter
from itertools import permutations
from math import factorial
from schub import length
from chains import T_tensor
from ls_side import I_table

for n in (3, 4, 5):
    perms = list(permutations(range(1, n+1)))
    Iv = Counter(); Nv = Counter(); tested = strict = 0
    fok = fbad = 0
    for u in perms:
        for w in perms:
            m = length(w)-length(u)
            if m < 1: continue
            T = T_tensor(u, w, n); I = I_table(u, w, n)
            if not T and not I: continue
            for k in range(1, n):
                if m > n-k: continue
                al = tuple(m if j == k-1 else 0 for j in range(n-1))
                iv = I.get(al, 0); nv = T.get(tuple([k]*m), 0)
                Iv[iv] += 1; Nv[nv] += 1; tested += 1
                assert iv <= nv, (u, w, al, iv, nv)
                if iv < nv: strict += 1
                if m >= 2:
                    if nv == factorial(m)*iv: fok += 1
                    else: fbad += 1
    print(f"n={n}: {tested} single-block tests")
    print(f"   LS/Pieri  I_(a e_k) values: {dict(sorted(Iv.items()))}   <- multiplicity-free iff keys subset of {{0,1}}")
    print(f"   Samuel    N_(w/u)(k^a) values: {dict(sorted(Nv.items()))}")
    print(f"   I < N in {strict} cases; block-factorial guess N(k^a)=a!*I fails {fbad} of {fok+fbad}")
