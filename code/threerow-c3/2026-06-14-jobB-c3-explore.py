#!/usr/bin/env python3
"""Job B exploration: look at c=3 M_j sequences and ratios to anchor a closed form."""
from fractions import Fraction as Fr
from job_jstar_engine import Mj_chain_all, C, v2, hook_f
from math import comb, factorial

def Mj_seq(lam):
    return Mj_chain_all(tuple(lam))

def dim3(alpha, beta, gamma):
    """f^{(alpha,beta,gamma)} via 3-row hook formula; 0 if not a partition."""
    if not (alpha >= beta >= gamma >= 0):
        return 0
    n = alpha+beta+gamma
    l1, l2, l3 = alpha+2, beta+1, gamma
    vand = (l1-l2)*(l1-l3)*(l2-l3)
    return factorial(n)*vand // (factorial(l1)*factorial(l2)*factorial(l3))

def finite_diff_degree(seq):
    """Return the degree of the polynomial interpolating seq (Fractions), or None if not poly."""
    s = [Fr(x) for x in seq]
    deg = 0
    cur = s[:]
    while len(cur) > 1 and any(x != 0 for x in cur):
        cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
        deg += 1
        if deg > len(seq)+2:
            return None
    return deg-1 if deg > 0 else 0

if __name__ == "__main__":
    print("=== c=3 sample sequences: M_j, dim_j=f^(a-j,b-j,3-j), ratio M_j/dim_j ===")
    for (a, b) in [(4, 3), (6, 3), (5, 4), (9, 6), (7, 4), (8, 5), (10, 7)]:
        c = 3
        lam = (a, b, c)
        if sum(lam) % 2 != 0:
            continue
        Ms = Mj_seq(lam)
        m = sum(lam)//2
        print(f"\n({a},{b},3) m={m}:")
        print(f"   M_j   = {Ms}")
        dims = [dim3(a-j, b-j, c-j) for j in range(m+1)]
        print(f"   dim_j = {dims}")
        ratios = []
        for j in range(m+1):
            if dims[j] != 0:
                ratios.append(Fr(Ms[j], dims[j]))
            else:
                ratios.append(None)
        print(f"   M/dim = {ratios}")

    # try: is M_j / [(a-b+1) * C(N, b-j)] a polynomial of degree 6 in j?
    print("\n=== test polynomiality of M_j / [(a-b+1) C(N,b-j)] over j in [0,b] ===")
    for (a, b) in [(6, 3), (9, 6), (8, 5), (10, 7)]:
        c = 3
        lam = (a, b, c); m = sum(lam)//2
        Ms = Mj_seq(lam)
        for Nname, N in [('a+2', a+2), ('a+1', a+1), ('m', m), ('a+3', a+3), ('2m-3j? no', a+2)]:
            quo = []
            ok = True
            for j in range(0, b+1):    # only j where b-j>=0
                base = (a-b+1)*comb(N, b-j)
                if base == 0:
                    ok = False; break
                quo.append(Fr(Ms[j], base))
            if not ok:
                continue
            deg = finite_diff_degree(quo)
            print(f"  ({a},{b},3) N={Nname}({N}): M/[(a-b+1)C(N,b-j)] poly-deg over j<=b = {deg}")
