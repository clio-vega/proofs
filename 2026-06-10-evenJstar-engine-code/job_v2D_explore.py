"""
Explore v2(D_j) where D_j = 2^j M_j = sum_b C(j,b)(-1)^b chi_b,
chi_b = chi^lambda(2^b 1^{2m-2b}).  We want a formula for v2(D_j).

val(j) = 2 v2(C(m,j)) + 2 v2(D_j) - j.
val0(j) := 2 v2(D_j) - j   (binomial-free part).
"""
import math
from job1_tie_census import partitions, M_vector, v2, chi_b_vector, mn_character

def Dvec(lam, m):
    chi = chi_b_vector(lam, m)
    D = []
    for j in range(m+1):
        s = 0
        for b in range(j+1):
            s += math.comb(j, b) * (-1)**b * chi[b]
        D.append(s)
    return D, chi

def show(lam, m):
    D, chi = Dvec(lam, m)
    Ms = M_vector(lam, m)
    # sanity: D_j = 2^j M_j
    for j in range(m+1):
        assert D[j] == (1<<j)*Ms[j], (lam,j,D[j],Ms[j])
    print(f"lam={lam} m={m}")
    print(f"  chi_b = {chi}")
    print(f"  D_j   = {D}")
    print("  j : v2(D) | floor(j/2) | v2(D)-floor(j/2) | v2(Cmj) | val0=2v2D-j | val")
    for j in range(m+1):
        if D[j]==0:
            print(f"  {j:2d}:  inf"); continue
        vD=v2(D[j]); fl=j//2; vb=v2(math.comb(m,j))
        val0=2*vD-j; val=2*vb+val0
        print(f"  {j:2d} : {vD:3d}   | {fl:3d}      | {vD-fl:3d}            | {vb:3d}    | {val0:3d}      | {val}")
    print()

if __name__=="__main__":
    for lam,m in [((4,1,1),3),((4,2,2),4),((6,3,1,1,1),6),((8,2,2),6),
                  ((6,6),6),((4,4,2,2),6),((3,3,1,1),4),((2,2),2)]:
        show(lam,m)
