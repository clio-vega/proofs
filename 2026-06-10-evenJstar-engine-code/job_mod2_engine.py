"""
THE MOD-2 ENGINE.
b-picture: Phi(z)=sum_k C(m,k) chi_k z^k = <s_lam,(p1^2+z p2)^m>.
Claim:  Phi(z) ≡ chi^lam(2^m) (1+z)^m  (mod 2).   [since p1^2≡p2 mod 2]
=> if chi^lam(2^m) odd, the min-coeff-valuation locus B* = submasks(m).

Test 1: verify Phi(z) ≡ chi(2^m)(1+z)^m mod 2 for all shapes.
Test 2 (general b-box): Phi(z)/2^e mod 2 == z^{j0}(1+z)^g ?  (e=min coeff val)
Test 3 (SHARP j-picture): within the parity class of J*, is the scaled even/odd
   sub-polynomial's bottom-mod-2 part also a shifted power of (1+y)?
"""
import math
from job1_tie_census import partitions, M_vector, v2, chi_b_vector

def poly_mod2(coeffs):
    """list of int coeffs -> set of exponents with odd coeff."""
    return frozenset(j for j,c in enumerate(coeffs) if c%2!=0)

def onepluszsup(m):
    """support of (1+z)^m mod 2 = submasks of m (Lucas)."""
    return frozenset(k for k in range(m+1) if (k & m)==k)

def is_shift_pow_1pz(support):
    """is `support` == j0 + submasks(g) for some j0,g?  (i.e. z^{j0}(1+z)^g mod2)"""
    if not support: return False,None,None
    S=sorted(support); j0=S[0]
    shifted=frozenset(s-j0 for s in S)
    g=max(shifted)
    return (shifted==onepluszsup(g)), j0, g

def main():
    # Test 1
    t1=0; t1ok=0
    for m in range(1,13):
        for lam in partitions(2*m):
            chi=chi_b_vector(lam,m)
            coeffs=[math.comb(m,k)*chi[k] for k in range(m+1)]
            lhs=poly_mod2(coeffs)
            chi2m=chi[m]   # chi^lam(2^m)
            rhs=onepluszsup(m) if chi2m%2!=0 else frozenset()
            t1+=1
            if lhs==rhs: t1ok+=1
    print(f"Test1  Phi ≡ chi(2^m)(1+z)^m mod2 : {t1ok}/{t1}")

    # Test 2: general b-box via Phi/2^e mod2 == shifted (1+z)^g
    t2=0; t2ok=0; bad=[]
    for m in range(1,13):
        for lam in partitions(2*m):
            chi=chi_b_vector(lam,m)
            coeffs=[math.comb(m,k)*chi[k] for k in range(m+1)]
            nz=[v2(c) for c in coeffs if c!=0]
            if not nz: continue
            e=min(nz)
            red=[ (c//(1<<e)) for c in coeffs ]   # divide by 2^e
            sup=poly_mod2(red)
            ok,j0,g=is_shift_pow_1pz(sup)
            t2+=1
            if ok: t2ok+=1
            else: bad.append((m,lam,sorted(sup)))
    print(f"Test2  Phi/2^e ≡ z^j0 (1+z)^g mod2 : {t2ok}/{t2}")
    for x in bad[:10]: print("   bad:",x)

if __name__=="__main__":
    main()
