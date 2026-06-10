"""
HYPOTHESIS H1 (b-picture):  the Newton-polygon minimum locus of
   T_lambda = sum_k C(m,k) chi_k (-i)^k        (G=0 iff T=0)
is EXACTLY the set of binary submasks of m with chi_k != 0:
   B* = { k : k AND m == k,  chi_k != 0 },
and beta_k = v2(C(m,k) chi_k) is CONSTANT = e on that set, strictly larger off it.

Test over ALL shapes (not just ties), all m<=12.  Also:
  - is v2(chi_k) constant on submasks?
  - distribution of |B*|.
"""
import math
from job1_tie_census import partitions, v2, chi_b_vector

def submasks(m):
    """all k with k&m==k, i.e. k submask of m, in increasing order."""
    out=[]; k=m
    # iterate submasks
    sub=m
    while True:
        out.append(sub)
        if sub==0: break
        sub=(sub-1)&m
    return sorted(out)

def test():
    tot=0
    H1_ok=0
    v2chi_const_on_sub=0
    Bsizes={}
    fails=[]
    for m in range(1,13):
        subs=set(submasks(m))
        for lam in partitions(2*m):
            chi=chi_b_vector(lam,m)
            beta={}
            for k in range(m+1):
                if chi[k]!=0:
                    beta[k]=v2(math.comb(m,k)*chi[k])
            if not beta: continue
            bmin=min(beta.values())
            B=set(k for k in beta if beta[k]==bmin)
            tot+=1
            Bsizes[len(B)]=Bsizes.get(len(B),0)+1
            # H1: B == submasks-with-chi-nonzero
            pred=set(k for k in subs if k<=m and chi[k]!=0)
            if B==pred: H1_ok+=1
            else: fails.append((m,lam,sorted(B),sorted(pred)))
            # is v2(chi_k) constant over submasks (chi!=0)?
            vs=set(v2(chi[k]) for k in subs if chi[k]!=0)
            if len(vs)<=1: v2chi_const_on_sub+=1
    print(f"shapes tested (m<=12): {tot}")
    print(f"  H1 (B* == submasks-of-m with chi!=0): {H1_ok}/{tot}")
    print(f"  v2(chi_k) constant over submasks:     {v2chi_const_on_sub}/{tot}")
    print(f"  |B*| distribution: {dict(sorted(Bsizes.items()))}")
    if fails:
        print(f"  H1 FAILURES ({len(fails)}); first 20:")
        for f in fails[:20]: print("   ",f)

if __name__=="__main__":
    test()
