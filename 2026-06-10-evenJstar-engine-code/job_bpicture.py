"""
B-PICTURE.  G = 2^{-m}(1+i)^m * Phi(-i),  Phi(-i)=sum_k C(m,k) chi_k (-i)^k.
So G=0 iff T:=sum_k C(m,k) chi_k (-i)^k = 0.
Term k has v_pi = 2 v2(C(m,k) chi_k).  Define beta_k=v2(C(m,k)chi_k),
B* = argmin beta_k.  Question: is |B*| even when >1?  Is it a 2-adic box?
Compare with the j-picture J* and mu.
"""
import math
from job1_tie_census import partitions, M_vector, v2, chi_b_vector

def jpicture(lam,m):
    Ms=M_vector(lam,m)
    val={}
    for j in range(m+1):
        if Ms[j]==0: continue
        val[j]=j+2*v2(math.comb(m,j)*Ms[j])
    mu=min(val.values()); J=sorted(j for j in val if val[j]==mu)
    return J,mu,val

def bpicture(lam,m):
    chi=chi_b_vector(lam,m)
    beta={}
    for k in range(m+1):
        if chi[k]==0: continue
        beta[k]=v2(math.comb(m,k)*chi[k])
    bmin=min(beta.values()); B=sorted(k for k in beta if beta[k]==bmin)
    return B,bmin,beta,chi

def is_pow2box(S):
    """Is sorted set S an affine box j0+{subset sums of distinct powers of 2}?"""
    S=sorted(S); n=len(S)
    if n&(n-1)!=0: return False,None  # not power of 2
    j0=S[0]; diffs=sorted(set(s-j0 for s in S))
    # generators = elements that are powers of two present
    gens=[d for d in diffs if d>0 and (d&(d-1))==0]
    # build box
    from itertools import chain, combinations
    box=set()
    for r in range(len(gens)+1):
        for combo in combinations(gens,r):
            box.add(j0+sum(combo))
    if box==set(S): return True,gens
    return False,None

if __name__=="__main__":
    print("Comparison j-picture vs b-picture on ties:")
    tot=0; bpow2=0; beven=0; bgt1=0
    bad=[]
    for m in range(2,13):
        for lam in partitions(2*m):
            J,mu,val=jpicture(lam,m)
            if len(J)<2: continue   # only ties
            B,bmin,beta,chi=bpicture(lam,m)
            tot+=1
            ok,gens=is_pow2box(B)
            if len(B)>1:
                bgt1+=1
                if len(B)%2==0: beven+=1
                if ok: bpow2+=1
                else: bad.append((m,lam,B))
    print(f"ties: {tot}")
    print(f"  |B*|>1: {bgt1};  of these |B*| even: {beven};  pow2-box: {bpow2}")
    if bad:
        print("  non-pow2box B* (first 15):")
        for x in bad[:15]: print("   ",x)
    # detailed look
    print("\nDetailed:")
    for lam,m in [((4,1,1),3),((4,2,2),4),((6,3,1,1,1),6),((8,2,2),6),((6,6),6),((2,2),2),((3,3,1,1),4)]:
        J,mu,val=jpicture(lam,m)
        B,bmin,beta,chi=bpicture(lam,m)
        print(f"lam={lam} m={m}: J*={J}(mu={mu})  B*={B}(bmin={bmin})  |J*|={len(J)} |B*|={len(B)}")
        print(f"   chi={chi}")
        print(f"   beta_k={[beta.get(k,'inf') for k in range(m+1)]}")
