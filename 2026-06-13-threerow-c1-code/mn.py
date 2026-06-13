from math import comb
from functools import lru_cache

# Murnaghan-Nakayama via beta-numbers (abacus). chi^lambda(mu).
# Represent partition as tuple sorted desc (no trailing zeros).

def beta_set(lam):
    k = len(lam)
    return [lam[i] + (k-1-i) for i in range(k)]  # distinct, desc

def from_beta(betas):
    betas = sorted(betas, reverse=True)
    k = len(betas)
    parts = [betas[i] - (k-1-i) for i in range(k)]
    return tuple(x for x in parts if x>0)

@lru_cache(maxsize=None)
def chi(lam, mu):
    # lam, mu tuples desc. mu = cycle type.
    if sum(lam)==0:
        return 1
    r = mu[0]
    rest = mu[1:]
    k = len(lam)
    betas = beta_set(lam)
    bset = set(betas)
    total = 0
    for i,bi in enumerate(betas):
        nb = bi - r
        if nb < 0 or nb in bset:
            continue
        # height = number of betas strictly between nb and bi
        height = sum(1 for b in betas if nb < b < bi)
        newbetas = [b for b in betas if b!=bi] + [nb]
        newlam = from_beta(newbetas)
        total += (-1)**height * chi(newlam, rest)
    return total

def fdim(mu):
    mu = tuple(x for x in mu if x>0)
    if not mu: return 1
    n=sum(mu)
    return chi(mu, (1,)*n)

def Mj(lam, j, m):
    lam=tuple(lam)
    s=0
    from fractions import Fraction
    acc=Fraction(0)
    for t in range(j+1):
        mu = (2,)*t + (1,)*(2*m-2*t)
        acc += comb(j,t)*(-1)**t * chi(lam, mu)
    val = acc / (2**j)
    assert val.denominator==1
    return int(val)

if __name__=="__main__":
    # quick sanity: f^lambda
    print("f(2,1)=",fdim((2,1)), "expect 2")
    print("f(3,2,1)=",fdim((3,2,1)),"expect 16")
    print("M_0 (3,2,1) m=3 =",Mj((3,2,1),0,3),"expect f=16")

def test_threerow():
    bad=0;tested=0; mism=[]
    for m in range(2,10):
        n=2*m
        for a in range(1,n+1):
            for b in range(1,a+1):
                c=n-a-b
                if c<1 or c>b: continue
                lam=(a,b,c)
                for j in range(0,c+1):
                    M=Mj(lam,j,m)
                    conj=fdim((a-j,b-j,c-j))
                    tested+=1
                    if M!=conj:
                        bad+=1
                        if len(mism)<40: mism.append((lam,m,j,M,conj))
                jhi=c+1
                if jhi<=m:
                    if Mj(lam,jhi,m)!=0:
                        mism.append(('ABOVE',lam,jhi,Mj(lam,jhi,m)))
    for x in mism: print("  ",x)
    print(f"tested={tested} mismatches={bad}")

if __name__=="__main__":
    test_threerow()
