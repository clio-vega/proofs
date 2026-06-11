"""
Sharp Newton polygon of G_lambda(i) in the M_j coordinate.

  G = sum_j C(m,j) i^j pi^j M_j,  pi=1+i,  M_j = <s_lam, p1^{2(m-j)} e2^j> >= 0.
  val(j) = j + 2 v2(C(m,j) M_j).
  V = min_j val(j),  J* = {j : val(j)=V}.

Leading pi-coefficient:  write C(m,j)M_j = 2^{a_j} u_j (u_j odd), a_j=v2(C(m,j)M_j).
  Each j in J* contributes u_j i^{j-a_j} to coeff of pi^V.
  C = sum_{j in J*} u_j i^{j-a_j}.   Claim: C = |J*| (mod pi), so C==0 mod pi <=> |J*| even.

We compute, per shape: M_j, val(j), J*, |J*|, C (Gaussian int), v_pi(C), and v_pi(G).
Census the distribution of |J*|.
"""
import sys
from sympy import Integer, I, expand, im, re
import symfunc as sf
from collections import Counter

def v2(n):
    n = int(n)
    if n == 0: return None
    k = 0
    while n % 2 == 0: n//=2; k+=1
    return k

def vpi(z):
    # v_pi of a Gaussian integer z, pi=1+i.  v_pi(z) = v2(N(z)) where N=norm.
    z = complex(z)
    a = int(round(z.real)); b = int(round(z.imag))
    if a==0 and b==0: return None
    N = a*a + b*b
    return v2(N)

def binom(n,k):
    from math import comb
    return comb(n,k)

def analyze(lam, m):
    slam = sf.schur_p(lam)
    p1sq = sf.power(sf.h1(),2)
    e2 = sf.e_n(2)
    M=[]
    for j in range(m+1):
        f = sf.mul(sf.power(p1sq, m-j), sf.power(e2, j))
        M.append(int(sf.inner(slam,f)))
    # val(j)
    val={}
    for j in range(m+1):
        if M[j]==0: continue
        val[j] = j + 2*v2(binom(m,j)*M[j])
    V = min(val.values())
    Jstar = sorted([j for j in val if val[j]==V])
    # leading coefficient C = sum_{j in J*} u_j i^{j-a_j}
    C = Integer(0)
    for j in Jstar:
        coef = binom(m,j)*M[j]
        a = v2(coef); u = coef // (2**a)
        C += u * (I**((j - a) % 4))
    C = expand(C)
    return dict(lam=lam,m=m,M=M,val=val,V=V,Jstar=Jstar,absJ=len(Jstar),
                C=C, vpiC=vpi(C))

def main():
    mmax = int(sys.argv[1]) if len(sys.argv)>1 else 7
    cnt = Counter()
    odd_ge3 = []
    leadcancel_vs_J = Counter()  # (|J*| parity, leading cancels?)
    parity_violation=[]
    for m in range(1,mmax+1):
        for lam in sf.partitions(2*m):
            d=analyze(lam,m)
            cnt[d['absJ']]+=1
            if d['absJ']%2==1 and d['absJ']>=3:
                odd_ge3.append((lam,m,d['absJ']))
            # check C = |J*| mod pi : C even-norm <=> |J*| even
            cancels = (d['vpiC'] is None) or (d['vpiC']>=1)
            if cancels != (d['absJ']%2==0):
                parity_violation.append((lam,m,d['absJ'],d['C']))
    print(f"|J*| distribution (m<={mmax}):", dict(sorted(cnt.items())))
    print(f"shapes with |J*| odd and >=3: {len(odd_ge3)}")
    for x in odd_ge3[:20]: print("   ODD>=3:",x)
    print(f"violations of 'C cancels mod pi <=> |J*| even': {len(parity_violation)}")
    for x in parity_violation[:10]: print("   VIOL:",x)

if __name__=="__main__":
    main()
