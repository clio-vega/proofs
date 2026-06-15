"""End-to-end verification of the COMPLETE Gap-2 proof chain (a odd, b even, a>=b>=4, a+b=2m-3).
Links:
 (LemA)  2 s2(k) <= k+1                                    all k>=0
 (LemC)  2 s2(j)+2 v2(j(j-1)) <= j+3                       j>=2
 (heart) 2 v2C(j,6)+2 s2(j) <= j-1                         j>=6   [<=> 2 s2(j-6)<=(j-6)+1]
 (NL1)   v2C(F+1,j)+v2(j(j-1)) >= v2(F)   F even,2<=j<=F+1  [subset id, self-contained]
 (i)     v2C(b+3,j) >= v2(b+2)+s2(j)-(j+3)/2               4<=j<=b   [NL1+LemC]
 (B)     v2C(a+4,j) >= v2(a-1)+v2(a+3)+s2(j)-(j+3)/2-1     6<=j<=a   [a-absorb+heart]
 (star2) tip: 2 v2 P2 >= 2 vP0+2s2(j)-(j+3)               6<=j<=b   [b-absorb+B]
 (star1) heavy: v2[C(a+4,j)C(b+3,j)H] >= vH0+g            4<=j<=b   [Pi1 trivial + Pi2 via a-absorb+(i)]
 (MAIN)  tildeDelta(j) >= 0                                4<=j<=b
"""
from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3)-6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def Qf(av,bv,jv): return (av-1)*(bv-2)*Hf(av,bv,jv)-720*comb(jv,6)

f=dict.fromkeys(["LemA","LemC","heart","NL1","i","B","star2","star1","MAIN"],0)
# standalone lemmas
for k in range(0,30000):
    if 2*s2(k)>k+1: f["LemA"]+=1
for j in range(2,30000):
    if 2*s2(j)+2*(v2(j)+v2(j-1))>j+3: f["LemC"]+=1
for j in range(6,30000):
    if 2*v2C(j,6)+2*s2(j)>j-1: f["heart"]+=1
for F in range(2,400,2):
    for j in range(2,F+2):
        if v2C(F+1,j)+v2(j)+v2(j-1) < v2(F): f["NL1"]+=1

# family lemmas + MAIN
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            a1,a3,b2 = v2(av-1),v2(av+3),v2(bv+2)
            vP0=v2(Qf(av,bv,0)); vH0=v2(Hf(av,bv,0))
            assert vP0==a1+a3+b2+v2(bv-2), (av,bv)   # sanity vP0
            assert vH0==a3+b2, (av,bv)
            for jv in range(4,bv+1):
                G2=2*s2(jv)-(jv+3)
                if 2*v2C(bv+3,jv) < 2*b2+G2: f["i"]+=1
                if jv>=6 and 2*v2C(av+4,jv) < 2*(a1+a3)+G2-2: f["B"]+=1
                if jv>=6:
                    P2=720*comb(av+4,jv)*comb(bv+3,jv)*comb(jv,6)
                    if 2*v2(P2) < 2*vP0+G2: f["star2"]+=1
                Hprod=comb(av+4,jv)*comb(bv+3,jv)*Hf(av,bv,jv)
                if 2*v2(Hprod) < 2*vH0+G2: f["star1"]+=1
                # MAIN tildeDelta = (j+3) -2s2 +2U,  U=v2P - vP0
                U=v2(comb(av+4,jv)*comb(bv+3,jv)*Qf(av,bv,jv))-vP0
                td=jv+3-2*s2(jv)+2*U
                if td<0: f["MAIN"]+=1
for k,v in f.items(): print(f"{v:6d} fails  {k}")
