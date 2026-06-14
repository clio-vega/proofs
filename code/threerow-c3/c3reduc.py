from math import comb
def Hf(av,bv,jv):
    return ( (av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4) )
def Qf(av,bv,jv):
    return (av-1)*(bv-2)*Hf(av,bv,jv)-720*comb(jv,6)
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

# a EVEN reduction target (L3''): after using j*C(a+4,j)=(a+4)C(a+3,j-1) to absorb v2(j),
# T3+v2(j) >= v2C(b+3,j)+v2Q3(j)-v2(b+3).  Claim (L3''): v2C(b+3,j)+v2Q3(j) >= v2(b+3)+1.
# But by symmetry could also reduce via b+3. Test BOTH and the symmetric min.
fA=fB=0; tightA=0
for mv in range(4,80):
    n=2*mv
    for av in range(2,n+1,2):
        for bv in range(3,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(1,bv+1):
                qj=Qf(av,bv,jv)
                LA=v2C(bv+3,jv)+v2(qj)-v2(bv+3)   # reduce via a+4
                LB=v2C(av+4,jv)+v2(qj)-v2(av+4)   # reduce via b+3
                if LA<1: fA+=1
                if LB<1: fB+=1
                if LA==1: tightA+=1
print(f"a EVEN (L3'') reduce-via-(a+4): v2C(b+3,j)+v2Q3(j)>=v2(b+3)+1 fails={fA} (tight {tightA})")
print(f"a EVEN reduce-via-(b+3): v2C(a+4,j)+v2Q3(j)>=v2(a+4)+1 fails={fB}")

# Understand v2 H(j) for a even (heavy factor odd so v2Q3(j)=v2H(j) for j<=5; tip for j>=6)
# Claim: for a even, v2 H(j) >= 1 for all j>=1 ?  and what is exact min?
import collections
mh=collections.defaultdict(lambda:10**9)
for mv in range(4,60):
    n=2*mv
    for av in range(2,n+1,2):
        for bv in range(3,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(1,min(bv,5)+1):
                mh[v2(jv)]=min(mh[v2(jv)],v2(Hf(av,bv,jv)))
print("a EVEN min v2 H(j) by v2(j) (j<=5):",dict(mh))

# a ODD: relative bound. tildeDelta(j)=j+3-2s2(j)+2*U(j), U=v2C(a+4,j)+v2C(b+3,j)+v2Q3(j)-v2Q3(0)
# Reduction: a+4,b+3 BOTH ODD now. v2Q3(0)=v2(a-1)+v2(a+3)+v2(b-2)+v2(b+2).
# Heavy (a-1)(b-2) even. Test: U(j) >= floor giving tildeDelta>=0 with ties {5,7,9} for (1,2).
# Try: is v2Q3(j) related to v2((a-1)(b-2)) + something?
fodd=0
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            for jv in range(4,bv+1):
                td=jv+3-2*s2(jv)+2*(v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0)))
                if td<0: fodd+=1
print(f"a ODD tildeDelta(j)>=0 for 4<=j<=b fails={fodd} (m<=79)")
