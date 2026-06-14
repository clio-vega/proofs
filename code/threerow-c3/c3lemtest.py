from math import comb
def Qf(av,bv,jv):
    H=( (av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4) )
    return (av-1)*(bv-2)*H-720*comb(jv,6)
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def v2C(N,k):
    if k<0 or k>N:return 10**9
    return s2(k)+s2(N-k)-s2(N)

# ===== a EVEN candidate lemmas =====
# (E1) v2 Q3(j) >= v2(j)+1   (the floor for the trivial-binomial bound)
# (E2) T3(j) >= 1 - v2(j)
# (E3) v2 H(j) >= 1 for all j>=1 (a even); and =1 exactly determining ties
fail_E1=fail_E2=0; tightE2=0
for mv in range(4,80):
    n=2*mv
    for av in range(2,n+1,2):
        for bv in range(3,av+1,2):
            if n-av-bv!=3:continue
            q0=Qf(av,bv,0)
            for jv in range(1,bv+1):
                qj=Qf(av,bv,jv)
                T3=v2C(av+4,jv)+v2C(bv+3,jv)+v2(qj)-v2(q0)
                if v2(qj)<v2(jv)+1: fail_E1+=1
                if T3<1-v2(jv): fail_E2+=1
                if T3==1-v2(jv): tightE2+=1
print(f"a EVEN: (E1) v2Q3(j)>=v2(j)+1 fails={fail_E1}; (E2) T3>=1-v2(j) fails={fail_E2} (tight {tightE2})")

# ===== a ODD candidate lemmas =====
# forced descent: Delta(1)=-1,Delta(2)=-2,Delta(3)=-3 always?
# relative bound for j>=4: tildeDelta(j)=Delta(j)+3 >= 0
def Delta(av,bv,jv):
    return jv-2*s2(jv)+2*(v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0)))
fd_fail=0; rel_fail=0
# candidate odd Compensation: U(j):=v2C(a+4,j)+v2C(b+3,j)+v2Q3(j)-v2Q3(0). find bound
import collections
minU=collections.defaultdict(lambda:10**9)
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            if Delta(av,bv,1)!=-1 or Delta(av,bv,2)!=-2 or Delta(av,bv,3)!=-3: fd_fail+=1
            for jv in range(4,bv+1):
                if Delta(av,bv,jv)+3<0: rel_fail+=1
                U=v2C(av+4,jv)+v2C(bv+3,jv)+v2(Qf(av,bv,jv))-v2(Qf(av,bv,0))
                minU[v2(jv)]=min(minU[v2(jv)],U)
print(f"a ODD: forced-descent D1=-1,D2=-2,D3=-3 fails={fd_fail}; relative tildeDelta(j>=4)>=0 fails={rel_fail}")
print("a ODD: min U(j) by v2(j) (j>=4):", dict(minU))
# what is Delta(j)+3 floor? Delta(j)+3 = j+3-2s2(j)+2U. Check =j+3-2s2(j)+2*(floor)
