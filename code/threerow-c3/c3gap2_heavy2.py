from math import comb
def v2(n):
    if n==0:return 10**9
    n=abs(n);c=0
    while n%2==0:n//=2;c+=1
    return c
def s2(n):return bin(int(n)).count('1')
def Hf(av,bv,jv):
    return ((av+3)*(av+4)*(bv+2)*(bv+3) -6*(av+3)*(bv+2)*comb(jv,1)
        -6*(av*bv+av+2*bv)*comb(jv,2)+36*comb(jv,3)+72*comb(jv,4))
def Gf(av,bv,jv): return (av+4)*(bv+3)-6*jv
def Ef(av,bv,jv): return (av*bv+av+2*bv)*comb(jv,2)-6*comb(jv,3)-12*comb(jv,4)

# Tests
t1=t2=t3=tPi2=0
import collections
mPi2=collections.defaultdict(lambda:10**9)
for mv in range(4,80):
    n=2*mv
    for av in range(1,n+1,2):
        for bv in range(4,av+1,2):
            if n-av-bv!=3:continue
            vH0=v2(av+3)+v2(bv+2)
            for jv in range(4,bv+1):
                Ca,Cb=comb(av+4,jv),comb(bv+3,jv)
                # T1: v2[Ca Cb H] >= vH0  (trivial-target candidate)
                if v2(Ca*Cb*Hf(av,bv,jv)) < vH0: t1+=1
                # check H decomposition identity
                if Hf(av,bv,jv)!=(av+3)*(bv+2)*Gf(av,bv,jv)-6*Ef(av,bv,jv): t2+=1
                # Pi2 = 6 Ca Cb E ; target v2 Pi2 >= vH0 + g(j) ; use x2: 2v2Pi2 >= 2vH0+2s2(j)-(j+3)
                Pi2=6*Ca*Cb*Ef(av,bv,jv)
                R=2*vH0+2*s2(jv)-(jv+3)
                if 2*v2(Pi2) < R: tPi2+=1
                mPi2[v2(jv)]=min(mPi2[v2(jv)], 2*v2(Pi2)-R)
print("T1 v2[Ca Cb H] >= vH0  fails:",t1)
print("H = (a+3)(b+2)G - 6E identity fails:",t2)
print("Pi2 target 2v2(6 Ca Cb E) >= 2vH0+2s2-(j+3) fails:",tPi2)
print("Pi2 margin (x2) by v2(j):",dict(sorted(mPi2.items())))
