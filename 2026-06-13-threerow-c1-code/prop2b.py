from mn import Mj
from math import comb
def C(n,k):
    if k<0 or k>n or n<0: return 0
    return comb(n,k)
def v2(x):
    x=abs(x); c=0
    while x%2==0:x//=2;c+=1
    return c
def Mfac(a,b,m,j):
    if j<=b-1:
        N=2*(m-j); num=C(N,b-j-1)*(a-b+1)*(b*(a+1)-j*(j-1)); den=(b-j)*(b-j+1)
        assert num%den==0
        return num//den
    if j==b: return 2*b*(m-b)
    if j==b+1: return b
    return 0
bad=0
for m in range(2,13):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        for j in range(0,m+1):
            if Mj((a,b,1),j,m)!=Mfac(a,b,m,j): bad+=1
print("full-support factored M_j mismatches=",bad)

# Prop-2 val for j<=b-1 (the bulk); val(j)=j+2 v2(C(m,j) M_j)
def valdir(a,b,m,j):
    M=Mj((a,b,1),j,m); cm=comb(m,j)
    if cm*M==0: return None
    return j+2*v2(cm*M)
# formula for j<=b-1:
def valform(a,b,m,j):
    K=b*(a+1)
    t=v2(comb(m,j))+v2(C(2*(m-j),b-j-1))+v2(K-j*(j-1))-v2((b-j)*(b-j+1))+v2(a-b+1)
    return j+2*t
bad2=0; cnt=0
for m in range(2,13):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        for j in range(0,b):  # j<=b-1
            vd=valdir(a,b,m,j)
            vf=valform(a,b,m,j)
            cnt+=1
            if vd!=vf:
                bad2+=1
                if bad2<8: print("VALMM",(a,b),m,j,vd,vf)
print(f"Prop-2 c=1 val formula (j<=b-1): tested {cnt}, mismatches={bad2}")
