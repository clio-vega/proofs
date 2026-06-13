from math import comb
def C(n,k):
    if k<0 or k>n or n<0: return 0
    return comb(n,k)
def v2(x):
    x=abs(x); c=0
    while x%2==0:x//=2;c+=1
    return c
def s2(n): return bin(n).count('1')
def Mfac(a,b,m,j):
    if j<=b-1:
        N=2*(m-j); return C(N,b-j-1)*(a-b+1)*(b*(a+1)-j*(j-1))//((b-j)*(b-j+1))
    if j==b: return 2*b*(m-b)
    if j==b+1: return b
    return 0
def valf(a,b,m,j):
    M=Mfac(a,b,m,j); cm=comb(m,j) if 0<=j<=m else 0
    if cm*M==0: return None
    return j+2*v2(cm*M)
def Dgen(a,b,j):
    return (j-4)-2*s2(j-2)+2*(v2(C(a+2,j))+v2(C(b-1,j-2))+v2(b+1)+v2(b*(a+1)-j*(j-1))-v2(a+1))
bad=0;cnt=0
for m in range(2,16):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        v0=valf(a,b,m,0)
        for j in range(2,b):
            vj=valf(a,b,m,j)
            if vj is None: continue
            cnt+=1
            if vj-v0!=Dgen(a,b,j): bad+=1
print(f"GENERAL c=1 Delta formula: tested {cnt}, mismatches={bad}")
# K-even: verify Delta(1)=-1, Delta(2)=-2, Delta(3)=-3 (forced descent), j0=3
desc_fail=0
for m in range(2,40):
    n=2*m
    for a in range(1,n+1,2):     # a odd
        b=n-a-1
        if b<1 or b>a or b%2==1: continue   # b even => K even
        v0=valf(a,b,m,0)
        for j in [1,2,3]:
            if j<=b-1 or True:
                vj=valf(a,b,m,j)
                if vj is None: continue
                if vj-v0 != -j: desc_fail+=1
print(f"K-even forced descent val(j)=val(0)-j for j=1,2,3: fails={desc_fail}")
