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
        N=2*(m-j); return C(N,b-j-1)*(a-b+1)*(b*(a+1)-j*(j-1))//((b-j)*(b-j+1))
    if j==b: return 2*b*(m-b)
    if j==b+1: return b
    return 0
def valf(a,b,m,j):
    M=Mfac(a,b,m,j); cm=comb(m,j) if 0<=j<=m else 0
    if cm*M==0: return None
    return j+2*v2(cm*M)
# b=1 (a even): claim val(2)-val(0)=2 v2(m); val(1)-val(0)=3+2v2(m)+... >0
ok1=True
for m in range(2,200):
    a=2*m-2; b=1
    if a<b: continue
    if a%2: continue
    d20=valf(a,b,m,2)-valf(a,b,m,0)
    if d20!=2*v2(m): ok1=False
print("b=1: val(2)-val(0)==2 v2(m) :",ok1)
# b=4 (a odd): claim val(5)-val(3)=2 v2(m-3); val(4)-val(3)=1+2 v2(m-3)
ok4=True
for m in range(5,200):
    a=2*m-5; b=4
    if a<b or a%2==0: continue
    d53=valf(a,b,m,5)-valf(a,b,m,3)
    d43=valf(a,b,m,4)-valf(a,b,m,3)
    if d53!=2*v2(m-3) or d43!=1+2*v2(m-3): ok4=False
print("b=4: val(5)-val(3)==2 v2(m-3) and val(4)-val(3)==1+2 v2(m-3) :",ok4)
# spurious same-parity boundary competitor j=b+1: val(b+1)>val(j0) for b not in {1,2,4}
ok_sp=True; mn=99
for m in range(2,81):
    n=2*m
    for a in range(1,n+1):
        b=n-a-1
        if b<1 or b>a: continue
        if b in (1,2,4): continue
        j0=0 if a%2==0 else 3
        vb1=valf(a,b,m,b+1); v0=valf(a,b,m,j0)
        if vb1 is not None and v0 is not None:
            mn=min(mn,vb1-v0)
            if vb1<=v0: ok_sp=False
print("spurious j=b+1 ruled out (val(b+1)>val(j0)) for b not in {1,2,4}, m<=80:",ok_sp,"min margin",mn)
