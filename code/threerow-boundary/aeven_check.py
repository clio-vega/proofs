from mn import Mj
from math import comb
def v2(n):
    n=abs(n); k=0
    if n==0: return 99
    while n%2==0: n//=2; k+=1
    return k
def vprod(start,cnt):  # v2 of product start..start+cnt-1
    s=0
    for x in range(start,start+cnt): s+=v2(x)
    return s
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)
bad=[]
for m in range(8,42):
    for b in range(7,m):
        a=2*m-3-b
        if a<b or a%2!=0: continue   # a even
        if b%2!=1: continue
        c=3; P=m-b-3
        if P<0 or b<6: continue
        lam=(a,b,c); v0=val(lam,0,m)
        if v0 is None: continue
        N2=a*(b+3)-(b**2+2*b+3); N1=a*b**2+5*a*b+6*a-b**3-b**2-4*b-6
        assert v2(N2)==1, (lam,v2(N2))   # claim v2(N2)=1 for a even
        Q2=(b+3)//2; W2=vprod(P+2,Q2)    # (P+2)...(P+1+Q2)
        Q1=(b+1)//2; W1=vprod(P+3,Q1)    # (P+3)...(P+2+Q1)
        D2pred=1+2*(W2-v2(P+2))
        D1pred=2*(W1+v2(N1)-1)
        for i,pred in [(2,D2pred),(1,D1pred)]:
            j=b+i
            if j>m: continue
            vj=val(lam,j,m)
            if vj is None: continue
            if vj-v0!=pred: bad.append((lam,m,i,vj-v0,pred))
        # also assert W1>=3
        assert W1>=3,(lam,W1)
print("a-even mismatches:",bad[:20],"count",len(bad))
print("OK: v2(N2)=1 always, W1>=3 always (asserts passed)")
