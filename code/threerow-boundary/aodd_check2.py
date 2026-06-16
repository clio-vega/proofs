from mn import Mj
from math import comb
def v2(n):
    n=abs(n); k=0
    if n==0: return 99
    while n%2==0: n//=2; k+=1
    return k
def vprod(start,cnt): return sum(v2(x) for x in range(start,start+cnt))
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)
bad=[]
# also directly check the FINAL bound inequalities used in the proof:
ineqfail=[]
for m in range(8,42):
    for b in range(10,m):
        if b%2!=0: continue
        a=2*m-3-b
        if a<b or a%2!=1: continue
        c=3; P=m-b-3
        if P<0: continue
        beta=b//2
        if beta<5: continue
        lam=(a,b,c); v0=val(lam,0,m)
        if v0 is None: continue
        N1=a*b**2+5*a*b+6*a-b**3-b**2-4*b-6
        N2=a*(b+3)-(b**2+2*b+3)
        U  = vprod(P+1, beta+2)
        Up = vprod(P+2, beta+1)
        Upp= vprod(P+3, beta)
        D3 = 2*(U - v2(P+beta+1) - v2(P+2)) - 3
        D2 = 2*(Up - 1 + v2(N2) - v2(a-1) - v2(P+2))
        D1 = 2*(Upp + v2(N1) - v2(a-1)) - 3
        for i,pred in [(1,D1),(2,D2),(3,D3)]:
            j=b+i
            if j>m: continue
            vj=val(lam,j,m)
            if vj is None: continue
            if vj-v0!=pred: bad.append((lam,m,i,vj-v0,pred))
        # the three KEY inequalities (two/single-factor Lemma F):
        if not (U - v2(P+beta+1) - v2(P+2) >= 1): ineqfail.append(("F2-top",lam))
        if not (Up - v2(P+beta+1) - v2(P+2) >= 1): ineqfail.append(("F2-sub2",lam))
        if not (Upp - v2(P+beta+1) >= 1):          ineqfail.append(("F1-sub1",lam))
        if not (v2(N1) >= 1):                       ineqfail.append(("N1even",lam))
print("a-odd form mismatches:", bad[:10], "count", len(bad))
print("Lemma-F inequality failures:", ineqfail[:10], "count", len(ineqfail))
