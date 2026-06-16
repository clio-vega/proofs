from mn import Mj
from math import comb
def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def s2(n):
    return bin(n).count('1')
def carries(x,y):  # = v2 C(x+y,x)
    return v2(comb(x+y,x))
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

bad=[]
for m in range(8,40):
    for b in range(3,m):
        a=2*m-3-b
        if a<b: continue
        c=3; P=m-b-3
        if P<0: continue
        lam=(a,b,c)
        v0=val(lam,0,m)
        if v0 is None: continue
        N2=a*(b+3)-(b**2+2*b+3)
        N1=a*b**2+5*a*b+6*a-b**3-b**2-4*b-6
        # direct formula for Delta(b+2):
        G2 = (b+2) - 2*s2(b+3) + 2*carries(2*(P+1),b+3) + 2*v2(N2) - 2*v2(a-1) - 2*v2(P+2)
        # Delta(b+1):
        G1 = (b-1) - 2*s2(b+1) + 2*carries(2*(P+2),b+1) + 2*v2(N1) - 2*v2(a-1)
        # Delta(b+3) general:
        G3 = (b+5) - 2*s2(b+5) + 2*carries(2*P,b+5) - 2*v2(a-1) - 2*v2(P+2)
        for i,Gi in [(1,G1),(2,G2),(3,G3)]:
            j=b+i
            if j>m: continue
            vj=val(lam,j,m)
            if vj is None:
                # M=0 boundary; skip
                continue
            D=vj-v0
            if D!=Gi:
                bad.append((lam,m,i,D,Gi))
print("mismatches:",bad[:20],"count",len(bad))
