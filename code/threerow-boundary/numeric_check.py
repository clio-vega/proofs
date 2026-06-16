from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(n)
    if n==0: return 99
    k=0
    while n%2==0: n//=2; k+=1
    return k
def val(lam,j,m):
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)

# enumerate c=3 shapes, box-interior (b large), check (A),(B) suffices-claims and sharp Δ
viol_A=[]; viol_B=[]; rows=[]
for m in range(8,46):
    for b in range(1,m):
        a=2*m-3-b
        if a<b: continue
        c=3
        if b<c: continue
        P=m-b-3
        if P<0: continue
        lam=(a,b,c)
        # box interior condition from note: a even b>=6 (j0=0); a odd b>=9 (j0=3)
        j0 = 0 if a%2==0 else 3
        if a%2==0 and b<6: continue
        if a%2==1 and b<9: continue
        theta = j0
        v0=val(lam,0,m)
        if v0 is None: continue
        for i in (1,2,3):
            j=b+i
            if j>m: continue
            vj=val(lam,j,m)
            if vj is None: continue
            D=vj-v0
            if not (D> -theta):
                rows.append(("BOUND VIOL",lam,m,i,D,theta))
        # check A: v2(N2)>=v2(P+1)-1
        N2=a*(b+3)-(b**2+2*b+3)
        if v2(N2) < v2(P+1)-1:
            viol_A.append((lam,m,v2(N2),v2(P+1)))
        # check B candidate: v2(N1)>=v2(P+2)-1
        N1=a*b**2+5*a*b+6*a-b**3-b**2-4*b-6
        if v2(N1) < v2(P+2)-1:
            viol_B.append((lam,m,v2(N1),v2(P+2)))
print("bound violations:", rows[:10], "count", len(rows))
print("A viol (v2(N2)>=v2(P+1)-1):", viol_A[:10], "count", len(viol_A))
print("B viol (v2(N1)>=v2(P+2)-1):", viol_B[:10], "count", len(viol_B))
