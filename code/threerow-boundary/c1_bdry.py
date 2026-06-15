from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n))
    if n==0: return None
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n):
    return bin(n).count('1')

# my hand formula for v2(f^lambda), lambda=(a,b,1):
# v2(f) = 2 - s2(m) + s2(a) + s2(b-1) + v2(a-b+1) - v2(a+2) - v2(b+1)
# val(b+1) = (b+1) + 2 v2(b) + 2[ s2(b+1)+s2(m-b-1)-s2(m) ]
# Delta = val(b+1)-val(0)
bad=0; tested=0
for m in range(2,40):
    a_b1=2*m-1
    for b in range(1,(a_b1)//2+1):
        a=a_b1-b
        if a<b: continue
        lam=(a,b,1)
        # direct
        f=fdim(lam)
        v0=2*v2(f)
        # my formula for v2(f)
        vf_formula = 2 - s2(m) + s2(a) + s2(b-1) + v2(a-b+1) - v2(a+2) - v2(b+1)
        if vf_formula != v2(f):
            print(f"  v2(f) MISMATCH lam={lam} m={m}: formula={vf_formula} actual={v2(f)}")
            bad+=1
        # boundary val(b+1)
        j=b+1
        if j>m: continue
        Mb1=Mj(lam,j,m)
        if Mb1!=b: print(f"  M_b+1 != b: lam={lam} M={Mb1}")
        valb1_direct = j + 2*v2(comb(m,j)*Mb1)
        valb1_formula = (b+1) + 2*v2(b) + 2*(s2(b+1)+s2(m-b-1)-s2(m))
        if valb1_direct!=valb1_formula:
            print(f"  val(b+1) MISMATCH lam={lam} m={m}: formula={valb1_formula} direct={valb1_direct}")
            bad+=1
        tested+=1
print(f"tested={tested} bad={bad}")
