from mn import Mj, fdim
from math import comb
def v2(n):
    n=abs(int(n)); 
    if n==0: return None
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(int(n)).count('1')

# verify clean Delta(b+1) formula
bad=0;tested=0
for m in range(2,60):
    for b in range(1,m):
        a=2*m-1-b
        if a<b: continue
        lam=(a,b,1)
        j=b+1
        if j>m: continue
        valb1 = j+2*v2(comb(m,j)*Mj(lam,j,m))
        val0 = 2*v2(fdim(lam))
        D_direct = valb1-val0
        D_form = (b-3)+2*v2(a+2)+2*(s2(m-b)-s2(a))
        if D_direct!=D_form:
            print(f"  MISMATCH lam={lam} m={m}: direct={D_direct} form={D_form}"); bad+=1
        tested+=1
print(f"Delta(b+1) clean formula: tested={tested} bad={bad}")

# val profile for a sample
def valj(lam,j,m):
    if j>m: return None
    M=Mj(lam,j,m)
    if M==0: return None
    return j+2*v2(comb(m,j)*M)
for lam,m in [((30,15,1),23),((31,14,1),23),((25,12,3),20),((28,13,3),22)]:
    prof=[(j,valj(lam,j,m)) for j in range(0,lam[1]+lam[2]+1)]
    V=min(v for j,v in prof if v is not None)
    print(f"\nlam={lam} m={m} b={lam[1]} V={V}")
    print("  j:val ", " ".join(f"{j}:{v}{'*' if v==V else ''}" for j,v in prof))
