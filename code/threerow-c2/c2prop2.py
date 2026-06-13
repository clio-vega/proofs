import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb

def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n):
    return bin(n).count('1')
def vC(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j) + j*(j-1)*(j-2)*(j-3)

def val_direct(a,b,j):
    m=(a+b+2)//2
    Mjv=Mj((a,b,2),j,m)
    if Mjv==0: return None
    return j+2*(vC(m,j)+v2(Mjv))

def Delta_formula(a,b,j):
    return j - 2*s2(j) + 2*vC(a+3,j) + 2*vC(b+2,j) + 2*(v2(Q(a,b,j))-v2(Q(a,b,0)))

bad=0; tested=0
for A in range(3,30):
    for B in range(2,A+1):
        if (A+B)%2!=0: continue
        for J in range(0,B+1):  # interior
            vd=val_direct(A,B,J); v0=val_direct(A,B,0)
            if vd is None or v0 is None: continue
            df=Delta_formula(A,B,J)
            tested+=1
            if (vd-v0)!=df:
                bad+=1
                if bad<10: print("BAD",(A,B),J,"direct",vd-v0,"formula",df)
print(f"Prop2(c=2) verify: tested={tested} bad={bad}")
