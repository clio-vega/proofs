import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from math import factorial
import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5f=sp.lambdify((a,b,j),sp.sympify(d['H5sym']),'math')
def H5v(A,B,J): return int(round(H5f(A,B,J)))
def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def vff(x,r):
    s=0
    for t in range(r): s+=v2(x-t)
    return s
def vpsi(A,B,J): return vff(A+2,J-4)+vff(B+1,J-4)+v2(H5v(A,B,J))

# j=6 aeven, kappa=8. Find failing residues (a even,b odd) with vpsi<8, list a few with v2 breakdown
J=6;K=8;M=2**K;off=M
found=[]
for ra in range(0,M,2):
    A=off+ra
    for rb in range(1,M,2):
        B=off+rb
        vp=vpsi(A,B,J)
        if vp<8:
            found.append((ra,rb,vp, vff(A+2,2),vff(B+1,2),v2(H5v(A,B,J))))
print(f"j=6 aeven failing residues count (vpsi<8): {len(found)}; first 8:")
for x in found[:8]: print("  ra,rb,vpsi,vA,vB,vH:",x)

# Take first failing (ra,rb); test VALID shapes a>=b>=10 with a≡ra,b≡rb mod M
if found:
    ra,rb=found[0][0],found[0][1]
    print(f"\nTesting valid shapes with a≡{ra}, b≡{rb} mod {M}:")
    cnt=0
    for k in range(0,8):
        B=rb+ k*M   # b ≡ rb
        if B<10: continue
        for la in range(0,6):
            A=ra+la*M
            if A<B or A<10: continue
            print(f"  a={A} b={B} (a>=b? {A>=B}): vpsi={vpsi(A,B,J)} Delta(6)={6-2*2+2*(vpsi(A,B,J)-2*(6-2))}")
            cnt+=1
            if cnt>=8: break
        if cnt>=8: break
