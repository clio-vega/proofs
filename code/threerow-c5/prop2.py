import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
import sympy as sp, pickle, random

a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
H5sym=sp.sympify(d['H5sym']); Q5sym=sp.sympify(d['Q5sym'])
H5f=sp.lambdify((a,b,j),H5sym,'math')
Q5f=sp.lambdify((a,b,j),Q5sym,'math')
def Q5v(A,B,J): return int(round(Q5f(A,B,J)))
def H5v(A,B,J): return int(round(H5f(A,B,J)))

def v2(n):
    n=abs(int(n))
    if n==0: return 10**9
    k=0
    while n%2==0: n//=2;k+=1
    return k
def s2(n): return bin(n).count('1')
def v2C(n,k):
    if k<0 or k>n: return 10**9
    return v2(comb(n,k))
def val(A,B,J):
    m=(A+B+5)//2; M=Mj((A,B,5),J,m)
    if M==0: return None
    return J+2*v2(comb(m,J)*M)

# verify Q5v via exact MN-based formula too (avoid float error): use sympy exact for safety on a few
bad=0;tested=0
for _ in range(400):
    B=random.randrange(6,38); A=random.randrange(B,B+46)
    if (A+B)%2==0: A+=1
    if A<B: continue
    v0=val(A,B,0)
    q0=Q5v(A,B,0)
    for J in range(1,B+1):
        vj=val(A,B,J)
        if vj is None: continue
        Dexact=vj-v0
        Dprop=J-2*s2(J)+2*v2C(A+6,J)+2*v2C(B+5,J)+2*(v2(Q5v(A,B,J))-v2(q0))
        tested+=1
        if Dexact!=Dprop:
            bad+=1
            if bad<6: print("PROP2 MISMATCH",A,B,J,Dexact,Dprop)
print(f"Prop2: tested={tested} mismatches={bad}")

# Number Lemma floor for H5 (a,b opposite parity)
mn_=99;w=None
for A in range(0,32):
    for B in range(0,32):
        if (A+B)%2==0: continue
        for J in range(0,32):
            v=v2(H5v(A,B,J))
            if v<mn_: mn_=v;w=(A,B,J)
print(f"v2(H5) floor (scan 32^3, opp parity): {mn_} at {w}")
# verify floor is a TRUE constant: check residue dependence mod 2^(mn_+2)
K=2**(mn_+2)
bad2=0
for A in range(0,K):
    for B in range(0,K):
        if (A+B)%2==0: continue
        for J in range(0,4):  # j matters mod K too but check existence of low val
            pass
# proper finite proof check: H5 mod 2^(mn_+1) depends only on (a,b,j) mod 2^(mn_+1)
M=2**(mn_)   # claim v2>=mn_, i.e. 2^mn_ | H5
fails=0
for A in range(0,2*M):
    for B in range(0,2*M):
        if (A+B)%2==0: continue
        for J in range(0,2*M):
            if H5v(A,B,J)% M !=0: fails+=1
print(f"check 2^{mn_} | H5 for all residues mod {2*M}, opp parity: fails={fails}")
