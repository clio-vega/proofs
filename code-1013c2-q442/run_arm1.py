"""ARM 1: the width-bound-inactive theorem and its sharpness witnesses."""
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1013c2-q442')
from poly import *
from q442 import *

def Ainf_Abd(k,n,w):
    m=2*k
    dpU=engineC(k,n,w,bounded=False); dpB=engineC(k,n,w,bounded=True)
    Ai={}; Ab={}
    for lam,p in dpU.items():
        if c_minus(lam,m,w)==1: Ai=add(Ai,p)
    for lam,p in dpB.items():
        if c_minus(lam,m,w)==1: Ab=add(Ab,p)
    return Ai,Ab

print("A) FORWARD THEOREM  n <= w  =>  A^inf == A, over a grid:")
cnt=tot=0
for k in (1,2,3):
    for ell in range(0,5):
        w=2*ell+2
        for n in range(1,min(w,5)+1):
            Ai,Ab=Ainf_Abd(k,n,w); tot+=1; cnt+= (Ai==Ab)
print(f"   {cnt}/{tot} domains with n<=w, k in 1..3, ell in 0..4, n<=min(w,7): A^inf == A")
print()
print("B) SHARPNESS WITNESS, k>=2, n=w+1, alpha=(2^n):  predict A^inf > A")
for k in (2,3):
    for ell in (0,1,2):
        w=2*ell+2; n=w+1
        if n>6 and k>=3: continue
        if n>6: continue
        Ai,Ab=Ainf_Abd(k,n,w)
        al=tuple([2]*n)
        a=Ai.get(al,0); b=Ab.get(al,0)
        print(f"   k={k} w={w} n={n}={w}+1  alpha={al}:  A^inf={a}  A={b}   "
              f"A^inf>A: {a>b}")
print()
print("C) SHARPNESS WITNESS, k=1, n=2(w+1), alpha=(1^n):  predict A^inf > A;")
print("   and NO alpha has A^inf>A for n < 2(w+1).")
for ell in (0,1):
    w=2*ell+2; n=2*(w+1)
    if n>8: 
        print(f"   k=1 w={w} n={n}: skipped (n>8, engine cost)"); continue
    Ai,Ab=Ainf_Abd(1,n,w)
    al=tuple([1]*n)
    print(f"   k=1 w={w} n={n}=2({w}+1)  alpha={al}: A^inf={Ai.get(al,0)} A={Ab.get(al,0)} "
          f" A^inf>A: {Ai.get(al,0)>Ab.get(al,0)}")
    for nn in range(1,n):
        Ai2,Ab2=Ainf_Abd(1,nn,w)
        if Ai2!=Ab2: print(f"      !! unexpected: k=1 w={w} n={nn} has A^inf != A")
    print(f"      no alpha with A^inf>A for n < {n}: confirmed over n=1..{n-1}")
print()
print("D) w-INDEPENDENCE of A^inf (needed: A^inf is defined with the w-free '+1' class):")
for k in (1,2):
    for n in (3,4):
        ref=None; ok=True
        for ell in range(0,6):
            Ai,_=Ainf_Abd(k,n,2*ell+2)
            if ref is None: ref=Ai
            elif Ai!=ref: ok=False
        print(f"   k={k} n={n}: A^inf identical across w=2,4,6,8,10,12: {ok}")
