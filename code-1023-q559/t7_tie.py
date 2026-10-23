import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from itertools import combinations, product
from psi import psi_alt, psi_headsum, mvec
from sampler import clean_w
import random
exec(open('t4_P4_blocks.py').read().split("rng = random.Random")[0].split('"""',2)[2])

rng=random.Random(606)
# A) is max(m) = -1 a genuine ORDER TIE between the two blocks?
tie=0; ties_cancel=0; examples=[]
tested=0
for trial in range(2000):
    r=rng.randint(2,3); n=rng.randint(r+1,r+2)
    theta=[rng.randint(-8,2) for _ in range(r)]
    m=mvec(theta)
    if len(set(m))!=r or max(m)!=-1: continue
    tested+=1
    if tested>25: break
    wq=clean_w(n,rng,lo=2,hi=12,dlo=1,dhi=5)
    wsym=[Q(x) for x in wq[:-1]]+[W]; Mp=wq[:-1]
    dA,cA=lead(block(wsym,theta,True)); dB,cB=lead(block(wsym,theta,False))
    assert dA==dB==-r, (dA,dB,-r)
    tie+=1
    s=sp.simplify(cA+cB)
    if s==0:
        ties_cancel+=1
        if len(examples)<5: examples.append((r,n,theta,m,[str(x) for x in wq]))
print("max(m) = -1 instances: %d, all with dA == dB == -r : %s" % (tie, tie==tested or tie>0))
print("   of these, cA + cB == 0 (leading term CANCELS):  %d/%d" % (ties_cancel, tie))
for e in examples: print("     CANCELLING:", e)

# B) how often does the sign (-1)^{r+j*} equal +1 ?  (diagnoses NEG CONTROL A2's 36/52)
rng2=random.Random(90210)
plus=minus=0
for trial in range(60):
    r=rng2.randint(2,3); n=rng2.randint(r+1,r+2)
    wq=clean_w(n,rng2,lo=2,hi=12,dlo=1,dhi=5)
    theta=[rng2.randint(-5,5) for _ in range(r)]
    m=mvec(theta)
    if len(set(m))!=r: continue
    jstar=max(range(r),key=lambda j:m[j])
    if (-1)**(r+jstar+1)==1: plus+=1
    else: minus+=1
print("A2 diagnosis: sign (-1)^{r+j*} = +1 in %d cases, = -1 in %d cases (total %d)" % (plus,minus,plus+minus))

# C) boundary: |M| < r  =>  Psi = 0 even with distinct m
bz=bn=0
for r in (2,3,4):
    for n in range(1,r):
        for theta in product(range(-2,3), repeat=r):
            m=mvec(list(theta))
            if len(set(m))!=r: continue
            w=clean_w(n,random.Random(n*100+r),lo=2,hi=9,dlo=1,dhi=3)
            bn+=1
            if psi_headsum(w,list(theta))==0: bz+=1
print("BOUNDARY |M| < r, distinct m: Psi = 0 in %d/%d" % (bz,bn))
