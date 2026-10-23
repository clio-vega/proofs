import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from itertools import combinations
from psi import psi_alt, psi_headsum, mvec, B as Bfun
from sampler import clean_w
import random
exec(open('t4_P4_blocks.py').read().split("rng = random.Random")[0].split('"""',2)[2])

# ---------- NEGATIVE CONTROLS on P4 ----------
rng = random.Random(90210)
cA1=cA2=cB1=cB2=0; nA1=nA2=nB1=nB2=0
for trial in range(60):
    r = rng.randint(2,3); n = rng.randint(r+1, r+2)
    wq = clean_w(n, rng, lo=2, hi=12, dlo=1, dhi=5)
    theta = [rng.randint(-5,5) for _ in range(r)]
    m = mvec(theta)
    if len(set(m)) != r: continue
    jstar = max(range(r), key=lambda j: m[j])
    wsym = [Q(x) for x in wq[:-1]] + [W]; Mp = wq[:-1]
    A = block(wsym, theta, True); dA,cAv = lead(A)
    Bb = block(wsym, theta, False); dB,cBv = lead(Bb)
    # A1: use a NON-argmax index jbad
    others=[j for j in range(r) if j!=jstar]
    jbad = max(others, key=lambda j: m[j])
    hb=[m[l] for l in range(r) if l!=jbad]; htb=[hb[l]-l for l in range(r-1)]
    pA1=(-1)**(r+jbad+1)*Q(psi_alt(Mp,htb)); pdA1=m[jbad]-r+1
    nA1+=1
    if dA==pdA1 and sp.simplify(cAv-pA1)==0: cA1+=1
    # A2: drop the (-1)^{r+j*} sign
    ht=[m[l] for l in range(r) if l!=jstar]; htt=[ht[l]-l for l in range(r-1)]
    pA2=Q(psi_alt(Mp,htt)); nA2+=1
    if dA==m[jstar]-r+1 and sp.simplify(cAv-pA2)==0: cA2+=1
    # B1: theta instead of theta+1
    pB1=(-1)**r*Q(psi_alt(Mp,list(theta))); nB1+=1
    if dB==-r and sp.simplify(cBv-pB1)==0: cB1+=1
    # B2: exponent -(r-1) instead of -r
    nB2+=1
    if dB==-(r-1): cB2+=1
print("NEG CONTROL A1 (non-argmax j):      hits %d/%d" % (cA1,nA1))
print("NEG CONTROL A2 (sign dropped):      hits %d/%d" % (cA2,nA2))
print("NEG CONTROL B1 (theta not theta+1): hits %d/%d" % (cB1,nB1))
print("NEG CONTROL B2 (exponent -(r-1)):   hits %d/%d" % (cB2,nB2))

# ---------- P2: does B(T) take both signs? ----------
rng = random.Random(4242)
both=0; trials=0; example=None
for trial in range(400):
    r=2; n=rng.randint(3,4)
    wq = clean_w(n, rng, lo=1, hi=12, dlo=1, dhi=5)
    signs=set()
    vals={}
    for T in combinations(range(n), r):
        b = Bfun(wq, T); vals[T]=b
        if b>0: signs.add(1)
        elif b<0: signs.add(-1)
    trials+=1
    if len(signs)>1:
        both+=1
        if example is None: example=(wq,dict(vals))
print("P2: r=2 instances where B(T) takes BOTH signs: %d/%d" % (both,trials))
if example: print("   witness w=%s  B=%s" % ([str(x) for x in example[0]], {k:str(v) for k,v in example[1].items()}))

# ---------- P3: reciprocal pair structure ----------
rng = random.Random(31337)
res=[]
for r in (1,2,3):
    for n in (r, r+1, r+2):
        if n < 2: continue
        allzero_cnt=0; tot=0; nonzero_seen=0
        for trial in range(40):
            base = clean_w(n-1, rng, lo=2, hi=12, dlo=1, dhi=5)
            w = list(base) + [1/base[0]]       # force reciprocal pair w_1 * w_n = 1
            if len(set(w)) != n: continue
            theta=[rng.randint(-4,4) for _ in range(r)]
            if len(set(mvec(theta)))!=r: continue
            tot+=1
            v = psi_headsum(w, theta)
            if v==0: allzero_cnt+=1
            else: nonzero_seen+=1
        res.append((r,n,tot,allzero_cnt,nonzero_seen))
print("P3: M containing a reciprocal pair w_a*w_b=1, distinct m:")
for r,n,tot,z,nz in res:
    print("   r=%d n=%d : %d tested, Psi=0 in %d, Psi!=0 in %d   [predicted all-zero iff n<=r+1]" % (r,n,tot,z,nz))
