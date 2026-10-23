"""Refutation arm: hunt for distinct m with Psi = 0 at a NON-degenerate w.
Pushes r to 3 and 4 and |M| past the 424-sample range of 10-22."""
import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from psi import psi_alt, psi_headsum, mvec
from sampler import clean_w
import random, itertools
rng = random.Random(1123581321)

# (a) exhaustive theta boxes at several r, n, 2 independent w per cell
cells=[]
for r,box,ns in [(2,range(-5,6),[2,3,4,5]), (3,range(-4,5),[3,4,5]), (4,range(-3,4),[4,5])]:
    for n in ns:
        cells.append((r,box,n))
tot=0; zeros=0; zero_examples=[]
sign_report={}
for r,box,n in cells:
    ws=[clean_w(n,rng,lo=2,hi=20,dlo=1,dhi=7) for _ in range(2)]
    cnt=0; z=0
    for theta in itertools.product(box, repeat=r):
        m=mvec(list(theta))
        if len(set(m))!=r: continue
        for w in ws:
            cnt+=1
            v=psi_alt(w,list(theta))
            if v==0:
                z+=1
                if len(zero_examples)<8: zero_examples.append((r,n,list(theta),m,[str(x) for x in w]))
    tot+=cnt; zeros+=z
    print("  r=%d n=%d box=[%d,%d]: %d evaluations (distinct m), Psi=0 in %d" % (r,n,min(box),max(box),cnt,z))
print("REFUTATION ARM TOTAL: %d evaluations, %d zeros" % (tot,zeros))
for e in zero_examples: print("   ZERO:", e)

# (b) POSITIVE CONTROL on the same instrument: repeated m must give 0 every time
pc=0; pcn=0
for r,box,n in cells:
    w=clean_w(n,rng,lo=2,hi=20,dlo=1,dhi=7)
    for theta in itertools.product(range(-3,4), repeat=r):
        m=mvec(list(theta))
        if len(set(m))==r: continue
        pcn+=1
        if psi_alt(w,list(theta))==0: pc+=1
print("POSITIVE CONTROL (repeated m => 0): %d/%d" % (pc,pcn))

# (c) cross-check psi_alt against the DEFINITION on a sample of the big sweep
rng2=random.Random(5)
agree=0; dis=0
for _ in range(120):
    r=rng2.randint(2,3); n=rng2.randint(r,r+2)
    w=clean_w(n,rng2,lo=2,hi=20,dlo=1,dhi=7); theta=[rng2.randint(-5,5) for _ in range(r)]
    if psi_alt(w,theta)==psi_headsum(w,theta): agree+=1
    else: dis+=1
print("definition-vs-alternant on sweep sampler: agree=%d disagree=%d" % (agree,dis))
