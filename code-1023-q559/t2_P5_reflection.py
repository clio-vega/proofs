import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, mvec
from sampler import clean_w
import random
rng = random.Random(5591)
ok=bad=0; fails=[]
for trial in range(300):
    r = rng.randint(1,3); n = rng.randint(r, r+2)
    w = clean_w(n, rng)
    theta = [rng.randint(-4,4) for _ in range(r)]
    m = mvec(theta)
    thetav = [-theta[j] - 2*(j+1) + 2 for j in range(r)]   # 1-indexed j
    assert mvec(thetav) == [-x for x in m], (mvec(thetav), m)
    winv = [1/x for x in w]
    e = F(1)
    for x in w: e *= x
    lhs = psi_headsum(winv, theta)
    rhs = (-1)**r * e**(2*r) * psi_headsum(w, thetav)
    if lhs == rhs: ok+=1
    else:
        bad+=1
        if len(fails)<3: fails.append((r,n,theta,lhs,rhs))
print("P5 reflection identity: OK=%d BAD=%d" % (ok,bad))
for f in fails: print("  FAIL", f)

# NEGATIVE CONTROL A: wrong exponent on e (e^2 instead of e^{2r})
okA=badA=0
for trial in range(60):
    r = rng.randint(2,3); n = rng.randint(r, r+2)
    w = clean_w(n, rng); theta=[rng.randint(-4,4) for _ in range(r)]
    thetav=[-theta[j]-2*(j+1)+2 for j in range(r)]
    e=F(1)
    for x in w: e*=x
    lhs = psi_headsum([1/x for x in w], theta)
    rhs = (-1)**r * e**2 * psi_headsum(w, thetav)
    if lhs==rhs: badA+=1
    else: okA+=1
print("NEG CONTROL A (e^2 not e^{2r}): disagree=%d agree(bad)=%d" % (okA,badA))

# NEGATIVE CONTROL B: m^vee = -m replaced by m^vee = m (no reflection of exponents)
okB=badB=0
for trial in range(60):
    r = rng.randint(2,3); n = rng.randint(r, r+2)
    w = clean_w(n, rng); theta=[rng.randint(-4,4) for _ in range(r)]
    e=F(1)
    for x in w: e*=x
    lhs = psi_headsum([1/x for x in w], theta)
    rhs = (-1)**r * e**(2*r) * psi_headsum(w, theta)
    if lhs==rhs: badB+=1
    else: okB+=1
print("NEG CONTROL B (theta^vee = theta): disagree=%d agree(bad)=%d" % (okB,badB))
