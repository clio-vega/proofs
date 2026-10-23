import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from psi import psi_headsum, mvec
from sampler import clean_w
import random
rng = random.Random(5591)
for trial in range(300):
    r = rng.randint(1,3); n = rng.randint(r, r+2); w = clean_w(n, rng)
    theta = [rng.randint(-4,4) for _ in range(r)]
hitsA=[]
for trial in range(60):
    r = rng.randint(2,3); n = rng.randint(r, r+2); w = clean_w(n, rng); theta=[rng.randint(-4,4) for _ in range(r)]
    thetav=[-theta[j]-2*(j+1)+2 for j in range(r)]
    e=F(1)
    for x in w: e*=x
    lhs = psi_headsum([1/x for x in w], theta); rhs=(-1)**r*e**2*psi_headsum(w,thetav)
    if lhs==rhs: hitsA.append((r,n,theta,mvec(theta),e,lhs,rhs,psi_headsum(w,thetav)))
hitsB=[]
for trial in range(60):
    r = rng.randint(2,3); n = rng.randint(r, r+2); w = clean_w(n, rng); theta=[rng.randint(-4,4) for _ in range(r)]
    e=F(1)
    for x in w: e*=x
    lhs = psi_headsum([1/x for x in w], theta); rhs=(-1)**r*e**(2*r)*psi_headsum(w,theta)
    if lhs==rhs: hitsB.append((r,n,theta,mvec(theta),e,lhs,psi_headsum(w,theta)))
print("CONTROL A hits (%d):" % len(hitsA))
for h in hitsA: print("   r=%d n=%d theta=%s m=%s e=%s  lhs=%s  Psi(thetav)=%s  m_repeat=%s  e2r_eq_e2=%s"
    % (h[0],h[1],h[2],h[3],h[4],h[5],h[7], len(set(h[3]))<h[0], h[4]**(2*h[0])==h[4]**2))
print("CONTROL B hits (%d):" % len(hitsB))
for h in hitsB: print("   r=%d n=%d theta=%s m=%s e=%s  lhs=%s  Psi(theta)=%s  m_repeat=%s"
    % (h[0],h[1],h[2],h[3],h[4],h[5],h[6], len(set(h[3]))<h[0]))
