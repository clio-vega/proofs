import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from itertools import combinations
import psi as P
from psi import psi_headsum, mvec
import random
random.seed(55900)
def rw(n): return [F(random.randint(2,40), random.randint(1,7)) for _ in range(n)]
# replay EXACTLY the control loop so the two hits are the same draws
for trial in range(400):
    r = random.randint(1,3); n = random.randint(r, r+2); w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
hits=[]
for trial in range(60):
    r = random.randint(2,3); n = random.randint(r, r+2); w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
    tot = F(0)
    for T in combinations(range(n), r):
        x=[w[i] for i in T]; tot += P.B(w,T)*P.alternant_ratio(x, list(theta))
    true = psi_headsum(w, theta)
    if tot == true:
        hits.append((r,n,theta,mvec(theta),tot,true))
for h in hits:
    r,n,theta,m,tot,true = h
    print("HIT r=%d n=%d theta=%s  m=%s  unshifted=%s  value=%s  distinct_m=%s distinct_theta=%s"
          % (r,n,theta,m,theta,tot, len(set(m))==r, len(set(theta))==r))
