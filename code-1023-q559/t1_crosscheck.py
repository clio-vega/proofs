import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, mvec
import random
random.seed(55900)

def rw(n):
    return [F(random.randint(2,40), random.randint(1,7)) for _ in range(n)]

ok = bad = 0
for trial in range(400):
    r = random.randint(1,3)
    n = random.randint(r, r+2)
    w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
    a = psi_headsum(w, theta); b = psi_alt(w, theta)
    if a == b: ok += 1
    else:
        bad += 1
        if bad < 4: print("MISMATCH", r, n, theta, a, b)
print("thm:alt crosscheck  OK=%d BAD=%d" % (ok, bad))

# negative control: drop the +j shift (use m_j = theta_j), should DISAGREE
ok2 = bad2 = 0
for trial in range(60):
    r = random.randint(2,3); n = random.randint(r, r+2)
    w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
    import psi as P
    # monkey: alternant with unshifted m
    from itertools import combinations
    tot = F(0)
    for T in combinations(range(n), r):
        x=[w[i] for i in T]; tot += P.B(w,T)*P.alternant_ratio(x, list(theta))
    if tot == psi_headsum(w, theta): bad2 += 1
    else: ok2 += 1
print("NEG CONTROL no-shift alternant: disagree=%d agree(bad)=%d" % (ok2, bad2))
