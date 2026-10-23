import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
from itertools import combinations
import psi as P
from psi import psi_headsum, psi_alt, mvec
import random
random.seed(55900)
def rw(n): return [F(random.randint(2,40), random.randint(1,7)) for _ in range(n)]
for trial in range(400):
    r = random.randint(1,3); n = random.randint(r, r+2); w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
for trial in range(60):
    r = random.randint(2,3); n = random.randint(r, r+2); w = rw(n)
    if len(set(w)) < n: continue
    theta = [random.randint(-4,4) for _ in range(r)]
    if (r,n,theta)==(2,3,[-1,3]):
        print("w =", w, " as floats", [float(x) for x in w])
        print("psi_headsum =", psi_headsum(w,theta), " psi_alt =", psi_alt(w,theta))
        print("e =", w[0]*w[1]*w[2])
        for i in range(3):
            print("  w_%d = %s   1-w^-2 = %s" % (i,w[i], 1-1/w[i]**2))
        for T in combinations(range(3),2):
            print("  T=%s B=%s alt=%s" % (T, P.B(w,T), P.alternant_ratio([w[i] for i in T], mvec(theta))))
        break
