"""Rank of M for the EXACT configuration used in the written proof.
Independent engine (mod-p linear algebra) for the same claim that
construction2.py establishes combinatorially."""
import sys, random
sys.path.insert(0, '/home/clio/projects/code/q460')
from capacity_modp import (Delta, ordered_matchings, rank_modp, projective_rows, P, inv)

def cols(a):
    E = [2*i for i in range(1, 2*a+1)]
    O = [2*i-1 for i in range(1, 2*a+2)] + [4*a+2]
    return E, O

def run(a, q, label):
    random.seed(11)
    E, O = cols(a)
    Ps, Qs = ordered_matchings(E), ordered_matchings(O)
    n = a + (a+1)
    xs = []
    while len(xs) < n:
        v = random.randrange(2, P)
        if all(v != w and v != q*w % P and w != q*v % P for w in xs): xs.append(v)
    M = []
    for Pm in Ps:
        row = []
        for Qm in Qs:
            val = 1
            for i in range(a):
                for j in range(a+1):
                    val = val * Delta(Pm[i], Qm[j], xs[i], xs[a+j], q) % P
            row.append(val)
        M.append(row)
    r, pr = rank_modp(M), projective_rows(M)
    import math
    full = math.factorial(2*a)//(2**a)
    print("  a=%d N=%d M=%d  matrix %4dx%-5d  rank=%-4d proj-rows=%-4d  (2a)!/2^a=%-4d  %s"
          % (a, 2*a+1, 4*a+2, len(Ps), len(Qs), r, pr, full,
             "FULL" if (r == full and pr == full) else "SHORT"))

q6 = None
for c in range(2, 10**7):
    if (c*c - c + 1) % P == 0: q6 = c; break
if q6 is None:
    rt = pow(P-3, (P+1)//4, P); q6 = (1+rt)*inv(2) % P
assert pow(q6,3,P) == P-1

print("PROOF CONFIGURATION: prefix on even cols, suffix on odd cols + 4a+2")
print("generic q:")
for a in (1,2,3): run(a, 123456789, "generic")
print("q^3 = -1:")
for a in (1,2,3): run(a, q6, "root of unity")
