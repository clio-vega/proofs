"""T4.  THEOREM R (reduction), the conjecture that reorganises Q565.

    (a)  w_a^2 = 1       =>  Psi_{(1^r),theta}(M) = Psi_{(1^r),theta}(M \ {w_a})
    (b)  w_a w_b = 1, a!=b =>  Psi_{(1^r),theta}(M) = Psi_{(1^r),theta}(M \ {w_a,w_b})

If true, the reciprocal locus is not a ZERO locus -- it is a RESTRICTION locus; Psi forgets
the pair.  Vanishing is then the degenerate case |M \ S| < r (empty sum).

Tested with BOTH independent implementations on the left and on the right.
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt

random.seed(1024)
VALS = [F(2), F(3), F(5), F(7), F(-2), F(-3), F(4), F(-5), F(11), F(-4), F(3,2), F(5,3), F(-7,2)]

def ok_set(w):
    return len(set(w)) == len(w) and all(x != 0 for x in w)

okA = badA = okB = badB = 0
failsA, failsB = [], []
thetalist = []
for r in range(1, 4):
    for base in ([tuple(range(r)), tuple(2*j for j in range(r)), tuple(1-j for j in range(r)),
                  tuple(3*j for j in range(r)), tuple(-j for j in range(r))]):
        m = [base[j]+j for j in range(r)]
        if len(set(m)) == r: thetalist.append((r, list(base)))

for (r, theta) in thetalist:
    for n in range(1, 7):
        for trial in range(6):
            # ---- (b) one reciprocal pair at positions 0,1 ----
            if n >= 2:
                v = random.choice([F(2), F(3), F(5), F(-2), F(-3), F(5,2)])
                rest = random.sample(VALS, n-2)
                w = [v, 1/v] + rest
                if ok_set(w):
                    wred = w[2:]
                    try:
                        L1, L2 = psi_headsum(w, theta), psi_alt(w, theta)
                        R1 = psi_headsum(wred, theta) if len(wred) >= 0 else F(0)
                        R2 = psi_alt(wred, theta)
                    except ZeroDivisionError:
                        pass
                    else:
                        if L1 == L2 == R1 == R2: okB += 1
                        else:
                            badB += 1
                            if len(failsB) < 4: failsB.append((n, r, theta, w, L1, L2, R1, R2))
            # ---- (a) w_0 = +-1 ----
            s = random.choice([F(1), F(-1)])
            rest = random.sample(VALS, n-1)
            w = [s] + rest
            if ok_set(w):
                wred = w[1:]
                try:
                    L1, L2 = psi_headsum(w, theta), psi_alt(w, theta)
                    R1, R2 = psi_headsum(wred, theta), psi_alt(wred, theta)
                except ZeroDivisionError:
                    pass
                else:
                    if L1 == L2 == R1 == R2: okA += 1
                    else:
                        badA += 1
                        if len(failsA) < 4: failsA.append((n, r, theta, w, L1, L2, R1, R2))

print(f"Theorem R(a)  w_a=+-1  reduction M -> M\\{{w_a}}      : {okA}/{okA+badA}")
for f in failsA: print("   FAIL", f)
print(f"Theorem R(b)  w_a w_b=1 reduction M -> M\\{{w_a,w_b}} : {okB}/{okB+badB}")
for f in failsB: print("   FAIL", f)

# ---- ITERATED: chain of pairs, the case prop:recip cannot see ----
print()
print("ITERATED reduction (chains).  n -> n-2mu-nu ; vanishing iff that is < r")
hits = 0; tot = 0
for r in range(1, 4):
    theta = list(range(r))
    for (mu, nu, extra) in [(1,0,0),(2,0,0),(2,0,1),(3,0,0),(1,1,0),(2,1,0),(1,1,1),(1,0,2),(2,1,1),(3,0,1)]:
        base = [F(2), F(3), F(5), F(7)][:mu]
        w = []
        for v in base: w += [v, 1/v]
        if nu >= 1: w.append(F(1))
        if nu >= 2: w.append(F(-1))
        w += [F(4), F(9), F(25)][:extra]
        n = len(w)
        if not ok_set(w): continue
        wred = [F(4), F(9), F(25)][:extra]
        try:
            L = psi_headsum(w, theta); Lr = psi_alt(w, theta)
            R = psi_alt(wred, theta) if len(wred) >= r else F(0)
        except ZeroDivisionError:
            continue
        tot += 1
        agree = (L == Lr == R)
        if agree: hits += 1
        print(f"  r={r} n={n} mu={mu} nu={nu} free={extra}: n-|S|={extra} "
              f"{'<' if extra<r else '>='} r  Psi={L}  reduced={R}  "
              f"{'ok' if agree else 'MISMATCH'}")
print(f"\niterated reduction: {hits}/{tot}")
