"""T6.  (a) diagnose the 1/72 half-hit of the determinant ansatz.
       (b) divisibility census: which reciprocal factors divide numerator(Psi)?
           Theorem R predicts  (w_a w_b - 1) | num  <=>  n <= r+1
                               (w_a^2   - 1) | num  <=>  n <= r
       (c) n = r: Psi = B(M) * s_lambda(M)  -> the ALTERNANT mechanism, separate from B.
       (d) EXPLICIT RATIONAL WITNESS: M over QQ, n >= r+2, NO reciprocal pair at all,
           m_j pairwise distinct, Psi^C = 0, every B(T) != 0.  (brief's P1 falsifier)
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B
import sympy as sp

# ---------- (a) diagnose the determinant half-hit ----------
print("(a) the 1/72 determinant hit, diagnosed")
random.seed(777)
VALS = [F(2), F(3), F(5), F(7), F(-2), F(-3), F(4), F(-5), F(3,2), F(-5,2), F(7,3)]
found = []
for n in range(2, 6):
    for trial in range(3):
        w = random.sample(VALS, n)
        for th in [(0,0),(0,1),(1,0),(1,1),(0,2),(2,0),(1,2),(2,1)]:
            m = [th[0], th[1]+1]
            if m[0] == m[1]: continue
            try: lhs = psi_headsum(w, list(th))
            except ZeroDivisionError: continue
            m1, m2 = m
            d = psi_headsum(w,[m1])*psi_headsum(w,[m2+1]) - psi_headsum(w,[m2])*psi_headsum(w,[m1+1])
            if d == lhs: found.append((n, th, m, w, lhs))
for f in found:
    n, th, m, w, lhs = f
    print(f"    n={n} theta={th} m={m}  Psi={lhs}   <- degenerate?  n<r+? ; |M|={n}, r=2")
if not found: print("    (none reproduced)")
print()

# ---------- (b) divisibility census ----------
print("(b) divisibility census of numerator(Psi) by reciprocal factors")
print(f"    {'n':>2} {'r':>2} {'theta':>10} | {'(w1w2-1)|num':>13} {'pred n<=r+1':>11} | {'(w1^2-1)|num':>13} {'pred n<=r':>9}")
def psi_sym(w, theta):
    n = len(w); r = len(theta); m = [theta[j]+j for j in range(r)]
    tot = sp.Integer(0)
    for T in itertools.combinations(range(n), r):
        x = [w[i] for i in T]
        Bt = sp.Integer(1)
        for xi in x: Bt *= (1 - 1/xi**2)
        for a2 in range(r):
            for b2 in range(a2+1, r): Bt *= (1 - 1/(x[a2]*x[b2]))
        for t in range(n):
            if t in T: continue
            for xi in x: Bt *= (1 - 1/(w[t]*xi))/(1 - w[t]/xi)
        Mx = sp.Matrix(r, r, lambda i, j: x[i]**m[j])
        V = sp.Integer(1)
        for i in range(r):
            for k in range(i+1, r): V *= (x[k]-x[i])
        tot += Bt*Mx.det()/V
    return sp.cancel(sp.together(tot))
cb_ok = cb_bad = 0
for n in range(1, 5):
    for r in range(1, min(n,3)+1):
        theta = tuple(range(r))
        w = sp.symbols(f'w1:{n+1}', nonzero=True)
        num, den = sp.fraction(psi_sym(list(w), list(theta)))
        num = sp.expand(num)
        d1 = (sp.rem(sp.Poly(num, w), sp.Poly(w[0]*w[1]-1, w)) == 0) if n >= 2 else None
        d2 = (sp.rem(sp.Poly(num, w), sp.Poly(w[0]**2-1, w)) == 0)
        p1, p2 = (n <= r+1), (n <= r)
        if n >= 2:
            cb_ok += (d1 == p1); cb_bad += (d1 != p1)
        cb_ok += (d2 == p2); cb_bad += (d2 != p2)
        print(f"    {n:>2} {r:>2} {str(theta):>10} | {str(d1):>13} {str(p1):>11} | {str(d2):>13} {str(p2):>9}")
print(f"    divisibility prediction from Theorem R: {cb_ok}/{cb_ok+cb_bad}")
print()

# ---------- (c) n = r : alternant mechanism ----------
print("(c) n=r: Psi = B(M)*s_lambda(M).  A Schur zero with every B(T) != 0 is a")
print("    vanishing mechanism that is NEITHER structural NOR a cancellation.")
# r=n=2, m=(0,3): det/V = x1^2+x1x2+x2^2 ; zero at x2/x1 a primitive cube root -> not rational.
# use m=(0,2): det/V = x1+x2 -> zero at w2=-w1, RATIONAL.
w = [F(3), F(-3)]
print(f"    r=n=2, theta=(0,1) so m=(0,2), det/V = x1+x2.  M={w}")
print(f"      B({{1,2}}) = {B(w,(0,1))}   (nonzero)")
print(f"      Psi(headsum)   = {psi_headsum(w,[0,1])}")
print(f"      Psi(alternant) = {psi_alt(w,[0,1])}")
print(f"      reciprocal pairs: {[(a,b) for a in range(2) for b in range(a,2) if w[a]*w[b]==1]}")
print()

# ---------- (d) explicit rational witness, n >= r+2, reciprocal-FREE ----------
print("(d) EXPLICIT RATIONAL WITNESS  (brief's pre-committed P1 falsifier)")
print("    r=1, n=3, m=1:  Psi = A_1/e,  A_1 = e_3 h_1 - e_2.  Search w3 in QQ.")
hits = []
for a_ in [F(x) for x in range(-9,10) if x != 0] + [F(3,2),F(5,2),F(-3,2),F(2,3),F(5,3),F(-5,2),F(7,2),F(-7,2),F(1,3),F(-1,3),F(4,3),F(-2,3)]:
    for b_ in [F(x) for x in range(-9,10) if x != 0] + [F(3,2),F(5,2),F(-3,2),F(2,3),F(5,3),F(-5,2),F(7,2),F(-7,2),F(1,3),F(-1,3),F(4,3),F(-2,3)]:
        if a_ == b_: continue
        # A_1 = e3*h1 - e2 = a b c (a+b+c) - (ab+ac+bc) = 0  -> solve quadratic in c
        # (ab)c^2 + (ab(a+b) - a - b) c - ab = 0
        A, Bq, C = a_*b_, a_*b_*(a_+b_) - a_ - b_, -a_*b_
        if A == 0: continue
        disc = Bq*Bq - 4*A*C
        if disc < 0: continue
        num_, den_ = disc.numerator, disc.denominator
        import math
        rn, rd = math.isqrt(num_), math.isqrt(den_)
        if rn*rn != num_ or rd*rd != den_: continue
        sq = F(rn, rd)
        for c_ in [(-Bq+sq)/(2*A), (-Bq-sq)/(2*A)]:
            if c_ == 0: continue
            w3 = [a_, b_, c_]
            if len(set(w3)) != 3: continue
            # reciprocal-free?
            if any(w3[i]*w3[j] == 1 for i in range(3) for j in range(i,3)): continue
            try: v1, v2 = psi_headsum(w3,[1]), psi_alt(w3,[1])
            except ZeroDivisionError: continue
            if v1 == 0 == v2:
                Bs = [B(w3,T) for T in itertools.combinations(range(3),1)]
                if all(x != 0 for x in Bs): hits.append((w3, Bs))
print(f"    rational reciprocal-free zeros found: {len(hits)}")
for (w3, Bs) in hits[:6]:
    print(f"      M = {w3}")
    print(f"        Psi(headsum)={psi_headsum(w3,[1])}  Psi(alternant)={psi_alt(w3,[1])}")
    print(f"        B(T) for the 3 singletons = {Bs}   (all nonzero)")
    print(f"        n=3 = r+2,  reciprocal support S = {{}} (empty),  n-|S|=3 >= r=1")
