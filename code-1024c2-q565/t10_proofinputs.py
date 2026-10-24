"""T10. POSITIVE CONTROL ON THE PROOF ITSELF  (point-evaluation form, no expansion).

Over the common denominator D = prod_{a<b<=L}(x_a^2-x_b^2) * prod_{a<=L,p}(c_p^2-x_a^2),
with exact cofactors
    (x_a^2-x_b^2)/(e_a x_a + e_b x_b) = e_a (x_a - e_a e_b x_b)
    (c_p^2-x_a^2)/(e_a x_a + c_p)     = (c_p - e_a x_a)
the numerator is the polynomial
    P = sum_eps eta^nu * prod_{a<b<=L} e_a (x_a - e_a e_b x_b) * prod_{a,p} (c_p - e_a x_a).

P is evaluated EXACTLY at random integer points (Fractions, no floats, no expansion).
 (a) P == 0 identically  iff  L >= s+2       [the identity]        -- 12 random points
 (b) deg_{x_1} P <= N-1                      [proof Step 1]        -- Lagrange interpolation
                                                                      in x_1 at N+3 nodes
 (c) P = 0 on every hyperplane x_a = +- x_b   [proof Step 2]        -- 8 random points each
(b),(c) are run only where P != 0, so they are not vacuous.
"""
import itertools, random
from fractions import Fraction as F

def Pval(L, s, nu, xs, cs):
    tot = F(0)
    for eps in itertools.product([1,-1], repeat=L):
        eta = 1
        for e in eps: eta *= e
        t = F(eta) if (nu % 2) else F(1)
        for a in range(L):
            for b in range(a+1, L): t *= eps[a]*(xs[a] - eps[a]*eps[b]*xs[b])
        for a in range(L):
            for p in range(s): t *= (cs[p] - eps[a]*xs[a])
        tot += t
    return tot

def rnd(k, lo=3, hi=400): return [F(random.randint(lo, hi)) for _ in range(k)]

def deg_in_x1(L, s, nu, maxdeg):
    """exact degree of P in x_1 via finite differences on maxdeg+2 nodes"""
    cs = rnd(s); rest = rnd(L-1) if L >= 2 else []
    nodes = [F(i) for i in range(1, maxdeg+3)]
    vals = [Pval(L, s, nu, [t]+rest, cs) for t in nodes]
    # Newton forward differences at unit-spaced nodes: degree = last nonzero difference order
    d = vals[:]; deg = -1
    for order in range(len(vals)):
        if any(x != 0 for x in d): deg = order
        d = [d[i+1]-d[i] for i in range(len(d)-1)]
        if not d: break
    return deg

random.seed(20261024)
print(f"{'L':>2} {'s':>2} {'N':>2} {'d':>2} | {'L>=s+2':>6} {'P==0':>5} | {'deg_x1 P':>8} {'<=N-1':>6} | {'P=0 on x_a=+-x_b':>17} | deg_x1 Delta")
agree = dis = degbad = divbad = 0
for L in range(1, 6):
    for s in range(0, 5):
        if L+s > 7: continue
        for delta in (0, 1):
            N = L+s; nu = N-1+delta
            zero = all(Pval(L, s, nu, rnd(L), rnd(s)) == 0 for _ in range(12))
            pred = (L >= s+2)
            agree += (zero == pred); dis += (zero != pred)
            if zero:
                print(f"{L:>2} {s:>2} {N:>2} {delta:>2} | {str(pred):>6} {'yes':>5} | {'--':>8} {'--':>6} | {'--':>17} | {2*(L-1)}")
                continue
            dg = deg_in_x1(L, s, nu, 2*(L-1)+2*s)
            okdeg = dg <= N-1; degbad += (not okdeg)
            allvan = True
            for a in range(L):
                for b in range(a+1, L):
                    for sgn in (1,-1):
                        for _ in range(8):
                            xs = rnd(L); cs = rnd(s); xs[b] = sgn*xs[a]
                            if Pval(L, s, nu, xs, cs) != 0: allvan = False
            divbad += (not allvan)
            print(f"{L:>2} {s:>2} {N:>2} {delta:>2} | {str(pred):>6} {'no':>5} | {dg:>8} {str(okdeg):>6} | {str(allvan):>17} | {2*(L-1)}")
print()
print(f"(a) identity  L>=s+2 <=> P==0 : {agree} cells agree, {dis} disagree   (L<=5, s<=4, N<=7, delta in 0,1)")
print(f"(b) Step 1 degree-bound violations : {degbad}")
print(f"(c) Step 2 divisibility violations : {divbad}")
print()
print("--- THE SQUEEZE (deg_x1 Delta vs the Step-1 bound), cell by cell ---")
mis = 0
for L in range(1, 7):
    for s in range(0, 6):
        if L+s > 9: continue
        forced = 2*(L-1) > L+s-1
        if forced != (L >= s+2): mis += 1
print(f"    2(L-1) > L+s-1  <=>  L >= s+2  :  mismatches over L<=6,s<=5 : {mis}")
