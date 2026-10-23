"""P4 refined: test the TWO BLOCKS of the w_n -> infinity expansion separately.

Block A = sum over T containing the last index n.   Predicted leading order
          m_{j*} - r + 1, leading coeff (-1)^{r+j*} Psi_{(1^{r-1}),hat-theta}(M').
Block B = sum over T avoiding n.                   Predicted leading order exactly -r,
          leading coeff (-1)^r Psi_{(1^r),theta+1}(M').
These are the two statements the proof actually uses.  Signs 1-indexed in j*.
"""
import sys; sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from itertools import combinations
from psi import psi_alt, mvec, B as Bfun, alternant_ratio
from sampler import clean_w
import random
W = sp.Symbol('W')

def Bsym(wsym, T):
    n = len(wsym); x = [wsym[i] for i in T]; out = sp.Integer(1)
    for xi in x: out *= (1 - 1/xi**2)
    for a in range(len(x)):
        for b in range(a+1, len(x)): out *= (1 - 1/(x[a]*x[b]))
    for t in range(n):
        if t in T: continue
        for xi in x: out *= (1 - 1/(wsym[t]*xi)) / (1 - wsym[t]/xi)
    return out

def alt_sym(x, m):
    r = len(m); M = sp.Matrix(r, r, lambda i,j: x[i]**m[j])
    den = sp.Integer(1)
    for i in range(r):
        for k in range(i+1, r): den *= (x[k]-x[i])
    return sp.cancel(M.det()/den)

def block(wsym, theta, contains_last):
    r = len(theta); n = len(wsym); m = mvec(theta); tot = sp.Integer(0)
    for T in combinations(range(n), r):
        if ((n-1) in T) != contains_last: continue
        tot += Bsym(wsym, T) * alt_sym([wsym[i] for i in T], m)
    return sp.cancel(sp.together(tot))

def lead(expr):
    num, den = sp.fraction(sp.cancel(expr))
    pn = sp.Poly(sp.expand(num), W); pd = sp.Poly(sp.expand(den), W)
    return pn.degree()-pd.degree(), sp.Rational(pn.LC(), pd.LC())

def Q(v): return sp.Rational(v.numerator, v.denominator)

rng = random.Random(777001)
okA=badA=okB=badB=0; fA=[]; fB=[]; tested=0
for trial in range(150):
    r = rng.randint(1,3); n = rng.randint(r+1, r+2)      # need n>r so block B is nonempty
    wq = clean_w(n, rng, lo=2, hi=12, dlo=1, dhi=5)
    theta = [rng.randint(-5,5) for _ in range(r)]
    m = mvec(theta)
    if len(set(m)) != r: continue
    tested += 1
    jstar = max(range(r), key=lambda j: m[j])
    hat_m = [m[l] for l in range(r) if l != jstar]
    hat_theta = [hat_m[l]-l for l in range(r-1)]
    wsym = [Q(x) for x in wq[:-1]] + [W]
    Mp = wq[:-1]
    # Block A
    A = block(wsym, theta, True)
    dA, cA = lead(A)
    pdA = m[jstar]-r+1
    vA = psi_alt(Mp, hat_theta) if r-1 > 0 else F(1)
    pcA = (-1)**(r+jstar+1) * Q(vA)
    if dA == pdA and sp.simplify(cA-pcA)==0: okA+=1
    else:
        badA+=1
        if len(fA)<5: fA.append((r,n,theta,m,jstar+1,dA,cA,pdA,pcA))
    # Block B
    Bb = block(wsym, theta, False)
    dB, cB = lead(Bb)
    pdB = -r
    pcB = (-1)**r * Q(psi_alt(Mp, [t+1 for t in theta]))
    if dB == pdB and sp.simplify(cB-pcB)==0: okB+=1
    else:
        badB+=1
        if len(fB)<5: fB.append((r,n,theta,m,dB,cB,pdB,pcB))
print("tested %d (distinct m, n>r)" % tested)
print("P4 BLOCK A (T contains n):  OK=%d BAD=%d" % (okA,badA))
for f in fA: print("   FAIL", f)
print("P4 BLOCK B (T avoids n):    OK=%d BAD=%d" % (okB,badB))
for f in fB: print("   FAIL", f)
