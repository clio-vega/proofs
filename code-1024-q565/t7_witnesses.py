"""T7.  The P1 falsifier, explicitly, over QQ -- and the cleanest possible form of it.

At r=1, theta=(0) the three-case closed form collapses:  m=0  =>  Psi = (e - 1/e)/e,
i.e.  Psi^C_{(1),(0)}(M) = 1 - e^{-2},  e = prod_t w_t,  for EVERY n.
So the zero locus there is exactly {e = +-1}: it depends on the PRODUCT only.
For n >= 3 that locus contains reciprocal-free rational points.  Hence P1 fails.

Then the same for r=2, by exact rational-root search on the numerator.
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B
import sympy as sp

def recipfree(w):
    n = len(w)
    return not any(w[i]*w[j] == 1 for i in range(n) for j in range(i, n))
def Sset(w):
    n = len(w); S = set()
    for i in range(n):
        for j in range(i, n):
            if w[i]*w[j] == 1: S |= {i, j}
    return S

print("=== (1)  Psi^C_{(1),(0)}(M) = 1 - e^{-2}  for every n  ===")
ok = bad = 0
for n in range(1, 8):
    for w in [[F(2),F(3),F(5),F(7),F(-2),F(-3),F(4)][:n],
              [F(3,2),F(5),F(-7,3),F(2),F(11),F(-5),F(4,3)][:n]]:
        if len(set(w)) != n: continue
        e = F(1)
        for x in w: e *= x
        lhs = psi_headsum(w, [0]); rhs = 1 - 1/e**2
        ok += (lhs == rhs); bad += (lhs != rhs)
print(f"    {ok}/{ok+bad}   (n=1..7, two families)")

print()
print("=== (2)  THE WITNESS:  r=1, n=3 = r+2, theta=(0), M = {2, 3, 1/6} ===")
W = [F(2), F(3), F(1,6)]
e = W[0]*W[1]*W[2]
print(f"    M = {W}      e = prod w_t = {e}")
print(f"    reciprocal pairs w_a w_b = 1 (a<=b): "
      f"{[(a+1,b+1) for a in range(3) for b in range(a,3) if W[a]*W[b]==1]}  -> NONE")
print(f"    reciprocal support S = {Sset(W)},  n - |S| = {3-len(Sset(W))} >= r = 1")
print(f"    Psi (head-injection sum, eq:Psidef) = {psi_headsum(W,[0])}")
print(f"    Psi (alternant expansion, thm:alt)  = {psi_alt(W,[0])}")
print(f"    B(T) over the three 1-subsets      = {[B(W,T) for T in itertools.combinations(range(3),1)]}")
print(f"    all B(T) nonzero: {all(B(W,T)!=0 for T in itertools.combinations(range(3),1))}")
print(f"    pairwise distinct m_j: vacuous at r=1 (single m_1 = 0)")

print()
print("=== (3)  a second witness at r=2 (so it is not an r=1 artefact) ===")
print("    r=2, n=4 = r+2, theta=(0,1) i.e. m=(0,2).  Fix w1,w2,w3 and look for")
print("    RATIONAL roots of numerator(Psi) in w4, then test reciprocal-freeness.")
w4 = sp.Symbol('w4', nonzero=True)
def psi_sym_num(wlist, theta):
    n = len(wlist); r = len(theta); m = [theta[j]+j for j in range(r)]
    tot = sp.Integer(0)
    for T in itertools.combinations(range(n), r):
        x = [wlist[i] for i in T]
        Bt = sp.Integer(1)
        for xi in x: Bt *= (1 - 1/xi**2)
        for a2 in range(r):
            for b2 in range(a2+1, r): Bt *= (1 - 1/(x[a2]*x[b2]))
        for t in range(n):
            if t in T: continue
            for xi in x: Bt *= (1 - 1/(wlist[t]*xi))/(1 - wlist[t]/xi)
        Mx = sp.Matrix(r, r, lambda i, j: x[i]**m[j])
        V = sp.Integer(1)
        for i in range(r):
            for k in range(i+1, r): V *= (x[k]-x[i])
        tot += Bt*Mx.det()/V
    return sp.fraction(sp.cancel(sp.together(tot)))[0]

found = []
import itertools as it
bases = [(2,3,5),(2,3,-5),(2,-3,5),(3,5,7),(2,5,7),(2,3,7),(1,2,3) if False else (4,5,7),
         (2,3,-7),(2,-5,7),(-2,3,5),(2,4,5),(3,4,5),(2,3,4),(2,6,7),(3,6,7)]
for (a_,b_,c_) in bases:
    wl = [sp.Rational(a_), sp.Rational(b_), sp.Rational(c_), w4]
    num = sp.expand(psi_sym_num(wl, [0,1]))
    P = sp.Poly(num, w4)
    for rt in sp.roots(P).keys():
        if not rt.is_rational: continue
        Wc = [F(a_), F(b_), F(c_), F(int(sp.numer(rt)), int(sp.denom(rt)))]
        if len(set(Wc)) != 4 or any(x == 0 for x in Wc): continue
        try:
            v1, v2 = psi_headsum(Wc,[0,1]), psi_alt(Wc,[0,1])
        except ZeroDivisionError: continue
        if v1 != 0 or v2 != 0: continue
        Bs = [B(Wc,T) for T in it.combinations(range(4),2)]
        found.append((Wc, recipfree(Wc), Sset(Wc), all(x != 0 for x in Bs), Bs))
print(f"    rational zeros found: {len(found)}")
for (Wc, rf, S, allnz, Bs) in found:
    print(f"      M = {Wc}   reciprocal-free: {rf}  S={S}  all B(T) nonzero: {allnz}")
    if rf and allnz:
        print(f"        Psi(headsum)={psi_headsum(Wc,[0,1])}  Psi(alternant)={psi_alt(Wc,[0,1])}")
        print(f"        m = (0,2) pairwise distinct; n=4 = r+2; n-|S| = {4-len(S)} >= r = 2")
        print(f"        B(T) over the six 2-subsets = {Bs}")

print()
print("=== (4) NEGATIVE CONTROL (must stay green): old Q559 sampler finds no zeros ===")
import random
random.seed(99)
VALS = [F(2),F(3),F(5),F(7),F(-2),F(-3),F(4),F(-5),F(11),F(-7),F(3,2),F(-5,3)]
zer = tot = 0
for r in range(1,4):
    for n in range(r, r+3):
        for _ in range(60):
            w = random.sample(VALS, n)
            # Q559 HYGIENE: exclude w=+-1, w_a w_b = 1, w_a = w_b
            if any(x*x == 1 for x in w): continue
            if any(w[i]*w[j] == 1 for i in range(n) for j in range(i+1,n)): continue
            if len(set(w)) != n: continue
            theta = list(range(r))
            try: v = psi_headsum(w, theta)
            except ZeroDivisionError: continue
            tot += 1; zer += (v == 0)
print(f"    hygienic sampler: {zer} zeros in {tot} draws   (expected 0; a zero here = engine regression)")
