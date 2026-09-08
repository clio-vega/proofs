"""CHECK 7 -- THE VARYING TEST FOR STEP 2.
Unpin a.  Take Htilde_{a,b}(z1) and Htilde_{1, a/b}(z2).  Then the total plethystic
shift cancels on the diagonal z1 = b z2 for EVERY a, and the diagonal sum is

      sum_l b^l Phi'_{l, N-l} f  =  h_N[(1+b) X] . f

The twist alphabets are (a - b) and (1 - a/b): both move with a.
The defect alphabet is (1 + b): it does not.   a=1 recovers Q99 Lemma 3.3.
"""
import sympy as sp
from vertex import *
from check6 import shift2, Phi
a, b = sp.symbols('a b')

cf1 = lambda k: a**k - b**k
cf2 = lambda k: 1 - (a/b)**k

ok = bad = 0
for deg in range(0, 4):
    for mu in parts(deg):
        f = {mu: sp.Integer(1)}
        d = deg
        for N in range(-1, 4):
            # sum over l with Phi_{l,N-l} f != 0 : needs l >= -d and N-l >= -d
            tot = {}
            for l in range(-d, N + d + 1):
                P = Phi(f, l, N-l, cf1, cf2)
                for k, v in P.items():
                    tot[k] = sp.expand(tot.get(k, 0) + b**l * v)
            tot = {k: sp.cancel(sp.together(v)) for k, v in tot.items()}
            tot = {k: v for k, v in tot.items() if v != 0}
            pred = pmul(h_alphabet(N, lambda k: 1 + b**k), f) if N >= 0 else {}
            D = sub(tot, pred)
            D = {k: sp.simplify(sp.cancel(v)) for k, v in D.items()}
            D = {k: v for k, v in D.items() if v != 0}
            if D: bad += 1; print("FAIL", mu, N, D)
            else: ok += 1
print("diagonal sum = h_N[(1+b)X] with a FREE:", ok, "agree,", bad, "mismatch")
print("  (the defect alphabet contains no 'a'; the two twist alphabets are a-b and 1-a/b)")

# explicit demonstration: fix the pole b=t, move the zero a.
t = sp.symbols('t')
print("\n  b = t fixed, a varied:")
for aval in [1, -1, 2, sp.Rational(1,2)]:
    tw1 = [sp.expand((aval**k - t**k)) for k in (1,2,3)]
    tw2 = [sp.expand(1 - (aval/t)**k) for k in (1,2,3)]
    df  = [sp.expand(1 + t**k) for k in (1,2,3)]
    print(f"   a={aval}:  twist1={tw1}  twist2={tw2}   defect={df}")
