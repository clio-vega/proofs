"""Job A (2026-06-19 CODE): independently de-risk the c=4 interior Number Lemma.

Validates the PROVE output 2026-06-19-c4-interior-number-lemma.md:
  Lemma N4:  16 | H  (v2 H >= 4, a CONSTANT floor, sharp)
and the deep-index absorption that closes the §6 residual Delta_hat(j) > 0, j>=8.

Strategy: tie H back to MN GROUND TRUTH.  Q_4 is MN-verified (c4verify.py);
we recover H := (Q_4 - P_8)/[(a-2)(b-3)] SYMBOLICALLY and compare to the
hardcoded Hval used in the proof's residue checks.  Only then do we census v2 H.
"""
import sympy as sp
import sys
sys.path.insert(0, '/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb

a, b, j = sp.symbols('a b j', integer=True)

# ---------------------------------------------------------------------------
# Q_4: the MN-verified octic numerator (from c4verify.py; itself checked vs MN).
# ---------------------------------------------------------------------------
Q4 = (a**4*b**4 + 6*a**4*b**3 - a**4*b**2 - 54*a**4*b - 72*a**4 + 10*a**3*b**4
 - 4*a**3*b**3*j**2 - 8*a**3*b**3*j + 60*a**3*b**3 - 24*a**3*b**2*j - 10*a**3*b**2
 + 28*a**3*b*j**2 + 80*a**3*b*j - 540*a**3*b + 24*a**3*j**2 + 192*a**3*j - 720*a**3
 + 23*a**2*b**4 - 12*a**2*b**3*j**2 - 48*a**2*b**3*j + 138*a**2*b**3 + 6*a**2*b**2*j**4
 - 12*a**2*b**2*j**3 + 30*a**2*b**2*j**2 - 144*a**2*b**2*j - 23*a**2*b**2 - 18*a**2*b*j**4
 + 60*a**2*b*j**3 - 6*a**2*b*j**2 + 504*a**2*b*j - 1242*a**2*b - 72*a**2*j**3 + 72*a**2*j**2
 + 1080*a**2*j - 1656*a**2 - 34*a*b**4 + 16*a*b**3*j**2 + 8*a*b**3*j - 204*a*b**3
 - 6*a*b**2*j**4 + 36*a*b**2*j**3 - 30*a*b**2*j**2 + 48*a*b**2*j + 34*a*b**2 - 4*a*b*j**6
 + 48*a*b*j**5 - 226*a*b*j**4 + 492*a*b*j**3 - 638*a*b*j**2 + 112*a*b*j + 1836*a*b
 + 12*a*j**6 - 144*a*j**5 + 732*a*j**4 - 1800*a*j**3 + 1752*a*j**2 - 984*a*j + 2448*a
 - 120*b**4 + 48*b**3*j**2 + 240*b**3*j - 720*b**3 - 12*b**2*j**4 - 24*b**2*j**3
 - 60*b**2*j**2 + 672*b**2*j + 120*b**2 + 8*b*j**6 - 96*b*j**5 + 524*b*j**4 - 1224*b*j**3
 + 1076*b*j**2 - 2880*b*j + 6480*b + j**8 - 28*j**7 + 298*j**6 - 1672*j**5 + 5305*j**4
 - 9244*j**3 + 9084*j**2 - 8928*j + 8640)

# ---------------------------------------------------------------------------
# (0) Re-verify Q_4 closed form against MN ground truth (a sample), so H rests
#     on something we re-checked here, not on trust.
# ---------------------------------------------------------------------------
def Mj_closed(av, bv, jv, mv):
    from fractions import Fraction
    N = 2*(mv-jv)
    num = comb(N, bv-jv) * (av-bv+1) * int(Q4.subs({a: av, b: bv, j: jv}))
    den = 24*(av+5-jv)
    for i in range(1, 5):
        den *= (bv+i-jv)
    return Fraction(num, den)

print("=== (0) Re-verify closed-form Q_4 vs Murnaghan-Nakayama ===")
bad = tested = 0
for mv in range(4, 15):
    n = 2*mv
    for av in range(4, n):
        bv = n-4-av
        if bv < 4 or bv > av:
            continue
        for jv in range(0, bv+1):
            M_mn = Mj((av, bv, 4), jv, mv)
            M_cf = Mj_closed(av, bv, jv, mv)
            assert M_cf.denominator == 1
            tested += 1
            if int(M_cf) != M_mn:
                bad += 1
print(f"    closed form vs MN: tested={tested}, mismatches={bad}  "
      f"({'OK' if bad == 0 else 'FAIL'})")

# ---------------------------------------------------------------------------
# (1) Recover H symbolically from MN-verified Q_4, and CROSS-CHECK against the
#     hardcoded Hval used by the proof's residue scripts (c4N4 / c4finite).
# ---------------------------------------------------------------------------
print("\n=== (1) Recover H = (Q_4 - P_8)/[(a-2)(b-3)] and cross-check ===")
P8 = sp.prod([j - t for t in range(8)])            # j(j-1)...(j-7) = 8! C(j,8)
heavyH = sp.expand(Q4 - P8)
H_quot, rem = sp.div(heavyH, sp.expand((a-2)*(b-3)), j)  # divide in j... use full div
# robust exact division over the polynomial ring:
H_recovered = sp.cancel(heavyH / ((a-2)*(b-3)))
H_recovered = sp.expand(H_recovered)
divides = sp.simplify(heavyH - (a-2)*(b-3)*H_recovered) == 0
print(f"    (a-2)(b-3) exactly divides Q_4 - P_8 ?  {divides}")

# The hardcoded H from c4N4.py / c4finite.py:
H_hard = (a**3*b**3 + 9*a**3*b**2 + 26*a**3*b + 24*a**3 + 12*a**2*b**3
 - 4*a**2*b**2*j**2 - 8*a**2*b**2*j + 108*a**2*b**2 - 12*a**2*b*j**2
 - 48*a**2*b*j + 312*a**2*b - 8*a**2*j**2 - 64*a**2*j + 288*a**2
 + 47*a*b**3 - 20*a*b**2*j**2 - 64*a*b**2*j + 423*a*b**2 + 6*a*b*j**4
 - 12*a*b*j**3 - 30*a*b*j**2 - 384*a*b*j + 1222*a*b + 24*a*j**3
 - 40*a*j**2 - 488*a*j + 1128*a + 60*b**3 - 24*b**2*j**2 - 120*b**2*j
 + 540*b**2 + 6*b*j**4 + 12*b*j**3 - 42*b*j**2 - 696*b*j + 1560*b
 - 4*j**6 + 48*j**5 - 244*j**4 + 648*j**3 - 664*j**2 - 648*j + 1440)
match = sp.expand(H_recovered - H_hard) == 0
print(f"    recovered H == hardcoded H (used in residue proofs) ?  {match}")

# H(0) factorization claimed in the note:
H0 = sp.factor(H_recovered.subs(j, 0))
print(f"    H(0) = {H0}")
print(f"      claimed (a+3)(a+4)(a+5)(b+2)(b+3)(b+4): "
      f"{sp.expand(H_recovered.subs(j,0) - (a+3)*(a+4)*(a+5)*(b+2)*(b+3)*(b+4)) == 0}")
print(f"    H is degree {sp.Poly(H_recovered, j).degree()} in j (claimed sextic=6)")

# P_8 sanity:
print(f"    P_8(j)=0 for j<=7 ?  {all(int(P8.subs(j,t)) == 0 for t in range(8))}; "
      f"P_8(8) = {int(P8.subs(j,8))} (=8!={3628800==int(P8.subs(j,8))})")
