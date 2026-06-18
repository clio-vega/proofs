import sympy as sp
import sys
sys.path.insert(0,'/home/clio/projects/proofs/code/threerow-c3')
from mn import Mj
from math import comb

a,b,j = sp.symbols('a b j', integer=True)

# Q_4 numerator (second factor) from c4factor.py output:
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

def Mj_closed(av,bv,jv,mv):
    from fractions import Fraction
    N = 2*(mv-jv)
    num = comb(N, bv-jv) * (av-bv+1) * int(Q4.subs({a:av,b:bv,j:jv}))
    den = 24*(av+5-jv)
    for i in range(1,5):
        den *= (bv+i-jv)
    val = Fraction(num,den)
    return val

# 1. Verify closed form vs Murnaghan-Nakayama
print("=== Verify Lemma 1 (closed form) vs MN ===")
bad=0; tested=0
for mv in range(4,17):
    n=2*mv
    for av in range(4,n):
        bv = n-4-av
        if bv<4 or bv>av: continue
        lam=(av,bv,4)
        for jv in range(0,bv+1):
            M_mn = Mj(lam,jv,mv)
            M_cf = Mj_closed(av,bv,jv,mv)
            assert M_cf.denominator==1, (lam,jv,M_cf)
            tested+=1
            if int(M_cf)!=M_mn:
                bad+=1
                if bad<10: print("  MISMATCH",lam,jv,M_mn,int(M_cf))
print(f"tested={tested} mismatches={bad}")

# 2. Structural decomposition: Q4 = heavy*H + P_8
print("\n=== Structural decomposition ===")
P8 = sp.prod([j-t for t in range(8)])   # j(j-1)...(j-7) = 8! C(j,8)
heavyH = sp.expand(Q4 - P8)
print("Q4 - P8 factored:")
print(sp.factor(heavyH))

# Q4(0):
K4 = sp.expand(Q4.subs(j,0))
print("\nQ4(0)=K4 factored:", sp.factor(K4))
