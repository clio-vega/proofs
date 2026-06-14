import sympy as sp
from math import comb
from mn import Mj

a,b,j = sp.symbols('a b j', integer=True)

# Closed form from c3factor.py
Q3 = (a**3*b**3 + 3*a**3*b**2 - 4*a**3*b - 12*a**3 + 6*a**2*b**3 - 3*a**2*b**2*j**2
      - 3*a**2*b**2*j + 18*a**2*b**2 + 3*a**2*b*j**2 - 3*a**2*b*j - 24*a**2*b
      + 6*a**2*j**2 + 18*a**2*j - 72*a**2 + 5*a*b**3 - 3*a*b**2*j**2 - 9*a*b**2*j
      + 15*a*b**2 + 3*a*b*j**4 - 12*a*b*j**3 + 24*a*b*j**2 - 15*a*b*j - 20*a*b
      - 6*a*j**4 + 24*a*j**3 - 36*a*j**2 + 66*a*j - 60*a - 12*b**3 + 6*b**2*j**2
      + 12*b**2*j - 36*b**2 - 3*b*j**4 + 12*b*j**3 - 27*b*j**2 + 18*b*j + 48*b
      - j**6 + 15*j**5 - 79*j**4 + 201*j**3 - 244*j**2 + 36*j + 144)

def Mj_closed(av,bv,jv,mv):
    N=2*(mv-jv)
    Bb=bv-jv
    if Bb<0: return None  # interior only
    base=comb(N,Bb)
    val = base*(av-bv+1)*Q3.subs({a:av,b:bv,j:jv}) / (6*(av-jv+4)*(bv-jv+1)*(bv-jv+2)*(bv-jv+3))
    return sp.Rational(val)

# verify interior 0<=j<=b
bad=0; tested=0
for mv in range(3,25):
    n=2*mv
    for av in range(1,n+1):
        for bv in range(3,av+1):       # c=3 needs b>=3
            cv=n-av-bv
            if cv!=3: continue
            for jv in range(0,bv+1):
                cl=Mj_closed(av,bv,jv,mv)
                ex=Mj((av,bv,3),jv,mv)
                tested+=1
                if cl is None or sp.Rational(ex)!=cl:
                    bad+=1
                    if bad<20: print("MISMATCH",(av,bv,3),"m",mv,"j",jv,"closed",cl,"exact",ex)
print(f"interior closed-form check: tested={tested} mismatches={bad}")

# Analyze sextic: Q3 as poly in j, K3=Q3(j=0), perturbation
K3 = sp.expand(Q3.subs(j,0))
pert = sp.expand(Q3 - K3)
print("\nK3 = Q3(j=0) =", sp.factor(K3))
print("\nQ3 - K3 collected in j:")
print(sp.collect(pert,j))
# does pert contain 720*C(j,6) = j(j-1)(j-2)(j-3)(j-4)(j-5) ?
P6 = sp.expand(j*(j-1)*(j-2)*(j-3)*(j-4)*(j-5))
# coefficient of pure-j^6..: the a,b-free part of pert
abfree = sp.expand(pert.subs({a:0,b:0}))
print("\na,b-free part of (Q3-K3):", abfree)
print("compare -P6 = ", sp.expand(-P6))
print("difference abfree-(-P6) =", sp.expand(abfree+P6))
