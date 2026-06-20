import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
from fractions import Fraction
import sympy as sp

a,b,j=sp.symbols('a b j', integer=True)

# Conjectured closed form denominator/prefactor:
# M_j = C(2(m-j), b-j) * (a-b+1) * Q5 / [ 120 * (a+6-j) * prod_{i=1..5}(b+i-j) ]
# => Q5 = M_j * 120 * (a+6-j)*prod(b+i-j) / [ C(2(m-j),b-j)*(a-b+1) ]

def Q5_value(A,B,J):
    m=(A+B+5)//2
    M=Mj((A,B,5),J,m)
    N=2*(m-J)
    denomfactor = 120*(A+6-J)
    for i in range(1,6): denomfactor*=(B+i-J)
    num = M*denomfactor
    den = comb(N,B-J)*(A-B+1)
    q=Fraction(num,den)
    assert q.denominator==1, (A,B,J,q)
    return int(q)

# sanity: Q5 integer for a bunch of shapes
import random
ok=True
for _ in range(200):
    # a+b odd, a>=b>=5+? need b-j>=0
    B=random.randrange(5,30)
    A=random.randrange(B, B+40)
    if (A+B)%2==0: A+=1   # need a+b odd
    if A<B: continue
    J=random.randrange(0, B+1)
    try:
        Q5_value(A,B,J)
    except AssertionError as e:
        ok=False; print("NONINT",e); break
print("Q5 integer test:", "PASS" if ok else "FAIL")

# Now interpolate Q5 as polynomial in j for FIXED (a,b) to find degree & tip.
# pick a representative (a,b)
A0,B0=41,30  # a+b=71 odd, b=30 large so j up to 30
pts=[(J, Q5_value(A0,B0,J)) for J in range(0,15)]
poly=sp.interpolate(pts, j)
poly=sp.expand(poly)
print("\nQ5(j) for (a,b)=(41,30): degree in j =", sp.degree(sp.Poly(poly,j)))
print("leading coeff:", sp.LC(sp.Poly(poly,j)))

# Check tip: Q5(j) - (a-3)(b-4)H5(j) should equal -10! C(j,10).
# For j=0..9, C(j,10)=0 => Q5(j) divisible by (a-3)(b-4); H5(j)=Q5(j)/((a-3)(b-4)).
print("\n(a-3)(b-4) =", (A0-3)*(B0-4))
for J in range(0,12):
    q=Q5_value(A0,B0,J)
    div=q/((A0-3)*(B0-4))
    print(f" j={J}: Q5={q}  Q5/((a-3)(b-4))={'%.4f'%div if div!=int(div) else int(div)}  C(j,10)={comb(J,10)}")
