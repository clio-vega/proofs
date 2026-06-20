import sys
sys.path.insert(0,'/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
from fractions import Fraction
import sympy as sp

a,b,j=sp.symbols('a b j')

def Q5_value(A,B,J):
    m=(A+B+5)//2
    M=Mj((A,B,5),J,m)
    N=2*(m-J)
    df=120*(A+6-J)
    for i in range(1,6): df*=(B+i-J)
    num=M*df; den=comb(N,B-J)*(A-B+1)
    q=Fraction(num,den); assert q.denominator==1
    return int(q)

def H5_value(A,B,J):  # for J<=9
    q=Q5_value(A,B,J)
    h=Fraction(q,(A-3)*(B-4)); assert h.denominator==1
    return int(h)

# interpolate H5(i)(a,b) for i=0..8 as poly in (a,b), deg<=5 each.
# grid of (A,B) with a+b odd, b>=10 (so b-i>=0 for i up to 8 and j-range fine), a>=b
def interp_ab(func, da=6, db=6):
    # collect points: need (da)*(db) points on a grid; ensure a+b odd
    pts=[]
    Bvals=[10,11,12,13,14,15,16][:db]
    Avals_base=[20,21,22,23,24,25,26,27][:da]
    # build value table
    # Use sympy interpolate over 2D via successive 1D interpolation
    # For fixed B, interpolate in A; then interpolate those poly-coeffs in B.
    polysA=[]
    for B0 in Bvals:
        ptsA=[]
        for A0 in Avals_base:
            A=A0
            if (A+B0)%2==0: A=A0  # need odd sum
            # ensure parity: a+b odd
            if (A+B0)%2==0: A=A0+1
            ptsA.append((A, func(A,B0)))
        # ensure distinct A
        xs=[p[0] for p in ptsA]
        if len(set(xs))<len(xs):
            ptsA=[(p[0]+2*k,func(p[0]+2*k,B0)) for k,p in enumerate(ptsA)]
        pA=sp.interpolate(ptsA, a)
        polysA.append((B0, sp.expand(pA)))
    # now interpolate coefficients in B
    P=sp.Poly(polysA[0][1], a)
    # get monomials in a present across all
    maxdeg=max(sp.Poly(p,a).degree() for _,p in polysA)
    result=0
    for deg in range(maxdeg+1):
        ptsB=[]
        for B0,p in polysA:
            coeff=sp.Poly(p,a).coeff_monomial(a**deg) if deg>0 else sp.Poly(p,a).coeff_monomial(1)
            ptsB.append((B0, coeff))
        pB=sp.interpolate(ptsB, b)
        result+= sp.expand(pB)*a**deg
    return sp.expand(result)

# Build H5 in C(j,k) basis: h_k = sum_{i=0}^k (-1)^(k-i) C(k,i) H5(i)
Hvals=[]
for i in range(0,9):
    Hi=interp_ab(lambda A,B,I=i: H5_value(A,B,I))
    Hvals.append(Hi)
    print(f"H5({i}) interpolated, deg_a={sp.degree(sp.Poly(Hi,a))}, deg_b={sp.degree(sp.Poly(Hi,b))}")

hcoeffs=[]
for k in range(0,9):
    hk=sum((-1)**(k-i)*comb(k,i)*Hvals[i] for i in range(0,k+1))
    hk=sp.expand(hk)
    hcoeffs.append(hk)

print("\nH5(0) =", sp.factor(hcoeffs[0]))
print("expected (a+3)(a+4)(a+5)(a+6)(b+2)(b+3)(b+4)(b+5)")
exp0=sp.expand((a+3)*(a+4)*(a+5)*(a+6)*(b+2)*(b+3)*(b+4)*(b+5))
print("H5(0) matches expected:", sp.expand(hcoeffs[0]-exp0)==0)

# Reconstruct H5(j) symbolic
def binom_poly(jj,k):
    r=sp.Integer(1)
    for t in range(k): r*= (jj-t)
    return r/sp.factorial(k)
H5sym=sum(hcoeffs[k]*binom_poly(j,k) for k in range(9))
H5sym=sp.expand(H5sym)

# Full Q5 reconstruction and verify vs MN
Q5sym=sp.expand((a-3)*(b-4)*H5sym - sp.factorial(10)*binom_poly(j,10))

import random
mism=0;tested=0
for _ in range(400):
    B=random.randrange(10,40); A=random.randrange(B,B+50)
    if (A+B)%2==0: A+=1
    if A<B: continue
    J=random.randrange(0,B+1)
    qexact=Q5_value(A,B,J)
    qsym=int(Q5sym.subs({a:A,b:B,j:J}))
    tested+=1
    if qexact!=qsym: mism+=1; 
    if mism<5 and qexact!=qsym: print("MISMATCH",A,B,J,qexact,qsym)
print(f"\nQ5 symbolic vs MN: tested={tested} mismatches={mism}")

# save H5 coeffs for later
import pickle
with open('/home/clio/projects/code/threerow-c5/H5sym.pkl','wb') as f:
    pickle.dump({'hcoeffs':[sp.srepr(h) for h in hcoeffs],'H5sym':sp.srepr(H5sym),'Q5sym':sp.srepr(Q5sym)}, f)
print("saved.")
for k in range(9):
    print(f"\nh_{k} = ", sp.factor(hcoeffs[k]))
