import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb
import sympy as sp

a,b,j=sp.symbols('a b j')
# q_j := M_j*(a+3-j)(b+2-j)(b+1-j)/C(N,b-j). Conjecture: polynomial in a,b,j. Fit it.
def qval(A,B,J):
    m=(A+B+2)//2
    N=2*(m-J); k=B-J
    if k<0 or k>N: return None
    cb=comb(N,k)
    if cb==0: return None
    from fractions import Fraction
    val=Fraction(Mj((A,B,2),J,m))*(A+3-J)*(B+2-J)*(B+1-J)/cb
    return val

# collect data
data=[]
for A in range(3,16):
    for B in range(2,A+1):
        if (A+B)%2!=0: continue  # a+b=2m-2 even
        for J in range(0,B+1):
            q=qval(A,B,J)
            if q is not None:
                data.append((A,B,J,q))
print("data points:",len(data))

# fit polynomial in a,b,j up to total degree 5
from itertools import product
monos=[]
for da in range(0,5):
    for db in range(0,5):
        for dj in range(0,5):
            if da+db+dj<=5:
                monos.append((da,db,dj))
import sympy
# build linear system over rationals
M=sp.Matrix([[sp.Rational(A)**da*sp.Rational(B)**db*sp.Rational(J)**dj for (da,db,dj) in monos] for (A,B,J,q) in data])
rhs=sp.Matrix([sp.Rational(q.numerator,q.denominator) for (A,B,J,q) in data])
sol=M.solve_least_squares(rhs) if M.rows>=len(monos) else None
# check exact
coeffs=sol
poly=sum(coeffs[i]*a**da*b**db*j**dj for i,(da,db,dj) in enumerate(monos))
poly=sp.expand(poly)
# verify
bad=0
for (A,B,J,q) in data:
    v=poly.subs({a:A,b:B,j:J})
    if sp.nsimplify(v-sp.Rational(q.numerator,q.denominator))!=0: bad+=1
print("fit residual bad=",bad)
print("q_j =",sp.factor(poly))
