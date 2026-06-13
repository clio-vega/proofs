import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
import sympy as sp
a,b,m=sp.symbols('a b m')
# Fit M_{b+1}, M_{b+2}, M_{b+3} as polynomials in a,b (with m=(a+b+2)/2)
def fitpoly(jfun, name, deg=4):
    data=[]
    for A in range(3,30):
        for B in range(2,A+1):
            if (A+B)%2: continue
            M=(A+B+2)//2
            J=jfun(B)
            if J>M: continue
            data.append((A,B,Mj((A,B,2),J,M)))
    # fit poly in a,b total deg <=deg
    monos=[(da,db) for da in range(deg+1) for db in range(deg+1) if da+db<=deg]
    Mat=sp.Matrix([[sp.Rational(A)**da*sp.Rational(B)**db for (da,db) in monos] for (A,B,_) in data])
    rhs=sp.Matrix([sp.Rational(v) for (_,_,v) in data])
    sol=Mat.solve_least_squares(rhs)
    poly=sp.expand(sum(sol[i]*a**da*b**db for i,(da,db) in enumerate(monos)))
    bad=sum(1 for (A,B,v) in data if poly.subs({a:A,b:B})!=v)
    print(f"{name}: bad={bad}  M = {sp.factor(poly)}")
    return poly
fitpoly(lambda B:B+1,"M_{b+1}")
fitpoly(lambda B:B+2,"M_{b+2}")
fitpoly(lambda B:B+3,"M_{b+3}",deg=2)
# b=2 explicit: M_0,M_2,M_4
print("\n--- b=2 family: lambda=(a,2,2), a=2m-4 even ---")
for A in range(4,40,2):
    M=(A+4)//2
    print(f"a={A} m={M}: M0={Mj((A,2,2),0,M)} M2={Mj((A,2,2),2,M)} M3={Mj((A,2,2),3,M)} M4={Mj((A,2,2),4,M)} M5={Mj((A,2,2),5,M) if 5<=M else '-'}")
