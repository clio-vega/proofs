import sympy as sp
from sympy import symbols, expand, Poly, I, im
m = symbols('m', real=True)

def ImP_coeff_b(b):
    s = 1 + I
    def rl(l):
        tot = 0
        for c in range(0, l//2 + 1):
            k = l - 2*c
            ff = 1
            for j in range(0, k + c):
                ff *= (m - j)
            tot += ff / (sp.factorial(k) * sp.factorial(c)) * s**k
        return sp.expand(tot)
    return expand(im(rl(b) - rl(b-1)))

def v2(n):
    n=int(n)
    if n==0: return None
    v=0
    while n%2==0: n//=2; v+=1
    return v

import math
for b in range(5, 25):
    Ib = Poly(ImP_coeff_b(b), m)
    R = (b-1)//2
    vRfact = v2(math.factorial(R))
    match=True
    fails=[]
    for mm in range(b, b+60):
        Iv = int(Ib.eval(mm))
        prod = 1
        for r in range(0,R+1): prod*=(mm-r)
        lhs = v2(Iv)             # v2(I_b(m)); None if 0
        rhs = v2(prod) - vRfact  # v2((m)_{R+1}) - v2(R!)
        if lhs != rhs:
            match=False
            fails.append((mm, lhs, rhs))
    tag = "MATCH (j=1 dominates -> offset = -v2(R!) const)" if match else f"differs at {fails[:4]}"
    print(f"b={b:2d} b%4={b%4} R={R} -v2(R!)={-vRfact}: {tag}")
