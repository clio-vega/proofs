import sympy as sp
from sympy import symbols, expand, Poly, I, im
import math, random
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
    if n==0: return 10**9
    v=0
    while n%2==0: n//=2; v+=1
    return v
for b in [8,12,16,20,24]:
    R=(b-1)//2; vR=v2(math.factorial(R))
    Ib=Poly(ImP_coeff_b(b),m)
    ok=True; fails=[]
    for _ in range(500):
        mm=random.randint(b,30000)
        Iv=int(Ib.eval(mm))
        prod=1
        for r in range(R+1): prod*=(mm-r)
        if v2(Iv)!=v2(prod)-vR: ok=False; fails.append(mm)
    print(f"b={b} (0 mod4): aggregate v2 identity {'OK' if ok else 'FAIL '+str(fails[:3])}")
