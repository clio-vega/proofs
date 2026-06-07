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
    if n==0: return 999
    v=0
    while n%2==0: n//=2; v+=1
    return v

for b in range(5, 25):
    Ib = Poly(ImP_coeff_b(b), m)
    R = (b-1)//2
    offsets=set()
    vals=[]
    for mm in range(b, b+40):
        Iv = int(Ib.eval(mm))
        if Iv==0:
            vals.append((mm,'ROOT')); continue
        prod = 1
        for r in range(0,R+1): prod*=(mm-r)
        off = v2(Iv)-v2(prod)
        offsets.add(off)
        vals.append((mm,off))
    print(f"b={b:2d} b%4={b%4}  v2(I_b)-v2(forcedprod) offsets over m in [b,b+40): {sorted(offsets)}")
