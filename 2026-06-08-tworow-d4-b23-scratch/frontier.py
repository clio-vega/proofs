import sympy as sp
m,u=sp.symbols('m u')
s=1+sp.I

def g_half(b):
    # exact g_b((2b-1)/2) as complex
    mm=sp.Rational(2*b-1,2)
    def falling(x,k):
        r=sp.Integer(1)
        for i in range(k): r*=(x-i)
        return r
    return sum(falling(mm,b-l)/(sp.factorial(b-2*l)*sp.factorial(l))*s**(b-2*l)
               for l in range(0,b//2+1) if b-2*l>=0)

print("=== exact 4^b * G_b((2b-1)/2) as Gaussian integer ===")
for b in range(1,14):
    G = sp.expand(g_half(b)-g_half(b-1)) if b>=1 else g_half(0)
    val = sp.expand(G*sp.Integer(4)**b)
    cb = sp.binomial(2*b,b)
    print(f"b={b:2d}: 4^b G_b(half) = {val}   ; C(2b,b)={cb}, ratio={sp.nsimplify(val/cb)}")
