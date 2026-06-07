import sympy as sp
m = sp.symbols('m', real=True)
I = sp.I

# Recompute I_b(m) (degree-b poly in m), exact.
def coeff_ImP(l):
    total = sp.Integer(0); c = 0
    while l - 2*c >= 0:
        a = l - 2*c; ac = a + c
        fall = sp.Integer(1)
        for j in range(ac): fall *= (m - j)
        total += fall/(sp.factorial(a)*sp.factorial(c)) * (1+I)**a
        c += 1
    return sp.im(sp.expand(total))

def Ib(b):
    return sp.expand(coeff_ImP(b) - coeff_ImP(b-1))

# polynomial binomial C(z,r) = z(z-1)...(z-r+1)/r!
def pbin(z, r):
    if r < 0: return sp.Integer(0)
    num = sp.Integer(1)
    for t in range(r): num *= (z - t)
    return num/sp.factorial(r)

# Krawtchouk (binary q=2): K_b(x;N)=sum_j (-1)^j C(x,j)C(N-x,b-j)
def Kraw(b, N, x):
    out = sp.Integer(0)
    for j in range(b+1):
        out += (-1)**j * pbin(x, j) * pbin(N - x, b - j)
    return sp.expand(out)

# Try to match I_b(m) to c * K_b(m; N) for integer N, rational c.
for b in range(5, 11):
    Ipoly = Ib(b)
    pI = sp.Poly(Ipoly, m)
    degI = pI.degree()
    print(f"\nb={b}, deg I_b={degI}")
    # search N
    found=False
    for N in range(0, 60):
        K = sp.Poly(Kraw(b, N, m), m)
        if K.degree()!=degI:
            continue
        # ratio of leading coeffs
        lc = pI.LC()/K.LC()
        if sp.simplify(sp.expand(Ipoly - lc*K.as_expr()))==0:
            print(f"   MATCH: I_b = ({lc}) * K_{b}(m; N={N})")
            found=True
            break
    if not found:
        print("   no binary-Krawtchouk K_b(m;N) match for N in 0..59")
