"""
Is the chain identity  log W = d Xi  genuinely false, or only false under new labels?
Clio, 2026-09-03. Built from Day 152 section 1 definitions ONLY.

Day 152 section 1 (verbatim):
    T^+ : u^alpha -> prod_i u_i^(alpha_i)   (RISING factorials)
    Psi^+(f) = T^+(f V)/V,   V = prod_{i<j}(u_i - u_j)
    F_P := Psi^+(exp(T e_2)) = sum_b Psi^+(e_2^b) T^b/b!
    tau  = automorphism induced by u_i -> u_i + 1
    H   := tau(F_P)/F_P
    d   := sum_i partial_{u_i}            <-- NOT partial_T
    wt(u_i)=+1, wt(t_i)=-1; ell^top_w(X)[T^n] = u-homogeneous part of degree n+w
    Xi  := ell^top_1(log F_P)              (Day 152 line 111)
    (P1): ell^top_0(H) = exp(d Xi),  i.e.  log W = d Xi     (Day 152 Theorem A)
"""
import sympy as sp

u1, u2, u3 = sp.symbols('u1 u2 u3')
U = (u1, u2, u3)
NB = 7                                   # build F_P through T^NB

V = sp.expand((u1 - u2) * (u1 - u3) * (u2 - u3))
e2 = u1*u2 + u1*u3 + u2*u3

def rising(x, m):
    return sp.prod([x + j for j in range(m)]) if m > 0 else sp.Integer(1)

def Tplus(poly):
    """u^alpha -> prod_i rising(u_i, alpha_i), extended linearly."""
    p = sp.Poly(sp.expand(poly), u1, u2, u3)
    out = sp.Integer(0)
    for alpha, co in p.terms():
        out += co * sp.prod([rising(U[i], alpha[i]) for i in range(3)])
    return sp.expand(out)

def PsiPlus(f):
    num = Tplus(sp.expand(f * V))
    qt, rem = sp.div(sp.Poly(num, u1), sp.Poly(V, u1))
    assert rem.is_zero, "Psi^+ : V does not divide"
    return sp.expand(qt.as_expr())

print("Building F_P = sum_b Psi^+(e_2^b) T^b/b!  ...")
FP = []
e2b = sp.Integer(1)
for b in range(NB + 1):
    FP.append(sp.expand(PsiPlus(e2b) / sp.factorial(b)))
    e2b = sp.expand(e2b * e2)
    print(f"   b={b} done", end="\r")
print("\n  [T^0]F_P =", FP[0], " [T^1]F_P =", sp.expand(FP[1]))

N = NB
def mul(a, b):
    return [sp.expand(sum(a[i]*b[n-i] for i in range(n+1))) for n in range(N+1)]
def inv(a):
    r = [sp.Integer(0)]*(N+1); r[0] = sp.Integer(1)
    for n in range(1, N+1):
        r[n] = sp.expand(-sum(a[j]*r[n-j] for j in range(1, n+1)))
    return r
def logser(a):
    ap = [(m+1)*a[m+1] for m in range(N)] + [sp.Integer(0)]
    d = mul(ap, inv(a))
    return [sp.Integer(0)] + [sp.expand(d[n-1]/n) for n in range(1, N+1)]
def homog(poly, deg):
    if sp.expand(poly) == 0: return sp.Integer(0)
    p = sp.Poly(sp.expand(poly), u1, u2, u3)
    return sp.expand(sum(co*u1**a*u2**b*u3**c for (a,b,c), co in p.terms() if a+b+c == deg))
def topdeg(poly):
    if sp.expand(poly) == 0: return None
    p = sp.Poly(sp.expand(poly), u1, u2, u3)
    return max(a+b+c for (a,b,c), co in p.terms() if co != 0)

X = logser(FP)                                    # X = log F_P
print("\n--- Fact II(c): deg_u [T^n] log F_P  (Day 152 claims = n+1, equality) ---")
for n in range(1, N+1):
    print(f"   n={n}: deg_u = {topdeg(X[n])}   n+1 = {n+1}   "
          f"{'OK' if topdeg(X[n])==n+1 else '*** MISMATCH ***'}")

Xi   = [homog(X[n], n+1) for n in range(N+1)]     # Xi  = ell^top_1(log F_P)
X0   = [homog(X[n], n)   for n in range(N+1)]     # X^(0) = ell^top_0(log F_P)

# H = tau(F_P)/F_P
tauFP = [sp.expand(f.subs({u1:u1+1, u2:u2+1, u3:u3+1}, simultaneous=True)) for f in FP]
H = mul(tauFP, inv(FP))
Wcal = [homog(H[n], n) for n in range(N+1)]       # ell^top_0(H)
print("\n  ell^top_0(H) coefficients (Day 152 line 372 pre-registers these):")
E1s, E2s, E3s = sp.symbols('E1 E2 E3')
def to_E(p):
    r = sp.symmetrize(sp.expand(p), U, formal=False) if hasattr(sp,'symmetrize') else None
    return sp.factor(sp.expand(p))
for n in range(4):
    print(f"    [T^{n}] ell_0(H) = {sp.factor(Wcal[n])}")

logW = logser(Wcal)
dXi  = [sp.expand(sum(sp.diff(Xi[n], v) for v in U)) for n in range(N+1)]

print("\n" + "="*78)
print("(P1)  log( ell^top_0 H )  =?=  d Xi,   d = sum_i partial_{u_i}   [3 VARIABLES, E_3 free]")
print("="*78)
p1_ok = True
for n in range(1, N+1):
    r = sp.expand(logW[n] - dXi[n])
    ok = (r == 0); p1_ok &= ok
    print(f"   n={n}: {'MATCH' if ok else 'DIFFERS: ' + str(sp.factor(r))}")
print(f"  => (P1) holds as written (partial = sum_i d/du_i), n=1..{N}: {p1_ok}")

print("\n" + "="*78)
print("Now the SAME identity with partial misread as partial_T :  log W =?= dXi/dT")
print("="*78)
dTXi = [sp.expand((n+1)*Xi[n+1]) for n in range(N)] + [sp.Integer(0)]
for n in range(1, 4):
    r = sp.expand(logW[n] - dTXi[n])
    print(f"   n={n}: {'MATCH' if r==0 else 'DIFFERS'}   "
          f"logW={sp.factor(logW[n])}   dXi/dT={sp.factor(dTXi[n])}")

print("\n" + "="*78)
print("RESTRICT TO u_3 = 0 and compare against Day 158's objects")
print("="*78)
z = {u3: 0}
print("  Does ell^top_0(H)|_{u3=0} equal Day 158's W = Y/(Tq)?")
# rebuild Day-158 W in u1,u2
E1 = u1+u2; E2 = u1*u2
Y = [sp.Integer(0)]*(N+1)
for _ in range(N+2):
    Ysq = [sp.expand(sum(Y[i]*Y[j-i] for i in range(j+1))) for j in range(N+1)]
    Y = [sp.Integer(0)] + [sp.expand((1 if n==1 else 0) + E1*Y[n-1] + E2*Ysq[n-1])
                            for n in range(1, N+1)]
qs = [sp.Integer(1)] + [sp.expand((-E1 if n==1 else 0) - 2*E2*Y[n-1]) for n in range(1, N+1)]
YoT = [sp.expand(Y[n+1]) for n in range(N)] + [sp.Integer(0)]
W158 = mul(YoT, inv(qs))
same = all(sp.expand(Wcal[n].subs(z) - W158[n]) == 0 for n in range(N))
print(f"    ell^top_0(H)|_u3=0  ==  Y/(Tq)  for n<={N-1}:  {same}")

print("\n  Day 158's claim: 'log W = d Xi is FALSE (checked at n=2)'.")
print("  Test the restricted identity, d taken as sum_i d/du_i THEN set u3=0:")
for n in range(1, 5):
    r = sp.expand(logW[n].subs(z) - dXi[n].subs(z))
    print(f"     n={n}: {'MATCH' if r==0 else 'DIFFERS: '+str(sp.factor(r))}")
print("\n  And with d/du_3 DROPPED (i.e. differentiating only the u3=0 slice):")
dXi_2var = [sp.expand(sp.diff(Xi[n].subs(z), u1) + sp.diff(Xi[n].subs(z), u2)) for n in range(N+1)]
for n in range(1, 5):
    r = sp.expand(logW[n].subs(z) - dXi_2var[n])
    print(f"     n={n}: {'MATCH' if r==0 else 'DIFFERS: '+str(sp.factor(r))}")
