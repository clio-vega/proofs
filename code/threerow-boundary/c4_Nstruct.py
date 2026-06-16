import sympy as sp
P,b = sp.symbols('P b')
a = 2*P+b+4
N3 = a*(b+4)-(b**2+3*b+8)
N2 = a**2*b**2+7*a**2*b+12*a**2-2*a*b**3-9*a*b**2-21*a*b-20*a+b**4+2*b**3+13*b**2+16*b-8
N1 = a**2*b**3+9*a**2*b**2+26*a**2*b+24*a**2-2*a*b**4-7*a*b**3-13*a*b**2-14*a*b+b**5-2*b**4+17*b**3-4*b**2-84*b-96
for name,N in [('N3',N3),('N2',N2),('N1',N1)]:
    Ne = sp.expand(N)
    print(name,"=",sp.factor(Ne))
    print("   expanded:",Ne)
    # substitute b=2B (b even) and b=2B+1 (b odd), factor out 2s
    B=sp.symbols('B')
    Nev = sp.expand(Ne.subs(b,2*B))
    Nod = sp.expand(Ne.subs(b,2*B+1))
    cev = sp.gcd(sp.Poly(Nev,P,B).coeffs()) if Nev!=0 else 0
    cod = sp.gcd(sp.Poly(Nod,P,B).coeffs()) if Nod!=0 else 0
    print("   b even (b=2B): content gcd =",cev, " ->", sp.factor(Nev))
    print("   b odd  (b=2B+1): content gcd =",cod, " ->", sp.factor(Nod))
    print()
