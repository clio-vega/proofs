import sympy as sp
a_,b_,j_ = sp.symbols('a b j', integer=True)
H = (a_**3*b_**3 + 9*a_**3*b_**2 + 26*a_**3*b_ + 24*a_**3 + 12*a_**2*b_**3
 - 4*a_**2*b_**2*j_**2 - 8*a_**2*b_**2*j_ + 108*a_**2*b_**2 - 12*a_**2*b_*j_**2
 - 48*a_**2*b_*j_ + 312*a_**2*b_ - 8*a_**2*j_**2 - 64*a_**2*j_ + 288*a_**2
 + 47*a_*b_**3 - 20*a_*b_**2*j_**2 - 64*a_*b_**2*j_ + 423*a_*b_**2 + 6*a_*b_*j_**4
 - 12*a_*b_*j_**3 - 30*a_*b_*j_**2 - 384*a_*b_*j_ + 1222*a_*b_ + 24*a_*j_**3
 - 40*a_*j_**2 - 488*a_*j_ + 1128*a_ + 60*b_**3 - 24*b_**2*j_**2 - 120*b_**2*j_
 + 540*b_**2 + 6*b_*j_**4 + 12*b_*j_**3 - 42*b_*j_**2 - 696*b_*j_ + 1560*b_
 - 4*j_**6 + 48*j_**5 - 244*j_**4 + 648*j_**3 - 664*j_**2 - 648*j_ + 1440)

def content2(expr, var):
    # gcd of integer coeffs' 2-adic val of polynomial in var (min v2 of coeffs)
    p=sp.Poly(sp.expand(expr), var)
    cs=[int(c) for c in p.all_coeffs()]
    def v2(n):
        n=abs(n)
        if n==0: return 99
        k=0
        while n%2==0:n//=2;k+=1
        return k
    return min(v2(c) for c in cs)

print("=== H(j)|_{a=-2}, |_{b=-1}, |_{both}  (factored) and 2-content ===")
for jj in range(3,8):
    Hj=sp.expand(H.subs(j_,jj))
    Ha=sp.expand(Hj.subs(a_,-2))
    Hb=sp.expand(Hj.subs(b_,-1))
    Hab=sp.expand(Ha.subs(b_,-1))
    print(f"\n j={jj}:")
    print("   H(j)|a=-2 =",sp.factor(Ha), " 2content(in b)=",content2(Ha,b_) if Ha!=0 else 'ZERO')
    print("   H(j)|b=-1 =",sp.factor(Hb), " 2content(in a)=",content2(Hb,a_) if Hb!=0 else 'ZERO')
    print("   H(j)|a=-2,b=-1 =",Hab," v2=",(lambda n:(99 if n==0 else (len(bin(abs(n)))-len(bin(abs(n)).rstrip('0'))) ))(int(Hab)))
