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

# Build the divisibility polynomials Phi_j and required 2^k.
# Delta(j) >0 reduces to 2^k | Phi_j(a,b) for all a==b mod 2.
def fall(top, r):  # symbolic falling factorial top*(top-1)*...*(top-r+1)
    e=sp.Integer(1)
    for t in range(r): e*= (top - t)
    return e

cases = {
 3: (sp.expand(H.subs(j_,3)), 3),                       # 8 | H(3)
 4: (sp.expand((a_+2)*(b_+1)*H.subs(j_,4)), 6),         # 64 | (a+2)(b+1)H(4)
 5: (sp.expand(fall(a_+2,2)*fall(b_+1,2)*H.subs(j_,5)), 6),  # 64
 6: (sp.expand(fall(a_+2,3)*fall(b_+1,3)*H.subs(j_,6)), 8),  # 256
 7: (sp.expand(fall(a_+2,4)*fall(b_+1,4)*H.subs(j_,7)), 8),  # 256
}

def residue_check(poly, k):
    M = 1<<k
    f = sp.lambdify((a_,b_), poly, 'math')
    worst = M
    for a in range(M):
        for b in range(a%2, M, 2):   # b == a mod 2
            v = int(f(a,b)) % M
            if v!=0:
                return False, (a,b,v)
    return True, None

for j,(poly,k) in cases.items():
    ok,info = residue_check(poly,k)
    print(f"j={j}: 2^{k} | Phi_j  for all a==b mod 2 ?  {'YES (PROVEN)' if ok else 'NO '+str(info)}")
print("\nAll finite-case divisibilities are complete residue checks => rigorous.")
