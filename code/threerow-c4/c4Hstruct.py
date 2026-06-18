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

# For j=3..7: H(j) as polynomial in a,b. Reduce mod 2 and mod higher powers; find content.
print("=== H(j) mod 2 (a,b same parity) ===")
for jj in range(3,8):
    Hj = sp.Poly(sp.expand(H.subs(j_,jj)), a_,b_)
    # reduce coeffs mod 2,4,8
    for mod in [2,4,8,16,32]:
        cs = {mono:int(c)%mod for mono,c in Hj.terms()}
        cs = {k:v for k,v in cs.items() if v!=0}
        if not cs:
            print(f"  H({jj}) ≡ 0 mod {mod}  (v2 H >= {mod.bit_length()-1} as polynomial)")
        else:
            print(f"  H({jj}) mod {mod}: {cs}")
            break
