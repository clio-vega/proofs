from identify_sweep import *

print("="*70)
print("TASK 1: Q_b in primitive integer form, b=5..20")
print("="*70)
data = {}
for b in range(5, 21):
    Ib, Q, forced = deflate(b)
    coeffs, prim = primitive_int_poly(Q)
    deg = sp.degree(prim, m)
    pc = Poly(prim, m).all_coeffs()
    lead = pc[0]
    const = pc[-1]
    data[b] = dict(prim=prim, coeffs=pc, lead=lead, const=const, forced=forced, deg=deg)
    print(f"\nb={b}: forced roots {forced},  deg Q_b = {deg}")
    print(f"  Q_b (primitive) = {prim}")
    print(f"  coeff list (high->low): {pc}")
    print(f"  leading coeff = {lead}   factorization: {factorint(abs(lead))}")
    print(f"  constant term = {const}  factorization: {factorint(abs(const)) if const!=0 else 'ZERO'}")

print("\n\n"+"="*70)
print("TASK 2: factor Q_b over Q, b=5..16  (check 4|b => (2m-(2b-1)) splits)")
print("="*70)
for b in range(5, 17):
    prim = data[b]['prim']
    fac = sp.factor(prim)
    factors = sp.Mul.make_args(fac)
    nontriv = [f for f in factors if not f.is_number]
    irred = (len(nontriv) == 1 and sp.degree(nontriv[0], m) == data[b]['deg'])
    print(f"\nb={b} (deg {data[b]['deg']}): factored = {fac}")
    print(f"   irreducible over Q? {irred}")
    if b % 4 == 0:
        lin = (2*m - (2*b-1))
        # does (2m-(2b-1)) divide prim?
        q, r = sp.div(Poly(prim, m), Poly(lin, m))
        print(f"   4|b: (2m-(2b-1))=(2m-{2*b-1}) divides Q_b? {r.as_expr()==0}")
