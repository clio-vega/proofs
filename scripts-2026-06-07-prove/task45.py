from identify_sweep import *

print("="*70)
print("TASK 4: fingerprint sequences")
print("="*70)
for b in range(5, 21):
    Ib, Q, forced = deflate(b)
    coeffs, prim = primitive_int_poly(Q)
    pc = [int(c) for c in Poly(prim, m).all_coeffs()]
    Ivals = [int(Ib.subs(m, b + k)) for k in range(0, 9)] if all(
        (Ib.subs(m, b+k)).is_integer for k in range(9)) else [Ib.subs(m, b+k) for k in range(9)]
    print(f"\nb={b}")
    print(f"  primitive Q_b coeffs (high->low): {pc}")
    print(f"  I_b(b..b+8) for m={b}..{b+8}: {Ivals}")

print("\n\n"+"="*70)
print("TASK 5: discriminant of Q_b (primitive), b=5..14, factored")
print("="*70)
for b in range(5, 15):
    Ib, Q, forced = deflate(b)
    coeffs, prim = primitive_int_poly(Q)
    disc = sp.discriminant(Poly(prim, m))
    disc = sp.Integer(disc)
    fac = factorint(abs(disc)) if disc != 0 else {}
    # perfect square test
    is_sq = sp.sqrt(disc).is_integer if disc > 0 else False
    print(f"\nb={b}: disc(Q_b) = {disc}")
    print(f"   |disc| factorization: {fac}")
    print(f"   sign: {'+' if disc>0 else '-' if disc<0 else '0'}, perfect square? {bool(is_sq)}")
