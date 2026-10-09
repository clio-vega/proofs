"""Controls for the Q380 census (brief Sections 4(b),(e)).

(1) The statistic must VARY with the claim.  Plant a mutation that makes a
    nontrivial fibre monomial and confirm the census detector fires.
(2) The entrywise SSYT check must fire on a corrupted tableau.
(3) Negative control: dominance failure -> empty SSYT.
"""
import kf
from kf import K, partitions, dominates, nstat, ssyt, is_ssyt, charge, reading_word

print("--- control 1: does the monomial detector vary with charge? ---")
lam, mu = (3,1), (1,1,1,1)
d, c = K(lam, mu)
print(f"  true   K_{{{lam},{mu}}} = {d}  (#SSYT={c}, #nonzero coeffs={len(d)})")
# PLANT: collapse charge to a constant on this fibre only
orig = kf.charge
kf.charge = lambda w: 7
d2, c2 = K(lam, mu)
kf.charge = orig
print(f"  planted (charge := const)      -> {d2}  #nonzero coeffs={len(d2)}")
print(f"  detector reads {len(d)} -> {len(d2)}: VARIES = {len(d)!=len(d2)}")
assert len(d) > 1 and len(d2) == 1 and c2 >= 2, "detector is constant - instrument dead"
print("  => the planted pair WOULD have entered the target set; the real one did not.")

print()
print("--- control 2: does the entrywise SSYT check fire? ---")
T = ssyt((3,1),(2,1,1))[0]
print("  genuine tableau", T, "-> is_ssyt =", is_ssyt(T,(3,1),(2,1,1)))
bad = (tuple(sorted(T[0], reverse=True)), T[1]) if len(set(T[0]))>1 else ((9,9,9),T[1])
print("  row-order corrupted", bad, "-> is_ssyt =", is_ssyt(bad,(3,1),(2,1,1)))
bad2 = ((1,1,2),(1,))
print("  column corrupted  ", bad2, "-> is_ssyt =", is_ssyt(bad2,(3,1),(2,1,1)))
bad3 = ((1,1,1),(2,))
print("  content corrupted ", bad3, "-> is_ssyt =", is_ssyt(bad3,(3,1),(2,1,1)))
assert is_ssyt(T,(3,1),(2,1,1)) and not is_ssyt(bad,(3,1),(2,1,1))
assert not is_ssyt(bad2,(3,1),(2,1,1)) and not is_ssyt(bad3,(3,1),(2,1,1))
print("  => entrywise check fires on all three corruptions.")

print()
print("--- control 3: negative control, dominance failure ---")
n_off = n_empty = 0
for N in range(1,8):
    for l in partitions(N):
        for m in partitions(N):
            if not dominates(l,m):
                n_off += 1
                if len(ssyt(l,m))==0: n_empty += 1
print(f"  {n_off} pairs off dominance, {n_empty} with SSYT empty  (must be equal)")
assert n_off == n_empty

print()
print("--- control 4: cardinality vs content (brief 4(e)) ---")
tot=0; dup=0
for N in range(1,9):
    for l in partitions(N):
        for m in partitions(N):
            Ts = ssyt(l,m)
            tot += len(Ts)
            if len(set(Ts)) != len(Ts): dup += 1
            for T in Ts:
                assert is_ssyt(T,l,m), (l,m,T)
print(f"  {tot} tableaux built for N<=8; all pass entrywise shape/content/row/col; duplicate sets: {dup}")

print()
print("--- control 5: two-row vacuity claim from the brief ---")
import sympy
bad=0; n=0
for N in range(2,13):
    for j in range(0, N//2+1):
        for k in range(0, N//2+1):
            lam=(N-j,j) if j>0 else (N,)
            mu=(N-k,k) if k>0 else (N,)
            if not dominates(lam,mu): continue
            d,c = K(lam,mu); n+=1
            if not (c==1 and d=={k-j:1}):
                bad+=1; print("   FAIL",lam,mu,d,c)
print(f"  K_({{N-j,j}}),({{N-k,k}}) = t^(k-j) with #SSYT=1 : {n} pairs, {bad} failures")
print("  => the two-row family is a SINGLETON fibre throughout: VACUOUS for Q380. Confirmed.")
