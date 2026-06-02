"""
TASK A: validate t=2 <-> Albion bridge on G_lambda data, n=2..10.
(i)  c(lam)=0 (2-core empty) <=> G_lam(-1) != 0
(ii) for c=0 shapes, G_lam(-1) == (-1)^{n(lam)} f^{q0} f^{q1} binom(w,|q0|)
     and check sign equals (-1)^{n(lam)}.
"""
from dcore_lib import (partitions, d_core_and_quotient, G_minus1,
                       num_syt, n_of_lambda)
from sympy import binomial

print("== TASK A (i): c(lam)=floor(|2-core|/2)=0  <=>  G_lam(-1) != 0 ==")
print("   (NOTE: c=0 means |2-core| in {0,1}: empty for even n, single cell for odd n)")
exceptions_i = []
n_shapes = 0
for n in range(2, 11):
    for lam in partitions(n):
        n_shapes += 1
        core, _ = d_core_and_quotient(lam, 2)
        c = sum(core) // 2   # floor(|2-core|/2)
        val = G_minus1(lam)
        if (c == 0) != (val != 0):
            exceptions_i.append((lam, sum(core), c, val))
print(f"  shapes checked: {n_shapes}")
print(f"  exceptions to (c=0 <=> G(-1)!=0): {len(exceptions_i)}")
for e in exceptions_i:
    print("   ", e)
print(f"  VERDICT (i): {'PASS' if not exceptions_i else 'FAIL'}")

print("\n== TASK A (ii): leading value (S(n,0)=1) for ALL c=0 shapes (both parities) ==")
# Even n (empty core): G_lam(-1) =?= (-1)^{n(lam)} f^q0 f^q1 binom(w,|q0|)
# Odd  n (single-cell core, f^core=1): same formula, w=(n-|core|)/2.
# Independent magnitude = f^q0 f^q1 binom(w,|q0|) (* f^core, which is 1 here).
fails_ii = []
signfails = []
count = 0
count_even = count_odd = 0
for n in range(2, 11):
    for lam in partitions(n):
        core, (q0, q1) = d_core_and_quotient(lam, 2)
        c = sum(core) // 2
        if c != 0:
            continue
        count += 1
        if n % 2 == 0:
            count_even += 1
        else:
            count_odd += 1
        val = G_minus1(lam)
        w = (n - sum(core)) // 2
        f0 = num_syt(q0); f1 = num_syt(q1); fcore = num_syt(core)
        base = f0 * f1 * int(binomial(w, sum(q0))) * fcore
        pred = ((-1) ** n_of_lambda(lam)) * base
        if val != pred:
            fails_ii.append((lam, val, pred))
        # independent magnitude / sign check
        if abs(val) != base:
            signfails.append(("MAG", lam, abs(val), base))
        else:
            sign = val // base if base != 0 else 0
            exp_sign = (-1) ** n_of_lambda(lam)
            if sign != exp_sign:
                signfails.append(("SIGN", lam, sign, exp_sign))
print(f"  c=0 shapes checked: {count}  (even-n empty-core: {count_even}, odd-n single-core: {count_odd})")
print(f"  full-formula mismatches: {len(fails_ii)}")
for f in fails_ii:
    print("   ", f)
print(f"  independent magnitude(=f0 f1 C(w,|q0|) f^core)/sign mismatches: {len(signfails)}")
for f in signfails:
    print("   ", f)
print(f"  sign always equals (-1)^n(lam): {'YES' if not signfails else 'NO'}")
print(f"  VERDICT (ii): {'PASS' if not fails_ii and not signfails else 'FAIL'}")

# small illustrative table
print("\n  sample (n, lam, q0, q1, |G(-1)|=f0 f1 C(w,|q0|), sign==(-1)^n(lam)):")
shown = 0
for n in [4, 6, 8]:
    for lam in partitions(n):
        core, (q0, q1) = d_core_and_quotient(lam, 2)
        if sum(core) != 0:
            continue
        val = G_minus1(lam); w = n // 2
        base = num_syt(q0) * num_syt(q1) * int(binomial(w, sum(q0)))
        if shown < 12:
            print(f"   n={n} lam={lam} q0={q0} q1={q1} base={base} G(-1)={val}")
            shown += 1
