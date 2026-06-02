"""
TASK B: d=3 character-side prediction table (Littlewood 3-core/3-quotient).
(i)  distribution of |3-core| and of floor(|3-core|/3) across shapes, n=3..10.
(ii) for 3-core-empty shapes: f^q0 f^q1 f^q2 * multinom(|q0|,|q1|,|q2|),
     cross-checked against #{standard 3-ribbon tableaux} = chi^lam(mu_3) up to sign,
     where mu_3 = class with all parts =3 (richest in 3's), n divisible by 3.
"""
from dcore_lib import (partitions, d_core_and_quotient, num_syt, n_of_lambda,
                       multinom, mn_char)
from collections import Counter

print("== TASK B (i): distribution of |3-core| and floor(|3-core|/3), n=3..10 ==")
coresize_dist = Counter()
ordexp_dist = Counter()
per_n = {}
for n in range(3, 11):
    cs = Counter(); oe = Counter()
    for lam in partitions(n):
        core, _ = d_core_and_quotient(lam, 3)
        s = sum(core)
        cs[s] += 1; oe[s // 3] += 1
        coresize_dist[s] += 1; ordexp_dist[s // 3] += 1
    per_n[n] = (cs, oe)

print("  |3-core| sizes observed (overall):", dict(sorted(coresize_dist.items())))
print("  floor(|3-core|/3) observed (overall):", dict(sorted(ordexp_dist.items())))
print("  (3-core sizes are 0,1,2,4,5,7,... i.e. NOT 3; 3 cells can always be removed)")
print("\n  per-n |3-core| distribution:")
for n in range(3, 11):
    cs, oe = per_n[n]
    print(f"   n={n}: coresize {dict(sorted(cs.items()))}  | floor/3 {dict(sorted(oe.items()))}")

print("\n== TASK B (ii): 3-core-empty leading-structure (only n divisible by 3) ==")
print("   pred = f^q0 f^q1 f^q2 * multinom(|q0|,|q1|,|q2|)")
print("   cross-check |chi^lam(3,3,...,3)| (= #standard 3-ribbon tableaux)")
fails = []
count = 0
samples = []
for n in range(3, 11):
    if n % 3 != 0:
        # empty 3-core requires n divisible by 3
        for lam in partitions(n):
            core, quo = d_core_and_quotient(lam, 3)
            if sum(core) == 0:
                fails.append(("UNEXPECTED-EMPTY", lam))
        continue
    k = n // 3
    mu3 = tuple([3] * k)  # all-3 cycle type
    for lam in partitions(n):
        core, (q0, q1, q2) = d_core_and_quotient(lam, 3)
        if sum(core) != 0:
            # chi on all-3 class must vanish if 3-core nonempty
            chi = mn_char(lam, mu3)
            if chi != 0:
                fails.append(("NONEMPTY-NONZERO-CHI", lam, chi))
            continue
        count += 1
        f0, f1, f2 = num_syt(q0), num_syt(q1), num_syt(q2)
        mult = multinom([len(q0) and sum(q0) or 0, sum(q1), sum(q2)])
        mult = multinom([sum(q0), sum(q1), sum(q2)])
        pred = f0 * f1 * f2 * mult
        chi = mn_char(lam, mu3)
        # cross-check: |chi| should equal pred (Littlewood); sign separate
        if abs(chi) != pred:
            fails.append(("MAG", lam, abs(chi), pred, (q0, q1, q2)))
        else:
            if len(samples) < 14:
                sgn = chi // pred if pred else 0
                samples.append((n, lam, (q0, q1, q2), pred, chi, sgn))
print(f"  3-core-empty shapes (n div by 3) checked: {count}")
print(f"  cross-check mismatches: {len(fails)}")
for f in fails[:20]:
    print("   ", f)
print(f"  VERDICT (ii): {'PASS' if not fails else 'FAIL'}")
print("\n  sample (n, lam, (q0,q1,q2), pred=f0f1f2*multinom, chi(3^k), sign):")
for s in samples:
    print("   ", s)
