"""STEP 0: census of  C(nu) := sum_i (L_i(nu)+R_i(nu))  on each slice.
Report counts, not a verdict.  Also record the observed value of C(nu) itself."""
import gen, sys
from collections import Counter

def census(m, nmax, dmax):
    stats = Counter()
    valforms = Counter()
    witnesses = []
    for n in range(m, nmax + 1):
        for (mu, lam) in gen.pairs(n, m, dmax):
            d = sum(lam) - sum(mu)
            for b in range(0, d + 1):
                nus = gen.slice_nus(mu, lam, n, m, b)
                if not nus:
                    continue
                stats['slices'] += 1
                vals = set()
                for nu in nus:
                    C = sum(L + R for (L, R) in gen.LR(nu, lam, n, m))
                    vals.add(C)
                    # the predicted value: S_nu + |lam|
                    valforms[C - sum(nu) - sum(lam)] += 1
                stats['summands'] += len(nus)
                if len(vals) == 1:
                    stats['const'] += 1
                else:
                    stats['nonconst'] += 1
                    if len(witnesses) < 5:
                        witnesses.append((n, mu, lam, b, sorted(vals)))
    return stats, valforms, witnesses

for (m, nmax, dmax) in [(2, 16, 20), (3, 14, 18), (4, 14, 18), (5, 12, 14), (6, 11, 12)]:
    st, vf, w = census(m, nmax, dmax)
    print(f"m={m}, {m}<=n<={nmax}, d<={dmax}: {st['slices']} nonempty slices, {st['summands']} summands")
    print(f"   C constant on slice: {st['const']}   non-constant: {st['nonconst']}")
    print(f"   C(nu) - S_nu - |lam|  distribution: {dict(vf)}")
    for x in w: print("   WITNESS", x)
    sys.stdout.flush()
