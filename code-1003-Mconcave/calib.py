"""Calibration: slice-sum formula vs independent chain enumeration, m=2,3,4.
Also the symmetry controls k(a,b)=k(b,a) and k(a,b)=k(d-a-b,b)."""
import gen, sys

def run(m, nmax, dmax):
    agree = dis = 0
    sym_ok = sym_bad = 0
    sym3_ok = sym3_bad = 0
    bad = []
    for n in range(m, nmax + 1):
        for (mu, lam) in gen.pairs(n, m, dmax):
            d = sum(lam) - sum(mu)
            K = gen.k_table_chains(mu, lam, n, m)
            for b in range(0, d + 1):
                off, co = gen.slice_sum(mu, lam, n, m, b)
                mine = {off + j: c for j, c in enumerate(co) if c}
                theirs = {a: v for (a, bb), v in K.items() if bb == b and v}
                if mine == theirs:
                    agree += 1
                else:
                    dis += 1
                    if len(bad) < 5:
                        bad.append((n, mu, lam, b, mine, theirs))
            # symmetries, from the INDEPENDENT table
            for (a, b), v in K.items():
                if K.get((b, a), 0) == v: sym_ok += 1
                else: sym_bad += 1
                if K.get((d - a - b, b), 0) == v: sym3_ok += 1
                else: sym3_bad += 1
    return agree, dis, bad, (sym_ok, sym_bad), (sym3_ok, sym3_bad)

for (m, nmax, dmax) in [(2, 7, 9), (3, 7, 8), (4, 7, 7)]:
    a, dd, bad, s2, s3 = run(m, nmax, dmax)
    print(f"m={m} n<={nmax} d<={dmax}: slices agree {a}, disagree {dd}")
    print(f"   k(a,b)=k(b,a): {s2[0]} ok / {s2[1]} bad ;  k(a,b)=k(d-a-b,b): {s3[0]} ok / {s3[1]} bad")
    for x in bad: print("   BAD", x)
