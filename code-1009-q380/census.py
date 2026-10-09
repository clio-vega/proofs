"""Q380 census.  For each N, over all pairs (lam, mu) of partitions of N:
counts, and the four populations the brief demands."""
import sys, json
from kf import K, partitions, dominates, nstat

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
rows = []
k1 = []          # pairs with #SSYT == 1
mono_nontriv = []  # pairs with #SSYT >= 2 and K a monomial   <-- the Q380 target set
monic_fail = []    # pairs where leading coeff != 1  (control on (6.5)(ii))
deg_fail = []      # pairs where deg != n(mu) - n(lam)
for N in range(1, NMAX + 1):
    ps = list(partitions(N))
    tot = nonempty = ge2 = mono = 0
    for lam in ps:
        for mu in ps:
            tot += 1
            if not dominates(lam, mu):
                # control: dominance failure must give empty SSYT
                d, c = K(lam, mu)
                assert c == 0, ("nonempty off dominance", lam, mu)
                continue
            d, c = K(lam, mu)
            assert c >= 1
            assert sum(d.values()) == c
            nonempty += 1
            if max(d) != nstat(mu) - nstat(lam):
                deg_fail.append((lam, mu))
            if d[max(d)] != 1:
                monic_fail.append((lam, mu, d))
            if c == 1:
                k1.append((N, lam, mu))
            if c >= 2:
                ge2 += 1
            if len(d) == 1:
                mono += 1
                if c >= 2:
                    mono_nontriv.append((lam, mu, d))
    rows.append((N, len(ps), tot, nonempty, ge2, mono))
    print(f"N={N:2d}  #parts={len(ps):4d}  pairs={tot:7d}  SSYT!=0={nonempty:7d}  "
          f"#SSYT>=2={ge2:7d}  K monomial={mono:7d}", flush=True)

print()
print("deg != n(mu)-n(lam) :", len(deg_fail))
print("leading coeff != 1  :", len(monic_fail))
print("### Q380 TARGET SET: #SSYT>=2 AND K monomial :", len(mono_nontriv))
if mono_nontriv:
    print(mono_nontriv[:20])
print()
print("#SSYT==1 pairs total:", len(k1))
json.dump({"k1": [[n, list(l), list(m)] for n, l, m in k1]}, open("k1.json", "w"))
