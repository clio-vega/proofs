"""How much of the converse does Theorem D cover?

Theorem D needs nu MINIMAL in Motegi's order among
    S(a,b) := { kappa : kappa >= lam, kappa >= mu, a,b in kappa,
                kappa minus {a,b} >= lam,  kappa minus {a,b} >= mu }
for SOME pair a<b.  We minimise over the SUPPORT BOUND (kappa >= lam,mu), a
superset of the true support -- more competitors, so a CONSERVATIVE criterion.
"""
import sys
from itertools import combinations
from collections import Counter

def subsets(N):
    out = []
    for k in range(N+1):
        out.extend(combinations(range(1, N+1), k))
    return out

def geq(M, I):
    s, r = len(M), len(I)
    if s < r:
        return False
    return all(M[s-r+j] >= I[j] for j in range(r))

def main(N):
    S = subsets(N)
    tot = 0; cov = 0; unc = []
    for lam in S:
        for mu in S:
            if mu < lam:
                continue
            for nu in S:
                if not (geq(nu, lam) and geq(nu, mu)):
                    continue
                if len(nu) - max(len(lam), len(mu)) < 2:
                    continue
                tot += 1
                covered = False
                for a, b in combinations(nu, 2):
                    Sab = [k for k in S
                           if geq(k, lam) and geq(k, mu) and a in k and b in k
                           and geq(tuple(x for x in k if x not in (a, b)), lam)
                           and geq(tuple(x for x in k if x not in (a, b)), mu)]
                    if nu in Sab and not any(k != nu and geq(nu, k) for k in Sab):
                        covered = True
                        break
                if covered:
                    cov += 1
                else:
                    unc.append((lam, mu, nu))
    print(f'N={N}: jump>=2 triples {tot};  Theorem D covers {cov}'
          f'  ({100*cov/tot:.1f}%);  uncovered {len(unc)}')
    if unc:
        c = Counter((len(l), len(m), len(n)) for l, m, n in unc)
        print(f'   uncovered profile (|lam|,|mu|,|nu|): {dict(sorted(c.items()))}')
        print(f'   all uncovered have lam=mu=empty?  '
              f'{all(l == () and m == () for l, m, _ in unc)}')
        print(f'   examples: {unc[:4]}')

if __name__ == '__main__':
    for N in range(2, 7):
        main(N)
