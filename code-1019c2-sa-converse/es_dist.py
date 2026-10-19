r"""How many terms does the residue sum actually have?

If |E_S| = 1 for every pole in the census, the census carries ZERO bits about
cancellation -- the measurement would be a plant drawn from inside its own
hypothesis class.  So: the distribution of |E_S| over the triples where the
predicate fires.  DOMAIN: all triples (lam,mu,nu) in (2^[N])^3, N = 2..6, all m.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'code-1017c2-coboundary'))
import q484_engine as E

for N in range(2, 7):
    P = E.subsets(N)
    dist = {}
    for m in range(1, N + 1):
        for lam in P:
            for mu in P:
                for nu in P:
                    c = 0
                    for R in P:
                        if m in R: continue
                        if E.geq(nu, tuple(sorted(R + (m,)))) and E.geq(R, lam) and E.geq(R, mu):
                            c += 1
                    if c:
                        dist[c] = dist.get(c, 0) + 1
    tot = sum(dist.values())
    multi = sum(v for k, v in dist.items() if k >= 2)
    print(f"N={N}: {tot} (m,lam,mu,nu) with E_S != {{}};  |E_S| distribution "
          f"{dict(sorted(dist.items()))};  |E_S|>=2 in {multi} ({100*multi/tot:.1f}%)")
