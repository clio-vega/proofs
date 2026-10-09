"""Two exact reductions for the Kostka number, and the irreducible core of the
K_{lam,mu}=1 locus.

R1 (first row):  if mu_1 = lam_1 then every 1 fills row 1 completely, so
     K_{lam,mu} = K_{(lam_2,lam_3,...),(mu_2,mu_3,...)}.
R2 (first column): if l(mu) = l(lam) = m then column 1 is forced to be
     (1,2,...,m), so K_{lam,mu} = K_{lam - col_1, mu - (1^m)}.
Note l(mu) >= l(lam) always when lam >= mu, so R2 is the opposite extreme to R1.
"""
from kf import K, partitions, dominates, ssyt


def nz(p):
    return tuple(x for x in p if x > 0)


def R1(lam, mu):
    return nz(lam[1:]), nz(mu[1:])


def R2(lam, mu):
    m = len(lam)
    return nz(tuple(x - 1 for x in lam)), nz(tuple(x - 1 for x in mu))


def kostka(lam, mu):
    return len(ssyt(lam, mu))


def irreducible(lam, mu):
    return len(lam) >= 2 and mu[0] < lam[0] and len(mu) > len(lam)


def core(lam, mu):
    """Reduce until irreducible or trivial."""
    while True:
        if len(lam) == 0:
            return lam, mu
        if len(lam) == 1:
            return lam, mu
        if mu[0] == lam[0]:
            lam, mu = R1(lam, mu); continue
        if len(mu) == len(lam):
            lam, mu = R2(lam, mu); continue
        return lam, mu


if __name__ == "__main__":
    # --- verify R1 and R2 are exact
    b1 = b2 = n1 = n2 = 0
    for N in range(1, 10):
        for lam in partitions(N):
            for mu in partitions(N):
                if not dominates(lam, mu):
                    continue
                k = kostka(lam, mu)
                if len(lam) >= 1 and mu[0] == lam[0]:
                    n1 += 1
                    l2, m2 = R1(lam, mu)
                    if kostka(l2, m2) != k:
                        b1 += 1; print("R1 FAIL", lam, mu, k, kostka(l2, m2))
                if len(mu) == len(lam):
                    n2 += 1
                    l2, m2 = R2(lam, mu)
                    if kostka(l2, m2) != k:
                        b2 += 1; print("R2 FAIL", lam, mu, k, kostka(l2, m2))
    print(f"R1 checked on {n1} pairs (N<=9): {b1} failures")
    print(f"R2 checked on {n2} pairs (N<=9): {b2} failures")

    # --- irreducible core of the K=1 locus
    print()
    print("IRREDUCIBLE pairs (l(lam)>=2, mu_1<lam_1, l(mu)>l(lam)) with K=1, N<=10:")
    cnt_irr = 0; cnt_irr1 = 0
    found = []
    for N in range(1, 11):
        for lam in partitions(N):
            for mu in partitions(N):
                if not dominates(lam, mu) or not irreducible(lam, mu):
                    continue
                cnt_irr += 1
                if kostka(lam, mu) == 1:
                    cnt_irr1 += 1
                    found.append((N, lam, mu))
    print(f"  {cnt_irr} irreducible dominance pairs for N<=10; {cnt_irr1} of them have K=1")
    for f in found:
        print("   ", f)


def R3(lam, mu):
    """Rectangular complementation (GL_r duality): with r = max(l(lam),l(mu)),
    c = lam_1, send (lam,mu) to (lam*, mu*), lam* = (c-lam_r,...,c-lam_1)."""
    r = max(len(lam), len(mu)); c = lam[0]
    L = tuple(lam) + (0,)*(r-len(lam)); M = tuple(mu) + (0,)*(r-len(mu))
    return nz(tuple(c-L[r-1-i] for i in range(r))), nz(tuple(c-M[r-1-i] for i in range(r)))
