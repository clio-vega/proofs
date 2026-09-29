"""B-EXC (Murota's SYMMETRIC exchange axiom) for
    J(lhat) = { alpha in N^ell : sort(alpha) dominated by lhat }.

Written from the DEFINITION, independently of the Lean files.  Three questions, all
armed by state/LEAN.md 2026-09-29 c1 section 3:

 (EX2)  B-EXC itself, widened past the 876317-triple range of
        mconvex_exchange_check.py:  for all alpha,beta in J and i with beta_i < alpha_i,
        does SOME j satisfy  alpha_j < beta_j  and  alpha - e_i + e_j in J  and
        beta + e_i - e_j in J ?

 (AUT)  Is the symmetric conjunct AUTOMATIC?  i.e. does EVERY j that works one-sidedly
        also work for beta?  If yes, today's theorem has no content beyond
        insupp_exchange.  If no, the witness is the NEGATIVE CONTROL that LEAN.md
        section 5 demands, and the construction must change.

 (NB)   Is the proposed proof's witness set the right one?  The derivation says the
        good j are exactly those in  B \\ A, where
            A = union of all alpha-tight S with i notin S      (maximal, alpha-tight)
            B = intersection of all beta-tight S with i in S   (minimal, beta-tight)
        and "tight" means alpha(S) = Lambda_{|S|}.  Check that
            { j : alpha_j < beta_j, both surgeries land in J }  ==  D cap (B \\ A)
        exactly -- not merely that it is nonempty.  A construction validated only by
        nonemptiness would pass with the wrong witness set.
"""
import itertools, sys

sys.setrecursionlimit(10000)


def partitions(d, maxlen, maxpart=None):
    if maxpart is None:
        maxpart = d
    if d == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for first in range(min(d, maxpart), 0, -1):
        for rest in partitions(d - first, maxlen - 1, first):
            yield (first,) + rest


def weak_compositions(d, ell):
    if ell == 0:
        if d == 0:
            yield ()
        return
    for first in range(d + 1):
        for rest in weak_compositions(d - first, ell - 1):
            yield (first,) + rest


def pad(nu, ell):
    return tuple(nu) + (0,) * (ell - len(nu))


def dominated(sigma, lhat, ell):
    s = t = 0
    for r in range(ell):
        s += sigma[r]
        t += lhat[r]
        if s > t:
            return False
    return True


def J_sorted(lhat, ell):
    """The PAPER's object: sort(alpha) dominated by lhat.  Not the subset form."""
    d = sum(lhat)
    L = pad(lhat, ell)
    return {a for a in weak_compositions(d, ell)
            if dominated(tuple(sorted(a, reverse=True)), L, ell)}


def step(alpha, i, j):
    """alpha - e_i + e_j."""
    a = list(alpha)
    a[i] -= 1
    a[j] += 1
    return tuple(a)


def lam(lhat, ell):
    L = pad(lhat, ell)
    Lam = [0] * (ell + 1)
    for r in range(ell):
        Lam[r + 1] = Lam[r] + L[r]
    return Lam


def tight_sets(alpha, Lam, ell):
    """All S subset [ell] with alpha(S) = Lambda_{|S|}, as frozensets."""
    out = []
    for k in range(ell + 1):
        for S in itertools.combinations(range(ell), k):
            if sum(alpha[c] for c in S) == Lam[k]:
                out.append(frozenset(S))
    return out


def setA(alpha, Lam, ell, i):
    """Union of all alpha-tight S avoiding i."""
    A = frozenset()
    for S in tight_sets(alpha, Lam, ell):
        if i not in S:
            A = A | S
    return A


def setB(beta, Lam, ell, i):
    """Intersection of all beta-tight S containing i.  ([ell] is always tight.)"""
    B = frozenset(range(ell))
    for S in tight_sets(beta, Lam, ell):
        if i in S:
            B = B & S
    return B


def run(dmax, ellmax, report_control=3):
    ex2_bad = ex2_tot = 0
    aut_bad = aut_tot = 0
    nb_bad = nb_tot = 0
    controls = []
    for ell in range(1, ellmax + 1):
        for d in range(0, dmax + 1):
            for lhat in partitions(d, ell):
                Jset = J_sorted(lhat, ell)
                Lam = lam(lhat, ell)
                for alpha in Jset:
                    Acache = {}
                    for beta in Jset:
                        for i in range(ell):
                            if not (beta[i] < alpha[i]):
                                continue
                            D = [j for j in range(ell) if alpha[j] < beta[j]]
                            one = [j for j in D if step(alpha, i, j) in Jset]
                            both = [j for j in one if step(beta, j, i) in Jset]
                            # (EX2)
                            ex2_tot += 1
                            if not both:
                                ex2_bad += 1
                                print("EX2 FAIL", lhat, ell, alpha, beta, i)
                            # (AUT): a one-sided j that fails symmetrically
                            for j in one:
                                aut_tot += 1
                                if j not in both:
                                    aut_bad += 1
                                    if len(controls) < report_control:
                                        controls.append((lhat, ell, alpha, beta, i, j))
                            # (NB): predicted witness set
                            if i not in Acache:
                                Acache[i] = setA(alpha, Lam, ell, i)
                            A = Acache[i]
                            B = setB(beta, Lam, ell, i)
                            pred = sorted(set(D) & (B - A))
                            nb_tot += 1
                            if pred != sorted(both):
                                nb_bad += 1
                                if nb_bad <= 5:
                                    print("NB FAIL", lhat, ell, alpha, beta, i,
                                          "pred", pred, "actual", sorted(both),
                                          "A", sorted(A), "B", sorted(B))
    print()
    print(f"range: |lhat| <= {dmax}, ell <= {ellmax}   (J built in the SORTED form)")
    print(f"(EX2) B-EXC, symmetric exchange : {ex2_bad} failures / {ex2_tot} triples (alpha,beta,i)")
    print(f"(AUT) one-sided j that also works for beta : "
          f"{aut_tot - aut_bad} / {aut_tot} quadruples (alpha,beta,i,j)")
    print(f"      => symmetric conjunct AUTOMATIC? {'YES' if aut_bad == 0 else 'NO'}"
          f"  ({aut_bad} one-sided j fail for beta)")
    print(f"(NB)  predicted witness set D cap (B\\A) == actual : "
          f"{nb_tot - nb_bad} / {nb_tot} triples, {nb_bad} mismatches")
    if controls:
        print()
        print("NEGATIVE CONTROL candidates (alpha_j<beta_j, alpha-e_i+e_j in J, beta+e_i-e_j NOT in J):")
        for (lhat, ell, alpha, beta, i, j) in controls:
            print(f"  lhat={lhat} ell={ell} alpha={alpha} beta={beta} i={i} j={j}")
            print(f"    alpha-e_i+e_j = {step(alpha,i,j)}  in J: {step(alpha,i,j) in J_sorted(lhat,ell)}")
            print(f"    beta +e_i-e_j = {step(beta,j,i)}  in J: {step(beta,j,i) in J_sorted(lhat,ell)}")
    return ex2_bad + nb_bad


if __name__ == "__main__":
    dmax = int(sys.argv[1]) if sys.argv[1:] else 11
    ellmax = int(sys.argv[2]) if sys.argv[2:] else 5
    sys.exit(1 if run(dmax, ellmax) else 0)
