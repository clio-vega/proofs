"""Q254: the Fink-Meszaros-St.Dizier (FMS) side of the comparison.

All definitions transcribed at first hand from
  projects/sources/1706.04935/forArxiv.tex
with line numbers:

  l.94   "By a diagram, we mean a subset D subset [n]^2 ... we view D as an
          ordered list of subsets D=(D_1,...,D_n) where D_j={i:(i,j) in D}"
  l.193  "Fix positive integers 1<=s_1<...<s_r<=n.  The sets {a_1,...,a_r}
          with a_1<...<a_r such that a_1<=s_1,...,a_r<=s_r are the bases of a
          matroid, called the Schubert matroid SM_n(s_1,...,s_r)"
  l.212  "For R,S subset [n], we say R<=S if #R=#S and the k-th least element
          of R does not exceed the k-th least element of S for each k"
  l.263-269 (inside the proof of Theorem 7): the support of chi_D is exactly
          {xi^C : C<=D},  xi^C_i = #{j : i in C_j},  and {S : S<=D_j} is
          exactly the basis set of SM_n(D_j).

So supp(chi_D) = sum over columns j of the indicator vectors of the bases of
SM_n(D_j).  That is what this module computes; it never needs the flagged Weyl
module itself.
"""
from itertools import combinations


def sm_bases(Dj, n):
    """Bases of the Schubert matroid SM_n(D_j), as sorted tuples.

    Dj = (s_1<...<s_r).  Bases: {a_1<...<a_r} subset [n] with a_k <= s_k.
    (FMS l.193.)  Equivalently {S subset [n] : S <= D_j} in FMS's order (l.212).
    """
    s = tuple(sorted(Dj))
    r = len(s)
    return [B for B in combinations(range(1, n + 1), r)
            if all(B[k] <= s[k] for k in range(r))]


def indicator(B, n):
    v = [0] * n
    for b in B:
        v[b - 1] = 1
    return tuple(v)


def supp_chi(D, n):
    """supp(chi_D) subset N^n for D = list of columns (subsets of [n])."""
    cur = {tuple([0] * n)}
    for Dj in D:
        if not Dj:
            continue
        vs = [indicator(B, n) for B in sm_bases(Dj, n)]
        cur = {tuple(a + b for a, b in zip(u, v)) for u in cur for v in vs}
    return cur


def active_rows(D):
    """m = max of the union of the columns (0 if D is empty)."""
    u = set()
    for Dj in D:
        u |= set(Dj)
    return max(u) if u else 0


def is_Sm_stable(S, m):
    """Is the support stable under permuting coordinates 1..m?

    Enough to test the adjacent transpositions (i,i+1), i<m: they generate S_m.
    """
    for i in range(m - 1):
        for a in S:
            b = list(a)
            b[i], b[i + 1] = b[i + 1], b[i]
            if tuple(b) not in S:
                return False
    return True


def is_bottom_justified(Dj, m):
    """D_j = {m-|D_j|+1, ..., m}?  (empty column counts as bottom-justified)"""
    if not Dj:
        return True
    k = len(Dj)
    return set(Dj) == set(range(m - k + 1, m + 1))


def lam_hat_from_columns(D):
    """If every column is bottom-justified, supp(chi_D) = {alpha: sort(alpha) <= lamhat}
    with lamhat' = sorted column sizes.  Return lamhat."""
    sizes = sorted((len(Dj) for Dj in D if Dj), reverse=True)
    if not sizes:
        return ()
    out = []
    for i in range(1, sizes[0] + 1):
        out.append(sum(1 for s in sizes if s >= i))
    # out is the conjugate of `sizes`; lamhat = conjugate of the column-size partition
    return tuple(out)
