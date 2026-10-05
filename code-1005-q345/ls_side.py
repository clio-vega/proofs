"""Side B: Lenart-Sottile's increasing labelled chains (math/0202090, Thm 2).

Built ONLY from the definitions in the paper:

  labelled Bruhat order: a directed edge  u --(k,b)--> w  for every Bruhat cover
  u <. w with u^{-1} w = (i,j), i < j, for each k with i <= k < j, and
  b = u(i) = w(j).   (j-i edges per cover.)

  a chain is INCREASING if its sequence of labels (k_1,b_1),...,(k_m,b_m) is
  strictly increasing in LEXICOGRAPHIC order on pairs.

  type alpha of a chain = the composition whose i-th part counts the
  occurrences of i as a FIRST coordinate.   x^gamma = x^alpha.

  I_alpha(u,w) = # increasing chains u -> w of type alpha.

NOTHING here calls monk_matrix or chains.T_tensor.  The only shared code is
`length` and `rmul_t` (the definition of the symmetric group), and the Schubert
table (independent ground truth, divided differences).
"""
from collections import Counter
from itertools import permutations
from schub import length, rmul_t


def labelled_edges(u, n):
    """all (k, b, w) with u --(k,b)--> w an edge of the labelled Bruhat order."""
    out = []
    lu = length(u)
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            w = rmul_t(u, i, j)
            if length(w) != lu + 1:
                continue
            b = u[i-1]                      # = u(i) = w(j)
            for k in range(i, j):           # i <= k < j
                out.append((k, b, w))
    return out


def increasing_chains(u, w, n):
    """yield every increasing labelled chain u -> w as a list of (k,b) labels."""
    target = length(w)
    def rec(v, last):
        if v == w:
            yield []
            return
        if length(v) >= target:
            return
        for (k, b, v2) in labelled_edges(v, n):
            if (k, b) <= last:              # strict lex increase required
                continue
            for rest in rec(v2, (k, b)):
                yield [(k, b)] + rest
    if length(u) > target:
        return
    yield from rec(u, (0, 0))               # (0,0) < every real label


def I_table(u, w, n):
    """Counter: type alpha (tuple of length n-1) -> I_alpha(u,w)."""
    out = Counter()
    for ch in increasing_chains(u, w, n):
        a = [0]*(n-1)
        for (k, b) in ch:
            a[k-1] += 1
        out[tuple(a)] += 1
    return out


def chain_lengths(u, w, n):
    """length distribution of the increasing chains -- a non-vacuity probe."""
    return Counter(len(ch) for ch in increasing_chains(u, w, n))
