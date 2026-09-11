"""
Fresh engine for Q147: ribbon operators with FREE height weights.

Nothing here hard-codes (-t)^N or t^ht.  The weight of a ribbon of height N is
an opaque symbol w[N] (shared across e) or w[e][N] (independent per e).

Two independent implementations of "add a connected e-ribbon":
  * ribbons_young : brute force over partitions mu of |lam|+e, test border strip
  * ribbons_maya  : bead move b -> b+e on the beta-set, height = beads jumped
They are cross-checked against each other.
"""
from itertools import combinations
import sympy as sp


# ---------- partitions ----------------------------------------------------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def contains(mu, lam):
    if len(lam) > len(mu):
        return False
    return all(mu[i] >= lam[i] for i in range(len(lam)))


def cells(lam):
    return {(i, j) for i, p in enumerate(lam) for j in range(p)}


# ---------- implementation 1: Young diagram ------------------------------
def ribbons_young(lam, e):
    """returns list of (mu, height) for connected e-ribbons added to lam."""
    out = []
    n = sum(lam) + e
    L = cells(lam)
    for mu in partitions(n):
        if not contains(mu, lam):
            continue
        S = cells(mu) - L
        if len(S) != e:
            continue
        # no 2x2 block
        if any((i, j) in S and (i + 1, j) in S and (i, j + 1) in S and (i + 1, j + 1) in S
               for (i, j) in S):
            continue
        # connected (edge-adjacency)
        start = next(iter(S))
        seen = {start}
        stack = [start]
        while stack:
            (i, j) = stack.pop()
            for (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
                if (a, b) in S and (a, b) not in seen:
                    seen.add((a, b)); stack.append((a, b))
        if len(seen) != e:
            continue
        rows = {i for (i, j) in S}
        out.append((mu, len(rows) - 1))
    return out


# ---------- implementation 2: Maya / beta-set ----------------------------
def maya(lam, L):
    """beads {lam_j - j : j>=1} intersected with [-L, +inf); returns frozenset."""
    ext = list(lam) + [0] * (L + 2)
    return frozenset(ext[j - 1] - j for j in range(1, L + 2) if ext[j - 1] - j >= -L)


def unmaya(M, L):
    beads = sorted([x for x in M], reverse=True)
    lam = []
    for j, b in enumerate(beads, start=1):
        p = b + j
        if p > 0:
            lam.append(p)
        else:
            break
    return tuple(lam)


def ribbons_maya(lam, e, L=None):
    if L is None:
        L = sum(lam) + e + 4
    M = maya(lam, L)
    out = []
    for b in sorted(M):
        if (b + e) in M:
            continue
        M2 = (M - {b}) | {b + e}
        ht = len([j for j in M if b < j < b + e])
        out.append((unmaya(M2, L), ht))
    return out


def _crosscheck(maxn=8, maxe=5):
    for n in range(0, maxn + 1):
        for lam in partitions(n):
            for e in range(1, maxe + 1):
                a = sorted(ribbons_young(lam, e))
                b = sorted(ribbons_maya(lam, e))
                assert a == b, (lam, e, a, b)
    return True


# ---------- the operator with free weights -------------------------------
def apply_R(vec, e, w):
    """vec: dict partition -> coefficient.  w: callable (e, N) -> weight."""
    out = {}
    for lam, c in vec.items():
        for mu, ht in ribbons_maya(lam, e):
            out[mu] = sp.expand(out.get(mu, 0) + c * w(e, ht))
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def commutator(lam, e, f, w):
    """[R_e, R_f] s_lam  as dict mu -> polynomial."""
    v = {tuple(lam): sp.Integer(1)}
    a = apply_R(apply_R(v, f, w), e, w)   # R_e R_f
    b = apply_R(apply_R(v, e, w), f, w)   # R_f R_e
    keys = set(a) | set(b)
    out = {}
    for k in keys:
        val = sp.expand(a.get(k, 0) - b.get(k, 0))
        if val != 0:
            out[k] = val
    return out


if __name__ == "__main__":
    print("cross-check young vs maya:", _crosscheck())
