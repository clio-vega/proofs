"""Core affine-symmetric-group machinery for Q220.

Conventions (fixed once, both papers reconciled):
  * residues live in Z/n, written 0..n-1; a subset S is identified with its
    n-periodic lift to Z.
  * u_S is the unique cyclically decreasing element with support S
    (Clio 2026-09-21, Lemma 2.2):  u_S(j) = j-1 if j-1 in S,
    else u_S(j) = min{ c >= j : c not in S }.
  * products act as (ab)(j) = a(b(j)), so u_S u_T means "apply u_T first".
  * Morse-Schilling's content(u) is exactly S.
"""

from itertools import combinations


def u_S(S, n):
    """Return the map j -> u_S(j) on Z, for S a frozenset of residues mod n."""
    Sl = set(S)
    def f(j):
        if (j - 1) % n in Sl:
            return j - 1
        c = j
        while c % n in Sl:
            c += 1
        return c
    return f


def window(f, n):
    return tuple(f(i) for i in range(1, n + 1))


def compose(f, g):
    return lambda j: f(g(j))


def length(win, n):
    """Shi's length formula for an affine permutation given by its window."""
    tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            tot += abs((win[j] - win[i]) // n)
    return tot


def is_affine(win, n):
    return sum(win) == n * (n + 1) // 2 and len(set(x % n for x in win)) == n


def runs(S, n):
    """Maximal cyclic runs of S, as lists of residues [m, m+1, ..., M] (mod n).

    Returns a list of lists.  If S is all of Z/n there are no runs (excluded:
    S is always a proper subset here).
    """
    Sl = set(S)
    assert len(Sl) < n
    out = []
    for m in sorted(Sl):
        if (m - 1) % n not in Sl:            # m is the bottom of a run
            r = [m]
            j = (m + 1) % n
            while j in Sl:
                r.append(j)
                j = (j + 1) % n
            out.append(r)
    return out


# ---------------------------------------------------------------- Clio Thm 4.1

def clio_letters(S, T, n):
    """The set E(S,T) of letters e that Theorem 4.1 removes from S.

    e in E iff, with m = bottom of the S-run containing e:
        m not in T,  m+1,...,e all in T,  e+1 not in T.
    (This is exactly 'the run of e is usable and e is its first index with
    e+1 not in T'.)  One letter per usable run.
    """
    Tl = set(T)
    out = set()
    for r in runs(S, n):
        m = r[0]
        if m in Tl:
            continue
        for e in r:
            if (e + 1) % n not in Tl:
                out.add(e)
                break
    return out


def clio_moves(S, T, n):
    """Theorem 4.1 moves as pairs (S', T')."""
    mv = set()
    for r in runs(S, n):
        m = r[0]
        if m in set(T):
            continue
        for e in r:
            if (e + 1) % n not in set(T):
                mv.add((frozenset(S) - {e}, frozenset(T) | {m}))
                break
    return mv


# ------------------------------------------------------- Morse-Schilling e~_1

def ms_order(x, n):
    """The order (7) of MS: x-1 > x-2 > ... > 0 > n-1 > ... > x+1.

    Returned INCREASING, i.e. [x+1, x+2, ..., x-1] read cyclically.
    """
    return [(x + 1 + i) % n for i in range(n - 1)]


def ms_pairing(S, T, x, n):
    """Return (L1, R1): unpaired elements of S resp. T in the uv-pairing wrt x.

    Implemented directly from MS's prescription: process b in content(u)
    in DECREASING order wrt (7); pair b with the smallest a > b in content(v)
    not yet used.
    """
    order = ms_order(x, n)                    # increasing
    pos = {r: i for i, r in enumerate(order)}
    assert x not in S and x not in T
    Ss = sorted(S, key=lambda r: pos[r], reverse=True)   # decreasing
    Tfree = sorted(T, key=lambda r: pos[r])              # increasing
    used = set()
    L1 = []
    pairs = {}
    for b in Ss:
        cand = [a for a in Tfree if pos[a] > pos[b] and a not in used]
        if cand:
            a = cand[0]                       # smallest a > b
            used.add(a)
            pairs[b] = a
        else:
            L1.append(b)
    R1 = [a for a in T if a not in used]
    return set(L1), set(R1), pairs


def ms_etilde(S, T, x, n):
    """MS e~_1 at parameter x.  Returns (b, S', T') or None if L1 is empty."""
    L1, R1, _ = ms_pairing(S, T, x, n)
    if not L1:
        return None
    pos = {r: i for i, r in enumerate(ms_order(x, n))}
    b = min(L1, key=lambda r: pos[r])         # min wrt (7)
    t = 0
    while (b - t - 1) % n in set(S):
        t += 1
    m = (b - t) % n
    return b, frozenset(S) - {b}, frozenset(T) | {m}


def X_set(S, T, n):
    return [x for x in range(n) if x not in S and x not in T]


# ------------------------------------------------------------------ enumeration

def additive_pairs(n):
    """All (S,T), S,T proper subsets of Z/n, with l(u_S u_T) = |S| + |T|.

    Yields (S, T, window of v).
    """
    subs = []
    for k in range(0, n):
        for c in combinations(range(n), k):
            subs.append(frozenset(c))
    wins = {}
    for S in subs:
        wins[S] = window(u_S(S, n), n)
    for S in subs:
        fS = u_S(S, n)
        for T in subs:
            fT = u_S(T, n)
            w = window(compose(fS, fT), n)
            if length(w, n) == len(S) + len(T):
                yield S, T, w
