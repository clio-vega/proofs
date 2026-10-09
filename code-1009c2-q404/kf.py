"""PARTIALLY ABANDONED -- read this before using anything here.

  * KF() / charge() / ssyt(): the Lascoux-Schuetzenberger charge engine.  DISCARDED.
    It was wrong twice (first cocharge, then a surviving convention error: the asymmetric
    test deg K_{lam,mu} = n(mu) - n(lam) with leading coefficient 1 failed on 133 of 471
    nonzero pairs, e.g. K_{(4,1),(2,2,1)} read 2*t**3).  All three SYMMETRIC tests passed
    throughout.  NOTHING in the 2026-10-09-c2 paper depends on it.  Ground truth for
    Hall-Littlewood lives in hl.py, built from the definition.

  * chi() / border_strips(): Murnaghan-Nakayama characters.  THESE ARE USED, as the
    independent mechanism for Check 1 (Y^mu_rho(0) = chi^mu_rho, 434 pairs n<=7,
    0 failures).  They are not implicated in the charge fault.
"""

import sympy as sp
from functools import lru_cache
from itertools import product as iproduct

t = sp.Symbol('t')

# ---------- partitions ----------
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

# ---------- SSYT of shape lam, content mu ----------
def ssyt(lam, mu):
    """All semistandard tableaux of shape lam with content mu (list of rows)."""
    cells = [(i, j) for i in range(len(lam)) for j in range(lam[i])]
    if sum(lam) != sum(mu): return
    k = len(mu)
    rows = [[None]*lam[i] for i in range(len(lam))]
    remaining = list(mu)
    def rec(idx):
        if idx == len(cells):
            yield [tuple(r) for r in rows]
            return
        i, j = cells[idx]
        lo = 1
        if j > 0: lo = max(lo, rows[i][j-1])          # weakly increasing along rows
        if i > 0: lo = max(lo, rows[i-1][j] + 1)      # strictly increasing down columns
        for v in range(lo, k+1):
            if remaining[v-1] == 0: continue
            rows[i][j] = v; remaining[v-1] -= 1
            yield from rec(idx+1)
            remaining[v-1] += 1; rows[i][j] = None
    yield from rec(0)

def reading_word(rows):
    """Macdonald's word: read rows from the TOP row down, each row RIGHT to LEFT.
       (tested below against known K values)"""
    w = []
    for r in rows:
        w.extend(reversed(r))
    return w

# ---------- charge ----------
def _standard_subword_positions(w):
    """Extract one standard subword from word w (list of values, positions implicit).
       Returns list of positions, in increasing letter order 1,2,3,...
       Rule: start at the RIGHTMOST 1.  Having chosen position p for letter i,
       choose for i+1 the first occurrence of i+1 strictly to the LEFT of p,
       scanning leftward and wrapping cyclically to the right end."""
    n = len(w)
    maxv = max(w)
    pos = {}
    # rightmost 1
    p = max(i for i in range(n) if w[i] == 1)
    pos[1] = p
    for v in range(2, maxv+1):
        # scan leftward from p-1, wrapping
        found = None
        for step in range(1, n+1):
            q = (p - step) % n
            if w[q] == v:
                found = q; break
        if found is None: return None
        pos[v] = found; p = found
    return [pos[v] for v in range(1, maxv+1)]

def charge_standard(positions):
    """positions[v-1] = position of letter v in the word.  index(1)=0;
       index(v+1) = index(v) + 1 if v+1 is to the RIGHT of v, else index(v)."""
    idx = [0]*len(positions)
    for v in range(1, len(positions)):
        # Macdonald III (6): index(v+1) = index(v) + 1 if v+1 lies to the LEFT of v,
        # else index(v).  (The opposite convention computes a different statistic --
        # measured 2026-10-09: it put leading coefficient 2 on K_{(4,1),(2,2,1)}.)
        idx[v] = idx[v-1] + (1 if positions[v] < positions[v-1] else 0)
    return sum(idx)

def charge(w):
    """Charge of a word of partition content."""
    w = list(w)
    total = 0
    while w:
        pos = _standard_subword_positions(w)
        if pos is None:
            raise ValueError("content not a partition: %r" % (w,))
        total += charge_standard(sorted_by_letter(pos))
        keep = [i for i in range(len(w)) if i not in set(pos)]
        w = [w[i] for i in keep]
    return total

def sorted_by_letter(pos):
    # pos already given in letter order 1,2,...  but positions are raw indices
    return pos

@lru_cache(maxsize=None)
def kostka_foulkes(lam, mu):
    """K_{lam,mu}(t) = sum over SSYT(lam,mu) of t^charge(word)."""
    out = sp.Integer(0)
    for rows in ssyt(lam, mu):
        out += t**charge(reading_word(rows))
    return sp.expand(out)

# ---------- the statistic above is COCHARGE (measured: K_{lam,lam} read t^{n(lam)}
# ---------- and K_{(n),mu} read 1 -- exactly swapped).  charge = n(content) - cocharge.
def n_stat(mu):
    return sum(i*m for i, m in enumerate(mu))

@lru_cache(maxsize=None)
def KF(lam, mu):
    """Kostka-Foulkes K_{lam,mu}(t) = sum_{T in SSYT(lam,mu)} t^{charge(T)}."""
    out = sp.Integer(0)
    for rows in ssyt(lam, mu):
        out += t**charge(reading_word(rows))
    return sp.expand(out)

# ---------- Murnaghan-Nakayama characters ----------
def border_strips(lam, r):
    """All (lam_minus, height) for removing a border strip of size r from lam.
       Uses the beta-number / first-column hook-length encoding."""
    res = []
    l = len(lam)
    # beta numbers: lam_i + (l - i) - 1 for i=0..l-1  (distinct, decreasing)
    beta = [lam[i] + (l - 1 - i) for i in range(l)]
    bset = set(beta)
    for b in beta:
        nb = b - r
        if nb >= 0 and nb not in bset:
            nbeta = sorted([x for x in beta if x != b] + [nb], reverse=True)
            # height = number of rows the strip occupies - 1
            #        = #{beta elements strictly between nb and b}
            ht = sum(1 for x in beta if nb < x < b)
            newlam = tuple(nbeta[i] - (l - 1 - i) for i in range(l))
            newlam = tuple(x for x in newlam if x > 0)
            if any(x < 0 for x in (nbeta[i] - (l - 1 - i) for i in range(l))):
                continue
            res.append((newlam, ht))
    return res

@lru_cache(maxsize=None)
def chi(lam, rho):
    """Irreducible S_n character chi^lam evaluated on cycle type rho."""
    if sum(lam) == 0: return 1
    if not rho: return 0
    r = rho[0]
    rest = rho[1:]
    tot = 0
    for nl, ht in border_strips(lam, r):
        tot += (-1)**ht * chi(nl, rest)
    return tot

# ---------- Green polynomials ----------
@lru_cache(maxsize=None)
def Y(mu, rho):
    """Y^mu_rho(t) = [P_mu] p_rho = sum_lam chi^lam_rho K_{lam,mu}(t)."""
    n = sum(mu)
    out = sp.Integer(0)
    for lam in partitions(n):
        c = chi(lam, rho)
        if c: out += c * KF(lam, mu)
    return sp.expand(out)
