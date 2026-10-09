"""Kostka-Foulkes engine, ENGINE A: Lascoux-Schutzenberger charge.

Conventions (Macdonald III (6.5) / Lascoux-Schutzenberger):
  K_{lam,mu}(t) = sum_{T in SSYT(lam,mu)} t^{charge(T)}
  charge of a word of partition content: decompose into standard subwords,
  charge = sum of charges of the standard subwords.
  For a standard word w (one copy each of 1..n):  index(1)=0, and for i>=2,
  index(i) = index(i-1) + 1 if i lies to the RIGHT of i-1 in w, else index(i-1).
  charge(w) = sum_i index(i).
Reading word of a tableau: read rows LEFT to RIGHT, from the BOTTOM row to the TOP row.
"""
import sys
from itertools import product
from functools import lru_cache
import sympy
from sympy import Integer, Poly, symbols

t = symbols('t')


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def nstat(lam):
    """n(lambda) = sum (i-1) lambda_i  (Macdonald I (1.5))."""
    return sum(i * p for i, p in enumerate(lam))


def dominates(lam, mu):
    """lam >= mu in dominance order (both partitions of the same n)."""
    sl = sm = 0
    L = max(len(lam), len(mu))
    for i in range(L):
        sl += lam[i] if i < len(lam) else 0
        sm += mu[i] if i < len(mu) else 0
        if sl < sm:
            return False
    return True


# ---------------------------------------------------------------- SSYT
def ssyt(lam, mu):
    """All semistandard Young tableaux of shape lam and content mu.
    Returned as tuples of row-tuples. Built by Gelfand-Tsetlin style
    horizontal-strip growth, then ENTRYWISE VERIFIED by the caller."""
    ell = len(mu)
    # chains emptyset = s0 <= s1 <= ... <= s_ell = lam, s_i/s_{i-1} horiz strip of size mu_i
    def hstrips(inner, target_size, bound):
        """all partitions outer with inner <= outer, outer/inner a horizontal strip
        of size target_size, outer contained in bound (=lam)."""
        res = []
        k = len(bound)
        inner = tuple(inner) + (0,) * (k - len(inner))
        def rec(i, cur, remaining):
            if i == k:
                if remaining == 0:
                    res.append(tuple(cur))
                return
            lo = inner[i]
            # horizontal strip: outer_i <= inner_{i-1}  (no two in same column)
            hi = bound[i]
            if i > 0:
                hi = min(hi, inner[i - 1])
            hi = min(hi, lo + remaining)
            for v in range(lo, hi + 1):
                if i + 1 < k and v < inner[i + 1]:
                    continue
                cur.append(v)
                rec(i + 1, cur, remaining - (v - lo))
                cur.pop()
        rec(0, [], target_size)
        return res

    chains = [[(0,) * len(lam)]]
    cur = [((0,) * len(lam),)]
    for i in range(ell):
        nxt = []
        for ch in cur:
            for out in hstrips(ch[-1], mu[i], lam):
                nxt.append(ch + (out,))
        cur = nxt
    out = []
    lam_pad = tuple(lam) + (0,) * 0
    for ch in cur:
        if ch[-1] != tuple(lam_pad):
            continue
        # build the tableau: entry i+1 in cells of ch[i+1]/ch[i]
        rows = [[] for _ in range(len(lam))]
        for i in range(ell):
            a, b = ch[i], ch[i + 1]
            for r in range(len(lam)):
                for _ in range(b[r] - a[r]):
                    rows[r].append(i + 1)
        out.append(tuple(tuple(r) for r in rows))
    return out


def is_ssyt(T, lam, mu):
    """ENTRYWISE verification: shape, content, weak rows, strict columns."""
    if tuple(len(r) for r in T) != tuple(lam):
        return False
    cnt = {}
    for r in T:
        for v in r:
            cnt[v] = cnt.get(v, 0) + 1
    want = {i + 1: m for i, m in enumerate(mu)}
    if cnt != want:
        return False
    for r in T:
        for a, b in zip(r, r[1:]):
            if a > b:
                return False
    for i in range(len(T) - 1):
        for j in range(len(T[i + 1])):
            if not (T[i][j] < T[i + 1][j]):
                return False
    return True


# ---------------------------------------------------------------- charge
def reading_word(T):
    """rows left-to-right, bottom row first."""
    w = []
    for r in reversed(T):
        w.extend(r)
    return tuple(w)


def charge_standard(w):
    """w a word with one copy each of 1..n (as a tuple). Returns charge."""
    pos = {v: i for i, v in enumerate(w)}
    n = len(w)
    idx = 0
    total = 0
    for v in range(2, n + 1):
        if pos[v] > pos[v - 1]:
            idx += 1
        total += idx
    return total


def charge(w):
    """charge of a word of partition content (Lascoux-Schutzenberger).
    Extract standard subwords: repeatedly, scan RIGHT to LEFT selecting the
    rightmost 1, then the first 2 strictly to its left cyclically, etc."""
    w = list(w)
    total = 0
    while w:
        n = max(w)
        # standard subword extraction
        chosen = []  # positions
        L = len(w)
        # find rightmost 1
        start = max(i for i in range(L) if w[i] == 1)
        chosen.append(start)
        cur = start
        for v in range(2, n + 1):
            # move LEFT cyclically from cur, take first occurrence of v
            found = None
            for d in range(1, L + 1):
                i = (cur - d) % L
                if w[i] == v and i not in chosen:
                    found = i
                    break
            assert found is not None, (w, v)
            chosen.append(found)
            cur = found
        sub = [w[i] for i in sorted(chosen)]
        total += charge_standard(tuple(sub))
        w = [w[i] for i in range(L) if i not in set(chosen)]
    return total


def K(lam, mu, check=True):
    """Kostka-Foulkes polynomial as a dict {exponent: coefficient}."""
    Ts = ssyt(lam, mu)
    if check:
        for T in Ts:
            assert is_ssyt(T, lam, mu), (T, lam, mu)
        assert len(set(Ts)) == len(Ts), "duplicate tableaux"
    d = {}
    for T in Ts:
        c = charge(reading_word(T))
        d[c] = d.get(c, 0) + 1
    return d, len(Ts)
