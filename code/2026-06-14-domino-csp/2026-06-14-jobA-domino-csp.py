#!/usr/bin/env python3
"""
JOB A -- does the descent-graded order law live on the Colmenarejo-Tenner-Thompson
(CTT, arXiv 2602.23343) domino-tableau CSP?

The order law (PROVED): ord_{q=-1} G_lambda(q) = floor(|2-core(lambda)|/2),
  G_lambda(q) = sum_{T in SYT(lambda)} q^{s(T)},   s(T) = sum_{i in Des(T)} w_i,
  w_i = 2i-1 if (n-i) odd else 0.

CTT CSP triple (from my reading log 2026-06-23):
  (DT(square), C = <g> of order n, X_n(q) = [n choose floor(n/2)]_q = sum_D q^{maj(D)}),
  g = SUBSET CYCLIC SHIFT on the binary-word image (NOT Schuetzenberger promotion).
  |DT| = C(n, floor(n/2)).

This script:
  (A1) Gaussian q-binomial + verify the CTT CSP on binary words W_{n,k} under cyclic shift.
  (A2) Enumerate standard domino tableaux (SDT) of candidate shapes; find the shape +
       maj/descent convention that reproduces [n choose floor(n/2)]_q.  (Grounds the object.)
  (A3) Compare Sum_{D in SDT(lambda)} q^{maj(D)} against my G_lambda(q):
       as polynomials, at q=-1 (order), and the 2-core-SIZE vs content guard.
  (A4) OC-trap guard on staircases delta_k -- but staircases are NOT rectangular,
       so report whether the CTT domain even contains the discriminators.

Pure Python (no SageMath in this container).
"""
from sympy import symbols, Poly, factor, expand, simplify, I, exp, pi, Rational, nsimplify
import sympy
from functools import lru_cache
from collections import Counter, defaultdict
from math import factorial

q = symbols('q')

# ===========================================================================
# (A1)  Gaussian q-binomial  and  CTT CSP on binary words
# ===========================================================================
def q_int(n):
    return sum(q**i for i in range(n))

def q_fact(n):
    r = Poly(1, q)
    for i in range(1, n+1):
        r = r * Poly(q_int(i), q)
    return r

def q_binom(n, k):
    """Gaussian binomial [n choose k]_q as a sympy Poly in q."""
    if k < 0 or k > n:
        return Poly(0, q)
    num = q_fact(n)
    den = q_fact(k) * q_fact(n-k)
    quo, rem = sympy.div(num.as_expr(), den.as_expr(), q)
    assert sympy.simplify(rem) == 0
    return Poly(quo, q)

def maj_word(w):
    """major index of a binary word: sum of positions i (1-indexed) with w[i-1]>w[i]."""
    return sum(i+1 for i in range(len(w)-1) if w[i] > w[i+1])

def words_binary(n, k):
    """all binary words length n with k ones."""
    from itertools import combinations
    res = []
    for ones in combinations(range(n), k):
        w = [0]*n
        for o in ones:
            w[o] = 1
        res.append(tuple(w))
    return res

def csp_check_binary(n):
    """Verify (W_{n,floor(n/2)}, cyclic shift order n, X_n(q)) is a CSP triple."""
    k = n//2
    W = words_binary(n, k)
    X = q_binom(n, k)
    # maj-generating function over W should equal X
    majgen = Poly(sum(q**maj_word(w) for w in W), q)
    polys_match = (majgen - X).is_zero
    # cyclic shift action: w -> shift left by 1 (cyclic).  fixed points of g^d
    def shift(w, d):
        return tuple(w[(i+d) % n] for i in range(n))
    zeta = exp(2*sympy.pi*I/n)
    csp_ok = True
    rows = []
    for d in range(n):
        fixed = sum(1 for w in W if shift(w, d) == w)
        Xval = X.as_expr().subs(q, zeta**d)
        Xval = sympy.nsimplify(sympy.expand(sympy.simplify(Xval)))
        Xval_round = complex(Xval)
        ok = abs(Xval_round.real - fixed) < 1e-6 and abs(Xval_round.imag) < 1e-6
        csp_ok = csp_ok and ok
        rows.append((d, fixed, Xval, ok))
    return polys_match, csp_ok, rows, X, len(W)

# ===========================================================================
# (A2)  Standard domino tableaux of a shape, and candidate maj statistics
# ===========================================================================
def domino_removals(lam):
    """All ways to remove one size-2 border strip (domino) from partition lam.
    Yields (mu, orient, rowtop) where orient in {'H','V'}, rowtop = top row (0-idx)."""
    lam = [p for p in lam if p > 0]
    L = len(lam)
    out = []
    # horizontal: remove last two cells of row r, both in row r
    for r in range(L):
        c = lam[r]
        if c >= 2:
            below = lam[r+1] if r+1 < L else 0
            if below <= c-2:           # removal keeps a Young diagram
                mu = lam[:]
                mu[r] = c-2
                out.append((tuple(p for p in mu if p > 0), 'H', r))
    # vertical: remove (r,c),(r+1,c) in column c=lam[r]=lam[r+1]
    for r in range(L-1):
        if lam[r] == lam[r+1]:
            c = lam[r]
            below2 = lam[r+2] if r+2 < L else 0
            if below2 <= c-1:
                mu = lam[:]
                mu[r] = c-1
                mu[r+1] = c-1
                out.append((tuple(p for p in mu if p > 0), 'V', r))
    return out

@lru_cache(maxsize=None)
def is_domino_tileable(lam):
    """A partition is tileable by dominoes iff its 2-core is empty."""
    lam = list(lam)
    while True:
        rem = domino_removals(tuple(lam))
        if not rem:
            break
        lam = list(rem[0][0])
    return len(lam) == 0

def standard_domino_tableaux(lam):
    """Enumerate SDT of shape lam as tuples of placements.
    A SDT = chain () = mu_0 < mu_1 < ... < mu_n = lam adding one domino each step.
    Return list of placements; placement k = (orient, rowtop) of domino labeled k+1.
    We build by REMOVING dominoes from lam (label = current size//2, decreasing)."""
    n = sum(lam)//2
    results = []
    def rec(cur, placements):
        if sum(cur) == 0:
            results.append(list(reversed(placements)))
            return
        seen = set()
        for mu, orient, rowtop in domino_removals(cur):
            key = (mu, orient, rowtop)
            if key in seen:
                continue
            seen.add(key)
            rec(mu, placements + [(orient, rowtop)])
    rec(tuple(lam), [])
    return results

# --- candidate descent statistics on a SDT (list of (orient,rowtop), label 1..n) ---
def maj_candidates(placements):
    """Return dict of candidate maj values under several descent conventions.
    placements[i] = (orient, rowtop) of domino labeled i+1."""
    n = len(placements)
    def rowtop(i):  # top row of domino i (1-indexed label)
        return placements[i-1][1]
    def rowbot(i):  # bottom row (rowtop, or rowtop+1 if vertical)
        o, rt = placements[i-1]
        return rt + (1 if o == 'V' else 0)
    cands = {}
    # convention 1: descent if top row of i+1 > top row of i ; maj = sum i
    cands['top_gt_sum_i'] = sum(i for i in range(1, n) if rowtop(i+1) > rowtop(i))
    # convention 2: descent if bottom row of i+1 > bottom row of i
    cands['bot_gt_sum_i'] = sum(i for i in range(1, n) if rowbot(i+1) > rowbot(i))
    # convention 3: descent if top row i+1 >= top row i AND not strictly left... use >=
    cands['top_ge_sum_i'] = sum(i for i in range(1, n) if rowtop(i+1) >= rowtop(i))
    # convention 4: descent if i+1 strictly below i (top), maj = sum of #descents (count)
    cands['top_gt_count'] = sum(1 for i in range(1, n) if rowtop(i+1) > rowtop(i))
    return cands

def sdt_majgen(lam, conv):
    sdts = standard_domino_tableaux(lam)
    return Poly(sum(q**maj_candidates(p)[conv] for p in sdts), q), len(sdts)

# ===========================================================================
# G_lambda(q) = sum_{T in SYT} q^{s(T)}  (my order-law polynomial)
# ===========================================================================
def syt_list(lam):
    cells = [(i, j) for i, p in enumerate(lam) for j in range(p)]
    n = len(cells)
    results = []
    target = list(lam)
    rowofval = [0]*(n+1)
    shape = [0]*(len(lam)+1)
    def rec(v):
        if v > n:
            results.append(tuple(rowofval[1:])); return
        for r in range(len(target)):
            if shape[r] < target[r] and (r == 0 or shape[r] < shape[r-1]):
                shape[r] += 1; rowofval[v] = r; rec(v+1); shape[r] -= 1
    rec(1)
    return results, n

def s_statistic(rowofval, n):
    s = 0
    for i in range(1, n):
        if rowofval[i] > rowofval[i-1]:
            s += (2*i-1) if ((n-i) % 2 == 1) else 0
    return s

def G_poly(lam):
    sl, n = syt_list(lam)
    return Poly(sum(q**s_statistic(rv, n) for rv in sl), q)

def root_mult(poly, root):
    if poly.is_zero:
        return None
    mult = 0; p = poly
    while p.eval(root) == 0 and not p.is_zero:
        quo, rem = sympy.div(p.as_expr(), (q-root), q)
        p = Poly(quo, q); mult += 1
        if mult > 500: break
    return mult

def two_core(lam):
    lam = list(lam)
    while True:
        rem = domino_removals(tuple(lam))
        if not rem: break
        lam = list(rem[0][0])
    return tuple(lam)

# ===========================================================================
# main
# ===========================================================================
def hdr(t):
    print("="*78); print(t); print("="*78)

def run():
    # ----- A1: ground the CTT object -----
    hdr("A1.  CTT CSP on binary words: X_n(q)=[n choose floor(n/2)]_q, cyclic shift")
    print(f"{'n':>3} {'|W|=C(n,k)':>11} {'polys match':>12} {'CSP holds':>10}   X_n(q)")
    for n in range(2, 9):
        pm, csp, rows, X, sz = csp_check_binary(n)
        print(f"{n:>3} {sz:>11} {str(pm):>12} {str(csp):>10}   {X.as_expr()}")
    print("\n  -> X_n(q) IS the maj generating function over binary words, and the")
    print("     cyclic-shift CSP holds.  CTT object grounded.")

    # ----- A2: which SHAPE + maj convention reproduces [n choose floor(n/2)]_q ? -----
    hdr("A2.  Standard domino tableaux: find shape+maj reproducing [n choose k]_q")
    convs = ['top_gt_sum_i', 'bot_gt_sum_i', 'top_ge_sum_i', 'top_gt_count']
    # candidate shapes: 2xN rectangle (n,n), and the n x n square (n^n) when even area
    print("Shape family A: two-row rectangle (N,N)  [the '2 x N' reading]")
    print(f"{'N':>3} {'|SDT|':>6} {'C(N,fl)':>8}  match-convention(s)")
    for N in range(1, 8):
        lam = (N, N)
        if not is_domino_tileable(lam):
            print(f"{N:>3}  not tileable"); continue
        sdts = standard_domino_tableaux(lam)
        target = q_binom(N, N//2)
        matches = []
        for cv in convs:
            mg = Poly(sum(q**maj_candidates(p)[cv] for p in sdts), q)
            if (mg - target).is_zero:
                matches.append(cv)
        from math import comb
        print(f"{N:>3} {len(sdts):>6} {comb(N, N//2):>8}  {matches}")

    print("\nShape family B: square (N^N) (N rows of length N), even area only")
    print(f"{'N':>3} {'area':>5} {'|SDT|':>7} {'C(N,fl)':>8}  match-convention(s)")
    from math import comb
    for N in range(2, 6):
        lam = tuple([N]*N)
        if (N*N) % 2 != 0 or not is_domino_tileable(lam):
            print(f"{N:>3} {N*N:>5}  odd/untileable"); continue
        sdts = standard_domino_tableaux(lam)
        target = q_binom(N, N//2)
        matches = []
        for cv in convs:
            mg = Poly(sum(q**maj_candidates(p)[cv] for p in sdts), q)
            if (mg - target).is_zero:
                matches.append(cv)
        print(f"{N:>3} {N*N:>5} {len(sdts):>7} {comb(N, N//2):>8}  {matches}")

    # ----- A3: compare SDT-maj polynomial to G_lambda(q) on shared shapes -----
    hdr("A3.  Compare  Sum_{SDT(lam)} q^{maj} vs G_lam(q)=Sum_{SYT} q^{s(T)}")
    print("For each tileable lam: |SYT| vs |SDT| (must match for a poly identity),")
    print("ord_{q=-1} of each, and floor(|2-core|/2).")
    print(f"{'lam':>10} {'|SYT|':>6} {'|SDT|':>6} {'ordG@-1':>8} {'2core':>8} {'fl/2':>5}")
    from itertools import count
    def partitions(n, mx=None):
        mx = n if mx is None else mx
        if n == 0:
            yield (); return
        for f in range(min(n, mx), 0, -1):
            for rest in partitions(n-f, f):
                yield (f,)+rest
    for n in range(2, 11, 2):
        for lam in partitions(n):
            if not is_domino_tileable(lam):
                continue
            sl, _ = syt_list(lam)
            nsyt = len(sl)
            sdts = standard_domino_tableaux(lam)
            G = G_poly(lam)
            ordG = root_mult(G, -1)
            core = two_core(lam)   # empty for tileable
            print(f"{str(lam):>10} {nsyt:>6} {len(sdts):>6} {str(ordG):>8} {str(core):>8} {sum(core)//2:>5}")

    print("\n  NOTE: tileable lam <=> empty 2-core <=> floor(|2-core|/2)=0 <=> ord_{q=-1}G=0.")
    print("  So on the domino-tileable shapes (CTT's domain) the order law is TRIVIAL (0).")

    # ----- A4: OC-trap guard -- are the discriminating staircases in CTT's domain? -----
    hdr("A4.  OC-trap guard: are the discriminating staircases delta_k tileable (in CTT domain)?")
    print(f"{'k':>3} {'delta_k':>16} {'tileable?':>10} {'2-core':>16} {'ord=fl(|core|/2)':>16}")
    for k in range(1, 8):
        dk = tuple(range(k, 0, -1))
        til = is_domino_tileable(dk)
        core = two_core(dk)
        print(f"{k:>3} {str(dk):>16} {str(til):>10} {str(core):>16} {sum(core)//2:>16}")

    # ----- A5: CTT/SDT poly reads the 2-QUOTIENT; order law reads the 2-CORE -----
    hdr("A5.  Complementarity: |SDT(lam)| = multinomial x f^{q0} f^{q1} (2-quotient)")
    print("If |SDT| = C(m,|q0|/?) * f^{q0} * f^{q1}, the domino-CSP reads only the 2-quotient,")
    print("which is the half of lam COMPLEMENTARY to the 2-core the order law reads.")
    print(f"{'lam':>10} {'|SDT|':>6} {'2-quotient (q0,q1)':>22} {'f^q0*f^q1*C(m,|q0|)':>20} {'ok':>4}")
    for n in range(2, 11, 2):
        for lam in partitions(n):
            if not is_domino_tileable(lam):
                continue
            sdts = standard_domino_tableaux(lam)
            q0, q1 = two_quotient(lam)
            m = n//2
            pred = comb(m, sum(q0)) * hook_f_local(q0) * hook_f_local(q1)
            ok = (pred == len(sdts))
            if not ok or n <= 6:
                print(f"{str(lam):>10} {len(sdts):>6} {str((q0,q1)):>22} {pred:>20} {str(ok):>4}")
    print("\n  -> |SDT(lam)| factors through the 2-QUOTIENT only.  The 2-core is empty on")
    print("     CTT's whole domain, so the CTT q-analogue is BLIND to the 2-core --")
    print("     which is precisely the invariant the order law reads.  Orthogonal halves.")

    # direct polynomial comparison where |SYT|=|SDT|
    hdr("A5b. Direct G_lam(q) vs CTT-SDT-maj(q) on tileable lam")
    print(f"{'lam':>10} {'G_lam(q)':>30} {'SDT-maj(q)':>26} {'equal?':>7}")
    for n in range(2, 9, 2):
        for lam in partitions(n):
            if not is_domino_tileable(lam):
                continue
            G = G_poly(lam)
            sdts = standard_domino_tableaux(lam)
            S = Poly(sum(q**maj_candidates(p)['top_gt_sum_i'] for p in sdts), q)
            eq = (G - S).is_zero
            if eq:
                print(f"{str(lam):>10} {str(G.as_expr()):>30} {str(S.as_expr()):>26} {str(eq):>7}")
    print("  (only the trivially-1-term / coincidental shapes coincide; see counts in A3)")

def two_quotient(lam):
    """2-quotient of partition lam via beta-numbers split by parity.
    Returns (q0, q1) partitions."""
    lam = [p for p in lam if p > 0]
    r = len(lam)
    # pad to even length
    if r % 2 == 1:
        lam = lam + [0]; r += 1
    beta = sorted([lam[i] + (r-1-i) for i in range(r)], reverse=True)
    even = sorted([b for b in beta if b % 2 == 0], reverse=True)
    odd = sorted([b for b in beta if b % 2 == 1], reverse=True)
    def to_part(bs):
        bs = sorted([b//2 for b in bs], reverse=True)
        k = len(bs)
        p = [bs[i] - (k-1-i) for i in range(k)]
        return tuple(x for x in p if x > 0)
    return to_part(even), to_part(odd)

def hook_f_local(lam):
    lam = [x for x in lam if x > 0]
    if not lam:
        return 1
    n = sum(lam); m0 = lam[0]
    conj = [sum(1 for rr in lam if rr > c) for c in range(m0)]
    prod = 1
    for i, row in enumerate(lam):
        for jc in range(row):
            prod *= (row-jc-1) + (conj[jc]-i-1) + 1
    return factorial(n)//prod

if __name__ == "__main__":
    run()
