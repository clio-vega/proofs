"""Q213: compare Clio's Theorem 4.1 exchange move with Morse-Schilling crystal operators.
Affine symmetric group S~_n via window notation [w(1),...,w(n)]."""
from itertools import combinations, chain

def identity(n): return list(range(1,n+1))

def right_s(w, i, n):
    """w * s_i.  i in 0..n-1.  s_i swaps positions i,i+1 (mod n)."""
    w = list(w)
    if i == 0:
        a, b = w[0], w[n-1]
        return [b-n] + w[1:n-1] + [a+n]
    else:
        w[i-1], w[i] = w[i], w[i-1]
        return w

def length(w, n):
    tot = 0
    for i in range(n):
        for j in range(i+1, n):
            tot += abs((w[j]-w[i])//n)
    return tot

def mult(w1, w2, n):
    """w1*w2 as functions: (w1 w2)(i) = w1(w2(i)).  window of w1w2."""
    def ev(w, i):
        # w(i) for any integer i
        q, r = divmod(i-1, n)
        return w[r] + q*n
    return [ev(w1, ev(w2, i)) for i in range(1, n+1)]

def rank(i, x, n):
    """position in decreasing order x-1 > x-2 > ... > 0 > n-1 > ... > x+1"""
    return (x-1-i) % n

def cyc_dec(S, n, x=None):
    """cyclically decreasing element u_S, S subset of Z/nZ, S != Z/nZ."""
    S = set(s % n for s in S)
    assert len(S) < n
    if x is None:
        x = min(set(range(n)) - S)
    order = sorted(S, key=lambda i: rank(i, x, n))   # decreasing wrt x
    w = identity(n)
    for i in order:
        w = right_s(w, i, n)
    return w

def check_cycdec_welldefined(n):
    for k in range(0, n):
        for S in combinations(range(n), k):
            avail = [x for x in range(n) if x not in S]
            ws = set(tuple(cyc_dec(S, n, x)) for x in avail)
            if len(ws) != 1: return (S, ws)
            if length(list(ws.pop()), n) != len(S): return ('len', S)
    return None

# ---------- MS pairing ----------
def ms_pairing(Su, Sv, x, n):
    """Su=content(u) (LEFT factor), Sv=content(v) (RIGHT).  Returns (L1,R1,pairs)."""
    uo = sorted(Su, key=lambda i: rank(i, x, n))     # decreasing order
    unpaired_v = set(Sv)
    pairs = []
    L1 = []
    for b in uo:
        # smallest a > b in content(v), unpaired: a with rank(a)<rank(b), maximal rank
        cands = [a for a in unpaired_v if rank(a, x, n) < rank(b, x, n)]
        if cands:
            a = max(cands, key=lambda i: rank(i, x, n))
            unpaired_v.discard(a); pairs.append((b, a))
        else:
            L1.append(b)
    return set(L1), set(unpaired_v), pairs

def ms_e1(Su, Sv, x, n):
    L1, R1, _ = ms_pairing(Su, Sv, x, n)
    if not L1: return None
    b = max(L1, key=lambda i: rank(i, x, n))          # min wrt order = largest rank
    t = 0
    while (b-t-1) % n in Su: t += 1
    return (Su - {b}, Sv | {(b-t) % n}, b, t, (b-t) % n)

def ms_f1(Su, Sv, x, n):
    L1, R1, _ = ms_pairing(Su, Sv, x, n)
    if not R1: return None
    a = min(R1, key=lambda i: rank(i, x, n))          # max wrt order = smallest rank
    s = 0
    while (a+s+1) % n in Sv: s += 1
    return (Su | {(a+s) % n}, Sv - {a}, a, s, (a+s) % n)

# ---------- Clio's Theorem 4.1 ----------
def runs(S, n):
    """maximal cyclic intervals [m,M] of S, as (m, list of elements in order m,m+1,...,M)"""
    S = set(S)
    if not S or len(S) == n: return []
    out = []
    for m in S:
        if (m-1) % n in S: continue
        elts = [m]
        j = (m+1) % n
        while j in S:
            elts.append(j); j = (j+1) % n
        out.append((m, elts))
    return out

def clio_moves(S, T, n):
    """all (run m, e, S', T') admissible under Thm 4.1"""
    out = []
    for m, elts in runs(S, n):
        if m in T: continue
        for e in elts:                      # order m, m+1, ..., M
            if (e+1) % n not in T:
                out.append((m, e, S - {e}, T | {m}))
                break                       # FIRST such e
    return out

def all_pairs(n):
    subsets = [frozenset(c) for k in range(n) for c in combinations(range(n), k)]
    for S in subsets:
        for T in subsets:
            uS = cyc_dec(S, n); uT = cyc_dec(T, n)
            w = mult(uS, uT, n)
            if length(w, n) == len(S)+len(T):
                yield set(S), set(T), tuple(w)
