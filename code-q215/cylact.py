"""The nilCoxeter action u_i on cylindric shapes, in the bead model of
2026-09-20-c1-cylindric-M-convexity.tex, and the bridge to affine Stanley
symmetric polynomials.

CONVENTION, DERIVED NOT TUNED.
  A cylindric shape of type (n,m) is x : Z -> Z strictly increasing with
  x_{i+m} = x_i + n.  Its bead set is A = {x_i} subset Z, with A + n = A and
  |A cap [0,n)| = m.  Let q_a = [a in A].
  Adding one box (per period) to row j means x_j -> x_j + 1, legal iff
  x_j + 1 not in A, i.e. (q_a, q_{a+1}) = (1,0) -> (0,1) at a = x_j.
  Lam math/0501335 sec 7 states his rule on the edge sequence p as
  (p_i, p_{i+1}) = (0,1) -> (1,0).  Matching the two forces  p = 1 - q
  (p is the HOLE indicator) and forces Lam's diagonal index to be my bead
  position, read in the same (increasing) direction.  Nothing is free except
  the ORIGIN of the labelling, and F~_w is invariant under rotation of the
  labels (Lam, Prop at l.1113: F~_w = F~_{p.w}).
"""
from itertools import combinations
from functools import lru_cache
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q209')
from affstan import identity, rmul_s, length, cyc_dec_elements

# ------------------------------------------------------------------ shapes
# A shape is the tuple (x_1,...,x_m) with x_1 < ... < x_m < x_1 + n.

def ext(x, i, n, m):
    """x_i for arbitrary i in Z."""
    q, r = divmod(i - 1, m)
    return x[r] + q * n

def is_shape(x, n, m):
    return len(x) == m and all(x[i] < x[i+1] for i in range(m-1)) and x[m-1] < x[0] + n

def normalise(x, n, m):
    """Rotate/translate the tuple so that x_1 <= ... <= x_m < x_1 + n holds with
    the same bead SET and the same indexing shift.  We keep the indexing, so all
    we do is re-sort into the canonical window representative starting at x[0]."""
    return tuple(x)

def beadset_res(x, n, m):
    return frozenset(v % n for v in x)

def u_act(x, i, n, m):
    """u_i . x .  Returns None for 0.  Moves every bead at residue i to i+1."""
    res = beadset_res(x, n, m)
    i %= n
    if i not in res or (i + 1) % n in res:
        return None
    return tuple(sorted(v + 1 if v % n == i else v for v in x))

def u_word(x, word, n, m):
    """Apply u_{word[0]} u_{word[1]} ... (leftmost acts LAST, operator order)."""
    cur = x
    for i in reversed(word):
        cur = u_act(cur, i, n, m)
        if cur is None:
            return None
    return cur

# ------------------------------------------------------- horizontal strips
def hstrip(nu, rho, n, m):
    """nu <<  rho : nu_i <= rho_i < nu_{i+1} for all i."""
    if any(rho[j] < nu[j] for j in range(m)):
        return False
    for j in range(m):
        if not (rho[j] < ext(nu, j + 2, n, m)):
            return False
    return True

def usize(x):
    return sum(x)

# ------------------------------------------------------ 321-avoidance
def reduced_words(w, n, maxlen=12):
    """All reduced words for w (as tuples of generators, LEFT-to-RIGHT so that
    w = s_{a_1} ... s_{a_L} acting on positions)."""
    L = length(w, n)
    if L > maxlen:
        raise ValueError("too long")
    out = []
    def rec(cur, acc):
        if length(cur, n) == 0:
            out.append(tuple(reversed(acc)))
            return
        for i in range(n):
            v = rmul_s(cur, i, n)
            if length(v, n) == length(cur, n) - 1:
                rec(v, acc + [i])
    rec(w, [])
    return out

def is_321_avoiding_words(w, n):
    """Full commutativity: no reduced word contains a braid FACTOR
    s_i s_{i+-1} s_i (three CONSECUTIVE letters).

    NB.  Lam math/0501335 sec 8 writes "subsequence"; taken literally that is a
    different and strictly stronger condition -- e.g. w = (-3,4,5) in S~_3 has
    the single reduced word 0 1 2 0, which contains 0,1,0 as a subsequence but
    no braid factor, and it has no 321 pattern in the sense of his own
    Proposition.  The proof of that Proposition ("so that w = v s_i s_{i+1} s_i u")
    shows the factor is what is meant.  This implementation uses the factor
    reading; is_321_avoiding_pattern is an independent check of it."""
    for word in reduced_words(w, n):
        for a in range(len(word) - 2):
            i, j, k = word[a], word[a+1], word[a+2]
            if i == k and (j == (i + 1) % n or j == (i - 1) % n) and n > 2:
                return False
    return True

def wval(w, i, n):
    """w(i) for i in Z, from the window (w(1),...,w(n))."""
    q, r = divmod(i - 1, n)
    return w[r] + q * n

def is_321_avoiding_pattern(w, n, K=3):
    """Lam's Prop: 321-avoiding iff no x<y<z in Z with w(x)>w(y)>w(z)."""
    rng = range(1 - K * n, (K + 1) * n + 1)
    vals = [(i, wval(w, i, n)) for i in rng]
    for a in range(len(vals)):
        for b in range(a+1, len(vals)):
            if vals[b][1] >= vals[a][1]:
                continue
            for c in range(b+1, len(vals)):
                if vals[c][1] < vals[b][1]:
                    return False
    return True

# ------------------------------------------------------ Lam's mu(w)
def code_inverse(w, n):
    """c'(w) with c'_{w(i)} = #{j < i : w(j) > w(i)}.  Returns (c'_1,...,c'_n)."""
    c = {}
    for i in range(1, n + 1):
        wi = wval(w, i, n)
        cnt = 0
        j = i - 1
        # w(j) > w(i) with j < i ; only finitely many since w(j) -> -inf as j -> -inf
        while True:
            wj = wval(w, j, n)
            if wj > wi:
                cnt += 1
            # stop when all further j have w(j) < w(i) for sure:
            # w(j) <= max_window + floor((j-1)/n)*n ; safe bound:
            if wi - wj > n * n + 2 * n:
                break
            j -= 1
        c[wi % n] = cnt
    return tuple(c[a % n] for a in range(1, n + 1))

def conjugate(part):
    part = [p for p in part if p > 0]
    if not part:
        return ()
    M = max(part)
    return tuple(sum(1 for p in part if p >= k) for k in range(1, M + 1))

def mu_of_w(w, n):
    """Lam thm:monomial: mu(w) = conjugate of the decreasing rearrangement of c'(w)."""
    return conjugate(sorted(code_inverse(w, n), reverse=True))

# ------------------------------------------------------ greedy lambda-hat
def greedy_chain(mu, lam, n, m, maxsteps=60):
    """g^0 = mu, g^t_i = min(g^{t-1}_{i+1} - 1, lam_i).  Returns (gammas, ellmin)."""
    g = tuple(mu)
    gammas = []
    for t in range(1, maxsteps + 1):
        gn = tuple(min(ext(g, j + 2, n, m) - 1, lam[j]) for j in range(m))
        gammas.append(usize(gn) - usize(g))
        g = gn
        if g == tuple(lam):
            return gammas, t
    return None, None

def lambda_hat(mu, lam, n, m):
    gammas, ellmin = greedy_chain(mu, lam, n, m)
    if gammas is None:
        return None, None
    return tuple(sorted((x for x in gammas if x > 0), reverse=True)), ellmin

# ------------------------------------------------ s^c in ell variables
def cyl_poly_support(mu, lam, n, m, ell):
    """supp of s^c_{lam/mu}(x_1..x_ell): set of weight vectors of chains, with
    multiplicities.  Returns dict alpha -> count."""
    layer = {tuple(mu): {(): 1}}
    for t in range(ell):
        new = {}
        for cur, paths in layer.items():
            for nxt in shapes_between(cur, lam, n, m):
                if not hstrip(cur, nxt, n, m):
                    continue
                k = usize(nxt) - usize(cur)
                tgt = new.setdefault(nxt, {})
                for a, cnt in paths.items():
                    key = a + (k,)
                    tgt[key] = tgt.get(key, 0) + cnt
        layer = new
    return layer.get(tuple(lam), {})

def shapes_between(nu, lam, n, m):
    """All cylindric shapes rho with nu subset rho subset lam (componentwise)."""
    ranges = [range(nu[j], lam[j] + 1) for j in range(m)]
    out = []
    def rec(j, acc):
        if j == m:
            t = tuple(acc)
            if is_shape(t, n, m):
                out.append(t)
            return
        for v in ranges[j]:
            if j and v <= acc[-1]:
                continue
            rec(j + 1, acc + [v])
    rec(0, [])
    return out
