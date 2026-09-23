"""
Affine nil-Temperley-Lieb algebra A_n acting on 01-words (Postnikov math/0205165,
ll.2647-2712), compared with height-graded ribbon addition R_e(t).

CONVENTIONS FIXED HERE (0-indexed sites 0..n-1 = Postnikov's 1..n minus one):
  a_i = E_{i,i+1}: moves a bead from site i to site i+1 (mod n).
        Nonzero iff site i occupied and site i+1 empty.
        a_{n-1} carries a factor q (Postnikov's a_n = q E_{n1}).
  Operators COMPOSE RIGHT-TO-LEFT (standard): in the word w = (w_1,...,w_m)
  stored as a python list, we apply w[0] FIRST.  (This matches Postnikov's
  printed e_2 = a_2 a_1 + ... once one reads his products left-to-right;
  see notes.)
State = frozenset of occupied sites.
"""
from itertools import combinations, permutations
from collections import defaultdict
import sympy as sp

q, t = sp.symbols('q t')


def states(n, k):
    return [frozenset(c) for c in combinations(range(n), k)]


def apply_a(i, S, n):
    """a_i on state S -> (newstate, coeff) or None."""
    if i not in S or (i + 1) % n in S:
        return None
    T = set(S); T.remove(i); T.add((i + 1) % n)
    c = q if i == n - 1 else sp.Integer(1)
    return frozenset(T), c


def apply_word(word, S, n):
    """word given in APPLICATION order (word[0] applied first)."""
    cur, coeff = S, sp.Integer(1)
    for i in word:
        r = apply_a(i, cur, n)
        if r is None:
            return None
        cur, c = r
        coeff *= c
    return cur, coeff


def op_from_word(word, n, k):
    """Matrix (dict state->dict state->coeff) of the monomial."""
    M = defaultdict(dict)
    for S in states(n, k):
        r = apply_word(word, S, n)
        if r is not None:
            T, c = r
            M[S][T] = M[S].get(T, 0) + c
    return M


def add_op(A, B, scalar=1):
    for S, row in B.items():
        for T, c in row.items():
            A[S][T] = sp.expand(A[S].get(T, 0) + scalar * c)
    return A


def clean(M):
    out = defaultdict(dict)
    for S, row in M.items():
        for T, c in row.items():
            c = sp.expand(c)
            if c != 0:
                out[S][T] = c
    return out


def mul_op(A, B, n, k):
    """A then B?  No: (B*A) as operators = apply A first.  We define
    compose(A,B) = 'apply A first, then B'."""
    out = defaultdict(dict)
    for S in states(n, k):
        for T, c in A.get(S, {}).items():
            for U, d in B.get(T, {}).items():
                out[S][U] = sp.expand(out[S].get(U, 0) + c * d)
    return clean(out)


def eq_op(A, B, n, k):
    A, B = clean(A), clean(B)
    keys = set(A) | set(B)
    for S in keys:
        ra, rb = A.get(S, {}), B.get(S, {})
        if set(k2 for k2, v in ra.items() if v != 0) != set(k2 for k2, v in rb.items() if v != 0):
            return False, S
        for T in ra:
            if sp.simplify(ra[T] - rb.get(T, 0)) != 0:
                return False, S
    return True, None


# ---------------------------------------------------------------
# Intervals, orientations, and the orientation-sum operator R_e(t)
# ---------------------------------------------------------------

def interval(i, e, n):
    return [(i + s) % n for s in range(e)]


def word_from_orientation(J, eps):
    """J = [j_0,...,j_{e-1}] consecutive in Z/nZ.
    eps[s] in {0,1} for s=0..e-2 orients the edge (j_s, j_{s+1}):
       eps[s]=0  ->  a_{j_s} applied BEFORE a_{j_{s+1}}
       eps[s]=1  ->  a_{j_{s+1}} applied BEFORE a_{j_s}
    Returns one linear extension (application order).  Any two linear
    extensions of the same orientation differ by swaps of commuting
    (non-adjacent) generators, hence give the same element of A_n."""
    e = len(J)
    succ = {s: [] for s in range(e)}
    indeg = {s: 0 for s in range(e)}
    for s in range(e - 1):
        if eps[s] == 0:
            succ[s].append(s + 1); indeg[s + 1] += 1
        else:
            succ[s + 1].append(s); indeg[s] += 1
    order, avail = [], [s for s in range(e) if indeg[s] == 0]
    while avail:
        s = avail.pop(0)
        order.append(s)
        for u in succ[s]:
            indeg[u] -= 1
            if indeg[u] == 0:
                avail.append(u)
    assert len(order) == e
    return [J[s] for s in order]


def R_orientation_sum(e, n, k):
    """  R_e(t) := sum over intervals J of length e, sum over orientations eps,
         t^{|eps|} a_{J,eps}   in A_n[t].  """
    M = defaultdict(dict)
    for i in range(n):
        J = interval(i, e, n)
        for eps in _bits(e - 1):
            w = word_from_orientation(J, eps)
            add_op(M, op_from_word(w, n, k), scalar=t ** sum(eps))
    return clean(M)


def _bits(m):
    if m == 0:
        yield ()
        return
    for x in range(2 ** m):
        yield tuple((x >> s) & 1 for s in range(m))


def cyclic_ribbon_adder(e, n, k):
    """Independent definition: move one bead forward by e steps on Z/nZ,
       weight t^{#beads strictly jumped over} * q^{#seam crossings}."""
    M = defaultdict(dict)
    for S in states(n, k):
        for j in S:
            tgt = (j + e) % n
            if tgt in S:
                continue
            h = sum(1 for s in range(1, e) if (j + s) % n in S)
            J = interval(j, e, n)
            wrap = 1 if (n - 1) in J else 0
            T = frozenset((set(S) - {j}) | {tgt})
            M[S][T] = sp.expand(M[S].get(T, 0) + t ** h * q ** wrap)
    return clean(M)


# ---------------------------------------------------------------
# General subsets: runs, orientations, and the operators E, H, Q
# ---------------------------------------------------------------

def runs_of(I, n):
    """maximal cyclic runs of a PROPER subset I of Z/nZ."""
    I = set(I)
    assert len(I) < n
    starts = [i for i in I if (i - 1) % n not in I]
    out = []
    for s in starts:
        r = [s]
        while (r[-1] + 1) % n in I:
            r.append((r[-1] + 1) % n)
        out.append(r)
    return out


def word_from_subset_orientation(I, n, eps):
    """eps: tuple of bits, one per internal edge, listed run by run (runs in the
    order returned by runs_of, edges left to right inside each run)."""
    R = runs_of(I, n)
    w, p = [], 0
    for r in R:
        m = len(r) - 1
        w += word_from_orientation(r, eps[p:p + m])
        p += m
    assert p == len(eps)
    return w            # runs are non-adjacent, so they commute: order irrelevant


def n_internal_edges(I, n):
    return len(I) - len(runs_of(I, n))


def proper_subsets(n):
    for r in range(0, n):
        for I in combinations(range(n), r):
            yield set(I)


def Q_op(n, k, z):
    """Q(z;t) = sum over proper subsets I, over orientations eps,
       t^{|eps|} z^{|I|} a_{I,eps}.   (I = empty gives 1.)"""
    M = defaultdict(dict)
    for S in states(n, k):
        M[S][S] = sp.Integer(1)
    for I in proper_subsets(n):
        if not I:
            continue
        m = n_internal_edges(I, n)
        for eps in _bits(m):
            w = word_from_subset_orientation(I, n, eps)
            add_op(M, op_from_word(w, n, k), scalar=t ** sum(eps) * z ** len(I))
    return clean(M)


def EH_op(n, k, z, which):
    """which='e' -> E(z) = 1 + sum_i e_i z^i ; which='h' -> H(z)."""
    M = defaultdict(dict)
    for S in states(n, k):
        M[S][S] = sp.Integer(1)
    for I in proper_subsets(n):
        if not I:
            continue
        m = n_internal_edges(I, n)
        eps = tuple(1 for _ in range(m)) if which == 'e' else tuple(0 for _ in range(m))
        w = word_from_subset_orientation(I, n, eps)
        add_op(M, op_from_word(w, n, k), scalar=z ** len(I))
    return clean(M)


def strip_adder(n, k, z, kind):
    """Independent check of E,H: kind='h' adds cylindric horizontal strips,
       'e' vertical, weight z^{size} q^{seam crossings}. Built from the
       ribbon description: disjoint non-adjacent ribbons of height 0 (resp max)."""
    raise NotImplementedError
