"""
Q460 -- the capacity computation.

LEMMA (capacity = rank).  Suppose the IK pairwise weight admits a memoryless
sequential factorization  W_sigma = prod_{a=1}^N omega(state(sigma(1..a-1)), sigma(a))
with state space S, |S| = D.  Fix a cut a and complementary column sets C, C^c.
Let P range over prefixes (matchings of C to rows 1..a) and Q over suffixes
(matchings of C^c to rows a+1..N).  Then

     M[P,Q] := prod_{i<=a<j} Delta_{P_i,Q_j}(x_i,x_j;z)

satisfies rank(M) <= D.

Proof.  W_{P.Q} = f(P) * g(Q) * M[P,Q] where f collects the within-prefix
factors and g the within-suffix factors.  The factorization gives
W_{P.Q} = alpha(P) * beta(state(P), Q).  Hence
M[P,Q] = (alpha(P)/f(P)) * (beta(state(P),Q)/g(Q)) = sum_{s in S} u_s(P) v_s(Q),
a sum of D rank-one matrices.  QED

So rank(M) is a certified LOWER BOUND on any state space.  We compute it.

Stratum: S = (1,1,...,1), M = 2N columns.  Every S_j = 1, so no two rows can
share a column: only the TIE-FREE cases A/B/C of eq:IK-delta1 occur -- exactly
the cases validated 13/13 against the paper's printed appendix expansions.
Column split: C = odd columns, C^c = even columns (interleaved, to expose the
nesting/crossing structure rather than the all-disjoint degeneration).
"""
import sympy as sp, itertools, sys, random
sys.path.insert(0, '/home/clio/projects/code/q460')
from ikdelta import Delta, q
import ikdelta as ik

def matchings(cols):
    """All perfect matchings of a column list, as tuples of sorted arcs."""
    cols = sorted(cols)
    if not cols: 
        yield ()
        return
    first, rest = cols[0], cols[1:]
    for k in range(len(rest)):
        partner = rest[k]
        remaining = rest[:k] + rest[k+1:]
        for m in matchings(remaining):
            yield ((first, partner),) + m

def ordered_matchings(cols):
    """Matchings of cols assigned to labelled rows: all permutations of arcs."""
    out = []
    for m in matchings(cols):
        for perm in itertools.permutations(m):
            out.append(perm)
    return out

def build_M(N, a, xs, qval, zvals=None):
    """The cut matrix M[P,Q], evaluated at given numeric q and x's (exact)."""
    cols = list(range(1, 2*N+1))
    C  = [c for c in cols if c % 2 == 1]
    Cc = [c for c in cols if c % 2 == 0]
    # prefix uses a arcs -> needs 2a columns from C ; take the first 2a of C
    assert 2*a <= len(C) and 2*(N-a) <= len(Cc) + (len(C)-2*a)
    Cpre  = C[:2*a]
    Csuf  = sorted(Cc + C[2*a:])
    Ps = ordered_matchings(Cpre)
    Qs = ordered_matchings(Csuf)
    sub = {q: qval}
    for i, v in enumerate(xs, start=1):
        sub[sp.Symbol('x%d' % i)] = v
    rows = []
    for P in Ps:
        row = []
        for Q in Qs:
            val = sp.Integer(1)
            for i in range(a):
                xi = sp.Symbol('x%d' % (i+1))
                for j in range(N-a):
                    xj = sp.Symbol('x%d' % (a+j+1))
                    val *= Delta(P[i], Q[j], xi, xj)
            row.append(sp.together(val).subs(sub))
        rows.append(row)
    return sp.Matrix(rows), len(Ps), len(Qs), Cpre, Csuf

def rank_at(N, a, qval, seed=1, label=""):
    random.seed(seed)
    # x's: distinct exact rationals, avoiding the poles x_j = x_i, q x_i
    xs = [sp.Rational(random.randint(2, 400), random.randint(1, 37)) for _ in range(N)]
    while len(set(xs)) < N:
        xs = [sp.Rational(random.randint(2, 400), random.randint(1, 37)) for _ in range(N)]
    M, nP, nQ, Cpre, Csuf = build_M(N, a, xs, qval)
    r = M.rank()
    print("  N=%d cut a=%d  q=%-10s  M is %dx%d   rank = %d   %s"
          % (N, a, qval, nP, nQ, r, label))
    return r, nP, nQ

if __name__ == "__main__":
    print("STRATUM: S=(1^{2N}); tie-free cases A/B/C only (validated 13/13).")
    print("Column split: prefix = first 2a odd columns, suffix = the rest.\n")
    print("GENERIC q (q = 5/3, an exact rational, not a root of unity):")
    res = {}
    for (N, a) in [(2,1), (3,1), (4,2), (5,2), (6,3)]:
        res[(N,a)] = rank_at(N, a, sp.Rational(5,3))
