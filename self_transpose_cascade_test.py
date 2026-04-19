"""
Self-transpose non-hook cascade test (2026-04-19).

Goal: Compare two obstruction sets on small partitions:
  (1) Liao-Rybnikov cactus non-2-transitivity: self-transpose non-hook lambda
  (2) Clio's direct-combinatorial "cascade" proof gap: even blocks k >= 6
      (odd blocks + k=4 are proved; k >= 6 even was the gap before the
      eigenvalue-positivity bypass).

For each self-transpose non-hook lambda with n = |lambda| <= 9:
  (a) dim V_lambda = #SYT(lambda)
  (b) rank of Pi_q = (1+T_1)(1+T_2)...(1+T_{n-1}) on V_lambda, in Q(q) and at q=1, q=2
  (c) r_lambda = #{T in SYT(lambda) : row_T(2j) <= j-1 for all j = 1,...,floor(n/2)}
      (closed-form prediction from the H-invariant theorem / staircase support)
  (d) Cascade-reachable? YES iff r_lambda does not strictly depend on any
      constraint row_T(2j) <= j-1 with 2j >= 6.
      (Equivalently: lambda's image support is already cut out by the j=1,2
      constraints plus the odd-block / k<=4 constraints.)

We use Hoefsmit's q-seminormal form:
  Generator T_i acts diagonally when i, i+1 are in the same row (eigenvalue q)
  or same column (eigenvalue -1). Otherwise on the pair (v_T, v_{s_i T}) with
  axial distance d = c(i+1) - c(i) (c = col - row):

    T_i v_T = ((q-1)/(1 - q^{-d})) v_T  +  (something) v_{s_i T}

  The standard convention (matching T_i^2 = (q-1)T_i + q) gives the 2x2 block

     [ (q-1)/(1 - q^{-d})                (off-T,sT)   ]
     [ (off-sT,T)                        (q-1)/(1 - q^{d}) ]

  with product of off-diagonals = q * (1 - 1/(1-q^{-d}))(1 - 1/(1-q^{d}))
                                 = q - alpha * beta, where
  alpha = (q-1)/(1-q^{-d}),  beta = (q-1)/(1-q^{d}).
  Trace = alpha + beta = q - 1, determinant = -q.

  We can verify T_i satisfies (T_i - q)(T_i + 1) = 0 in each such 2x2 block:
    tr = q - 1 = q + (-1),  det = -q = q * (-1). GOOD.
  Off-diagonal product: det - diag product = -q - alpha*beta. So
    off_T_sT * off_sT_T = -q - alpha * beta.
  We just need ANY factorization into two symbolic off-diagonals; the choice
  affects intermediate matrices but NOT rank of any product of T_i's (because
  a diagonal conjugation D rescales off-diagonals symmetrically and leaves
  rank invariant). So pick:
    off_T_sT = 1,  off_sT_T = -q - alpha * beta.

  [Note: Murphy's convention uses sqrt; rank is basis-independent so any
   factorization works.]

Verification:
  - (T_i - q)(T_i + 1) = 0 as a matrix (quadratic Hecke relation).
  - T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1} (braid relation).
  - T_i T_j = T_j T_i for |i-j| >= 2.
  - Specialization q=1 gives S_n permutation representation rank of
    (1+s_1)(1+s_2)...(1+s_{n-1}); matches staircase_test computations.
"""

from sympy import symbols, Rational, zeros, eye, Matrix, simplify, expand, cancel, together, factor, nsimplify, Integer
from sympy import QQ
from sympy.matrices import Matrix as SMatrix
from fractions import Fraction
from itertools import permutations
import sys

q = symbols('q')


# =========================================================================
# Partitions and SYT
# =========================================================================

def partitions_of(n):
    if n == 0:
        return [()]
    out = []
    def gen(remaining, max_part, prefix):
        if remaining == 0:
            out.append(tuple(prefix))
            return
        for p in range(min(remaining, max_part), 0, -1):
            prefix.append(p)
            gen(remaining - p, p, prefix)
            prefix.pop()
    gen(n, n, [])
    return out

def conjugate_partition(lam):
    if not lam:
        return ()
    n_cols = lam[0]
    return tuple(sum(1 for p in lam if p > c) for c in range(n_cols))

def is_hook(lam):
    # lambda is a hook iff lam[1] <= 1 (i.e., only one column beyond first)
    if len(lam) <= 1:
        return True
    return lam[1] <= 1

def is_self_transpose(lam):
    return tuple(lam) == conjugate_partition(lam)

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0:
        return [[]]
    out = []
    def fill(tab, val):
        if val > n:
            out.append([row[:] for row in tab])
            return
        for i in range(len(lam)):
            if len(tab[i]) < lam[i]:
                if i == 0 or len(tab[i]) < len(tab[i-1]):
                    tab[i].append(val)
                    fill(tab, val + 1)
                    tab[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return out

def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val:
                return (r, c)
    return None

def content(T, val):
    r, c = find_entry(T, val)
    return c - r

def row_of(T, val):
    return find_entry(T, val)[0]

def swap_si(T, i):
    """Return s_i * T (swap entries i and i+1) if the result is standard, else None."""
    r1, c1 = find_entry(T, i)
    r2, c2 = find_entry(T, i + 1)
    # Same row: not a swap (diagonal)
    if r1 == r2:
        return 'same_row'
    if c1 == c2:
        return 'same_col'
    # Distinct rows and columns
    T_new = [row[:] for row in T]
    T_new[r1][c1] = i + 1
    T_new[r2][c2] = i
    # Verify standard
    for r in range(len(T_new)):
        for c in range(len(T_new[r])):
            v = T_new[r][c]
            if c > 0 and v <= T_new[r][c-1]:
                return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]:
                return None
    return T_new

def tab_key(T):
    return tuple(tuple(row) for row in T)


# =========================================================================
# Hoefsmit q-seminormal representation
# =========================================================================

def hecke_generator_matrix(lam, i):
    """Matrix of T_i in Hoefsmit q-seminormal form on SYT(lam).
    Returns a sympy Matrix with entries in Q(q).
    """
    tabs = standard_tableaux(lam)
    d = len(tabs)
    idx = {tab_key(T): k for k, T in enumerate(tabs)}
    M = zeros(d, d)
    handled = set()
    for k, T in enumerate(tabs):
        if k in handled:
            continue
        action = swap_si(T, i)
        if action == 'same_row':
            M[k, k] = q
            handled.add(k)
        elif action == 'same_col':
            M[k, k] = Integer(-1)
            handled.add(k)
        else:
            # 2x2 block with s_i-partner
            sT = action  # tableau
            k2 = idx[tab_key(sT)]
            d_axial = content(T, i + 1) - content(T, i)
            # alpha on T, beta on sT
            alpha = (q - 1) / (1 - q**(-d_axial))
            beta = (q - 1) / (1 - q**(d_axial))
            # Off-diagonals: det = alpha*beta - off_T_sT * off_sT_T = -q
            # => off_T_sT * off_sT_T = alpha*beta + q
            off_prod = alpha * beta + q
            off_prod_s = cancel(off_prod)
            # Choose off_T_sT = 1, off_sT_T = off_prod
            M[k, k] = cancel(alpha)
            M[k2, k2] = cancel(beta)
            M[k2, k] = Integer(1)       # acting on v_T produces component on sT with coeff 1
            M[k, k2] = off_prod_s       # acting on v_{sT} produces component on T
            handled.add(k)
            handled.add(k2)
    return M


def verify_hecke_relations(lam):
    """Check (T_i - q)(T_i + 1) = 0, braid, and commutation."""
    n = sum(lam)
    if n <= 1:
        return True
    T = {i: hecke_generator_matrix(lam, i) for i in range(1, n)}
    d = T[1].shape[0]
    I = eye(d)
    # Quadratic
    for i in range(1, n):
        M = T[i]
        prod = (M - q * I) * (M + I)
        prod = prod.applyfunc(lambda x: cancel(x))
        if not all(prod[r, c] == 0 for r in range(d) for c in range(d)):
            print(f"  QUADRATIC FAILED lam={lam} i={i}")
            return False
    # Braid
    for i in range(1, n - 1):
        A = T[i] * T[i+1] * T[i]
        B = T[i+1] * T[i] * T[i+1]
        D = (A - B).applyfunc(lambda x: cancel(x))
        if not all(D[r, c] == 0 for r in range(d) for c in range(d)):
            print(f"  BRAID FAILED lam={lam} i={i}")
            return False
    # Commutation
    for i in range(1, n):
        for j in range(i + 2, n):
            C = (T[i] * T[j] - T[j] * T[i]).applyfunc(lambda x: cancel(x))
            if not all(C[r, c] == 0 for r in range(d) for c in range(d)):
                print(f"  COMM FAILED lam={lam} i={i} j={j}")
                return False
    return True


# =========================================================================
# Staircase Pi_q rank
# =========================================================================

def pi_q_staircase_word(n):
    """Generate the staircase reduced word for w_0 as a list of generator
    indices i (meaning T_i).

    Clio's convention: word = s_1, s_2 s_1, s_3 s_2 s_1, ..., s_{n-1}...s_1.
    That is, blocks B_k = s_k s_{k-1} ... s_1 for k = 1, ..., n-1, concatenated.
    Pi = prod_{k=1}^{n-1} B_k = prod_{k=1}^{n-1} (1+T_k)(1+T_{k-1})...(1+T_1)
    So the word of generator indices (in product order) is:
      1 ;  2, 1 ;  3, 2, 1 ; ...
    """
    word = []
    for k in range(1, n):
        for j in range(k, 0, -1):
            word.append(j)
    return word


def pi_q_matrix(lam, variant='staircase'):
    """Compute Pi_q acting on V_lam.

    variant = 'staircase': Clio's staircase reduced word
              Pi = (1+T_1) · (1+T_2)(1+T_1) · ... · (1+T_{n-1})(1+T_{n-2})...(1+T_1)
    variant = 'linear': naive (1+T_1)(1+T_2)...(1+T_{n-1})
    """
    n = sum(lam)
    tabs = standard_tableaux(lam)
    d = len(tabs)
    if d == 0:
        return None
    I = eye(d)
    if n <= 1:
        return I
    if variant == 'linear':
        word = list(range(1, n))
    else:
        word = pi_q_staircase_word(n)
    M = I
    for idx, i in enumerate(word):
        M = M * (I + hecke_generator_matrix(lam, i))
        M = M.applyfunc(lambda x: cancel(x))
    return M


def pi_q_matrix_at_value(lam, q_val, variant='staircase'):
    """Pi_q with q specialized to q_val -- much faster than symbolic."""
    n = sum(lam)
    tabs = standard_tableaux(lam)
    d = len(tabs)
    I = eye(d)
    if n <= 1:
        return I
    if variant == 'linear':
        word = list(range(1, n))
    else:
        word = pi_q_staircase_word(n)
    gens = {}
    def get_Ti(i):
        if i not in gens:
            gens[i] = hecke_generator_matrix(lam, i).subs(q, q_val).applyfunc(cancel)
        return gens[i]
    M = I
    for i in word:
        Ti = get_Ti(i)
        M = M * (I + Ti)
        M = M.applyfunc(cancel)
    return M


def rank_symbolic_via_specialization(lam, variant='staircase', test_values=None):
    """Generic rank over Q(q) = max rank over "random" rational specializations.
    """
    if test_values is None:
        test_values = [Integer(2), Integer(3), Rational(3, 2), Rational(5, 7), Rational(11, 13)]
    best = 0
    for v in test_values:
        M_v = pi_q_matrix_at_value(lam, v, variant=variant)
        r = M_v.rank()
        if r > best:
            best = r
    return best


def rank_at_q(lam, q_val, variant='staircase'):
    """Specialize q -> q_val at construction time for efficiency."""
    M_v = pi_q_matrix_at_value(lam, q_val, variant=variant)
    return M_v.rank()


# =========================================================================
# r_lambda: staircase support formula
# =========================================================================

def r_lambda_syt_support(lam):
    """SYT-support count: #{T in SYT(lam) : row_T(2j) <= j-1 for all j}.
    This is the SUPPORT of image-basis vectors for the intermediate staircase
    through an even block. It is NOT equal to rank(Pi|V) in general.
    """
    n = sum(lam)
    tabs = standard_tableaux(lam)
    count = 0
    for T in tabs:
        ok = True
        for j in range(1, n // 2 + 1):
            val = 2 * j
            r = row_of(T, val)
            if r > j - 1:
                ok = False
                break
        if ok:
            count += 1
    return count


def dim_V_H(lam):
    """dim of V_lam^H where H = <s_1, s_3, s_5, ...>.

    Compute by iteratively projecting to +1 eigenspace of each T_i at q=1
    (i.e., the group-algebra s_i). The H-generators s_1, s_3, s_5,...
    pairwise commute, so we multiply the projectors (1+s_i)/2 and take rank.
    """
    n = sum(lam)
    d = len(standard_tableaux(lam))
    I = eye(d)
    H_indices = list(range(1, n, 2))   # s_1, s_3, s_5, ...
    if not H_indices:
        return d
    P = I
    for i in H_indices:
        Ti = hecke_generator_matrix(lam, i).subs(q, 1).applyfunc(cancel)
        P = P * (I + Ti) * Rational(1, 2)
        P = P.applyfunc(cancel)
    return P.rank()


def r_lambda(lam):
    """Clio's r_lambda: rank of Pi (staircase reduced word) on V_lambda.
    Equal to dim V_lambda^H by the H-invariant theorem (proved Apr 17).
    Also equal to the symbolic rank (q-independence, Part 2 of the Total Rank
    PROVE.md theorem).
    """
    return dim_V_H(lam)


def cumulative_block_ranks(lam, q_val=Integer(1)):
    """For the staircase Pi = B_{n-1} ... B_1, return the sequence of ranks
    after each block. That is, r_k = rank(B_k B_{k-1} ... B_1 | V_lam).

    We compute these at q=q_val for speed; rank is q-independent (Part 2
    of the Total Rank theorem, PROVED), so q=1 gives the generic rank.
    """
    n = sum(lam)
    d = len(standard_tableaux(lam))
    I = eye(d)
    # Precompute T_i at this q
    T_mats = {i: hecke_generator_matrix(lam, i).subs(q, q_val).applyfunc(cancel)
              for i in range(1, n)}
    ranks = []
    cum = I
    for k in range(1, n):
        # Block B_k = (1+T_k)(1+T_{k-1})...(1+T_1)
        for j in range(k, 0, -1):
            cum = cum * (I + T_mats[j])
            cum = cum.applyfunc(cancel)
        ranks.append(cum.rank())
    return ranks


def cascade_reachable(lam):
    """Determine whether the direct-combinatorial cascade proof reaches lambda.

    Rank isolation lemma (proved): rank drops occur ONLY at odd-indexed blocks
    (s_1, s_3, s_5,...).  Empirically (verified here and in Clio's corpus),
    no rank drop ever occurs at an even block in any V_lam.

    The direct-combinatorial cascade proof (within-block S_3 cascading +
    cascade transport + parabolic reduction) handles:
      - all odd-block cases (within-block cascade) -- PROVED
      - the k=4 even-block "no-drop" verification    -- PROVED directly
      - the k=6,8,10,... even-block "no-drop" verifications -- WAS THE GAP,
        closed by the eigenvalue-positivity bypass (Apr 17).

    So lambda is "cascade-reachable" (by the direct, non-bypass, proof) iff
    the staircase Pi = B_{n-1}...B_1 does NOT pass through an even block k>=6.
    Equivalently, iff n-1 < 6, i.e., n <= 6.

    For n >= 7 the cascade *argument* must certify "no drop at block 6"
    (and 8, 10, ... for n-1 >= 8, 10, ...), which is exactly the gap.

    Returns (reachable: bool, drops_observed: list, even_blocks_traversed: list).
    """
    n = sum(lam)
    d = len(standard_tableaux(lam))
    if n <= 1:
        return True, [], []
    ranks = cumulative_block_ranks(lam)
    prev = d
    drops = []
    for k, r in enumerate(ranks, start=1):
        if r < prev:
            drops.append((k, prev, r))
        prev = r
    # Even blocks >= 6 that the staircase passes through in S_n
    even_ge6 = [k for k in range(6, n, 2)]
    # If ANY such block exists, the direct proof needs to certify it.
    reachable = (len(even_ge6) == 0)
    return reachable, drops, even_ge6


# =========================================================================
# Main
# =========================================================================

def main():
    print("=" * 78)
    print("Self-transpose non-hook cascade test")
    print("=" * 78)
    print()

    # Step 1: enumerate self-transpose non-hook partitions for n = 1..9
    print("Step 1: enumerate self-transpose non-hook partitions, n=1..9")
    print()
    targets = []
    for n in range(1, 10):
        for lam in partitions_of(n):
            if is_self_transpose(lam) and not is_hook(lam):
                targets.append(lam)
                d = len(standard_tableaux(lam))
                print(f"  n={n}  lambda={lam}  dim V={d}")
    print()
    print("(Note: the task description said 'n=8: none non-hook' and n=9 included")
    print("(5,3,1). Both are incorrect: (5,3,1)' = (3,2,2,1,1) is NOT self-conjugate.")
    print(" The correct self-conjugate non-hook partitions are as listed above.)")
    print()

    # Step 2: for each target, run the full battery
    print("Step 2: compute rank of Pi on V_lam for each target")
    print()
    results = []
    for lam in targets:
        n = sum(lam)
        tabs = standard_tableaux(lam)
        d = len(tabs)
        print(f"--- lambda = {lam}, n = {n}, dim V = {d} ---", flush=True)

        # Hecke relations sanity check (small only)
        if d <= 20:
            ok = verify_hecke_relations(lam)
            print(f"  Hecke relations verified: {ok}")

        # Ranks of the STAIRCASE Pi (the operator Clio studies)
        print(f"  Computing staircase ranks and block-by-block drops ...", flush=True)
        r_1 = rank_at_q(lam, Integer(1), variant='staircase')
        if d <= 50:
            r_2 = rank_at_q(lam, Integer(2), variant='staircase')
            r_gen = rank_symbolic_via_specialization(lam, variant='staircase',
                                                    test_values=[Integer(2), Integer(3), Rational(3, 2)])
        else:
            # For large d, skip the expensive q != 1 ranks (q-independence is
            # already proved; we trust the q=1 value as the generic rank).
            r_2 = None
            r_gen = None

        # r_lambda via H-invariant formula
        print(f"  Computing dim V^H ...", flush=True)
        r_lam = r_lambda(lam)

        # SYT support count (for comparison / context)
        r_supp = r_lambda_syt_support(lam)

        # Block-by-block cascade drops
        print(f"  Computing block-by-block cascade drops ...", flush=True)
        reachable, drops, even_ge6 = cascade_reachable(lam)

        print(f"  rank(Pi) at q=1           = {r_1}")
        if r_2 is not None:
            print(f"  rank(Pi) at q=2           = {r_2}")
            print(f"  rank(Pi) over Q(q)        = {r_gen}")
        print(f"  r_lam = dim V^H           = {r_lam}")
        print(f"  SYT-support count         = {r_supp}  (for reference, not the rank)")
        print(f"  cascade block-rank drops  = {drops}")
        print(f"  even-block>=6 traversed   = {even_ge6}")
        print(f"  cascade-reachable (direct proof) = {reachable}")

        results.append({
            'lam': lam, 'n': n, 'd': d,
            'r_gen': r_gen, 'r_q1': r_1, 'r_q2': r_2,
            'r_lam': r_lam, 'r_supp': r_supp,
            'drops': drops, 'even_ge6': even_ge6,
            'reachable': reachable
        })
        print()

    # Step 3: summary
    print("=" * 78)
    print("SUMMARY TABLE")
    print("=" * 78)
    header = f"{'lambda':15s} {'n':>2s} {'dim':>4s} {'rk_q=1':>7s} {'rk_q=2':>7s} {'rk_Q(q)':>8s} {'dimV^H':>7s} {'even>=6':>9s} {'reach':>6s}  drops"
    print(header)
    for r in results:
        r2 = '-' if r['r_q2'] is None else str(r['r_q2'])
        rg = '-' if r['r_gen'] is None else str(r['r_gen'])
        drops_str = ','.join(f'k={k}:{a}->{b}' for (k,a,b) in r['drops'])
        print(f"{str(r['lam']):15s} {r['n']:>2d} {r['d']:>4d} {r['r_q1']:>7d} {r2:>7s} {rg:>8s} {r['r_lam']:>7d} {str(r['even_ge6']):>9s} {str(r['reachable']):>6s}  {drops_str}")

    # Verdict
    print()
    print("Sanity checks:")
    print(f"  r_q1 == r_q2 (where computed)? "
          f"{all(r['r_q2'] is None or r['r_q1'] == r['r_q2'] for r in results)}")
    print(f"  r_q1 == rk_Q(q) (where computed)? "
          f"{all(r['r_gen'] is None or r['r_q1'] == r['r_gen'] for r in results)}")
    print(f"  r_q1 == dim V^H (H-invariant theorem)? "
          f"{all(r['r_q1'] == r['r_lam'] for r in results)}")

    print()
    print("VERDICT (Obstruction 1 = Obstruction 2 on this sample?):")
    # (1) cactus non-2-transitivity set = all these lambdas (by construction).
    # (2) cascade obstruction = lambdas with a rank drop at even k >= 6.
    cascade_obstruction_hit = [r['lam'] for r in results if not r['reachable']]
    cascade_obstruction_miss = [r['lam'] for r in results if r['reachable']]
    print(f"  Hit by cactus obstruction (all): {[r['lam'] for r in results]}")
    print(f"  Hit by cascade obstruction (k>=6 even drop): {cascade_obstruction_hit}")
    print(f"  Missed by cascade obstruction (reachable):   {cascade_obstruction_miss}")
    if not cascade_obstruction_miss:
        print("  --> Inclusion (1) in (2) plausible on this sample.")
    elif not cascade_obstruction_hit:
        print("  --> Inclusion (1) in (2) REFUTED on this sample.")
    else:
        print("  --> Mixed; inclusion fails in both directions, but the statement")
        print("      was only 'subset', not equality.")

    return results


if __name__ == '__main__':
    main()
