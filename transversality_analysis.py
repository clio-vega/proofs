"""
Transversality Analysis: Why does im(Pi_{k-1}) avoid -1(s_k) for even k?

The SYT Support Characterization says:
  T ∈ S (survives) iff row_T(2j) ≤ j-1 for all j = 1,...,floor(k/2)

Key puzzle: S contains T with rho(T,k) = -1 (pure -1(s_k) eigenvectors),
yet im(Pi_{k-1}) ∩ -1(s_k) = {0} for all irreps and even k.

This script investigates:
1. Principal angles between im(Pi_{k-1}) and -1(s_k)
2. Projection of problematic SYT |T> onto im(Pi_{k-1})
3. Whether s_k+1 vs -1 components of image satisfy dominance
4. T' structure: if T ∈ S with rho(T,k)=-1, does T' ∉ S? (T' doesn't exist when rho=-1)
5. Block structure analysis of 2x2 s_k blocks
"""

import numpy as np
from itertools import combinations
import sys

# ====== SYT infrastructure (self-contained) ======

def partitions(n):
    def _gen(n, max_val):
        if n == 0:
            yield ()
            return
        for i in range(min(n, max_val), 0, -1):
            for rest in _gen(n - i, i):
                yield (i,) + rest
    yield from _gen(n, n)

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0:
        return [[]]
    results = []
    def fill(tableau, val):
        if val > n:
            results.append([list(row) for row in tableau])
            return
        for i in range(len(lam)):
            if len(tableau[i]) < lam[i]:
                if i == 0 or len(tableau[i]) < len(tableau[i-1]):
                    tableau[i].append(val)
                    fill(tableau, val + 1)
                    tableau[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return results

def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val:
                return (r, c)
    return None

def row_of(T, val):
    for r, row in enumerate(T):
        if val in row:
            return r
    return None

def axial_distance(T, k):
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    return (c2 - r2) - (c1 - r1)

def tab_key(T):
    return tuple(tuple(row) for row in T)

def swap_entries(T, k):
    T_new = [list(row) for row in T]
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    T_new[r1][c1] = k + 1
    T_new[r2][c2] = k
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c-1]:
                return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]:
                return None
    return T_new

def satisfies_support_condition(T, k):
    """Check row_T(2j) <= j-1 for all j=1,...,floor(k/2)."""
    for j in range(1, k // 2 + 1):
        val = 2 * j
        r = row_of(T, val)
        if r is None:
            continue
        if r > j - 1:
            return False
    return True

# ====== Seminormal representation matrices ======

def seminormal_matrix(lam, k):
    tabs = standard_tableaux(lam)
    d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
    M = np.zeros((d, d))
    for T in tabs:
        i = tab_idx[tab_key(T)]
        rho = axial_distance(T, k)
        M[i, i] = 1.0 / rho
        T_sw = swap_entries(T, k)
        if T_sw is not None:
            key = tab_key(T_sw)
            if key in tab_idx:
                j = tab_idx[key]
                M[i, j] = np.sqrt(1.0 - 1.0/rho**2)
    return M, tabs

def projection_plus(lam, k):
    """P_k = (I + s_k)/2, projects onto +1(s_k) eigenspace."""
    M, tabs = seminormal_matrix(lam, k)
    d = len(tabs)
    return (np.eye(d) + M) / 2.0, tabs

def eigenspace(M, eigenval, tol=1e-10):
    d = M.shape[0]
    A = M - eigenval * np.eye(d)
    U, s, Vt = np.linalg.svd(A)
    null_mask = s < tol
    if not np.any(null_mask):
        return np.zeros((d, 0))
    return Vt[null_mask].T

def image_basis(M, tol=1e-10):
    if M.shape[1] == 0:
        return np.zeros((M.shape[0], 0))
    U, s, Vt = np.linalg.svd(M, full_matrices=True)
    rank = np.sum(s > tol)
    if rank == 0:
        return np.zeros((M.shape[0], 0))
    return U[:, :rank]

def subspace_intersection_dim(A, B, tol=1e-8):
    if A.shape[1] == 0 or B.shape[1] == 0:
        return 0
    combined = np.hstack([A, -B])
    _, s, _ = np.linalg.svd(combined)
    return max(0, combined.shape[1] - np.sum(s > tol))

def principal_angles(A, B, tol=1e-10):
    """Compute principal angles between column spaces of A and B."""
    if A.shape[1] == 0 or B.shape[1] == 0:
        return np.array([])
    QA, _ = np.linalg.qr(A)
    QB, _ = np.linalg.qr(B)
    # Use only orthonormal columns
    QA = QA[:, :A.shape[1]]
    QB = QB[:, :B.shape[1]]
    M_cross = QA.T @ QB
    s = np.linalg.svd(M_cross, compute_uv=False)
    s = np.clip(s, -1, 1)
    angles = np.arccos(s)
    return np.sort(angles)

def project_onto_subspace(v, basis):
    """Project v onto column space of basis."""
    if basis.shape[1] == 0:
        return np.zeros_like(v), 0.0
    Q, _ = np.linalg.qr(basis)
    Q = Q[:, :basis.shape[1]]
    comp = Q @ (Q.T @ v)
    frac = np.linalg.norm(comp) / (np.linalg.norm(v) + 1e-14)
    return comp, frac

def staircase_product(lam, up_to_block):
    """Product of staircase projections through blocks 1,...,up_to_block."""
    tabs = standard_tableaux(lam)
    d = len(tabs)
    result = np.eye(d)
    for block in range(1, up_to_block + 1):
        for k_in_block in range(block, 0, -1):
            P, _ = projection_plus(lam, k_in_block)
            result = P @ result
    return result, tabs

# ====== Main analysis ======

def analyze_transversality(n, k):
    """
    Analyze transversality for S_n, staircase block k (even).
    Returns dict with all findings.
    """
    assert k % 2 == 0, "k must be even"
    assert k <= n, "k must be <= n"

    print(f"\n{'='*60}")
    print(f"S_{n}, k={k} (analyzing im(Pi_{{k-1}}) vs -1(s_k))")
    print(f"{'='*60}")

    summary = {'n': n, 'k': k, 'irreps': []}

    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue

        # Build Pi_{k-1} and its image
        Pi, _ = staircase_product(lam, k - 1)
        im = image_basis(Pi)
        rank_im = im.shape[1]

        if rank_im == 0:
            continue

        # Build s_k matrix and its -1 eigenspace
        # Need k < n for s_k to act on SYT of size n
        if k >= n:
            continue
        Sk, _ = seminormal_matrix(lam, k)
        minus1_sk = eigenspace(Sk, -1.0)
        plus1_sk = eigenspace(Sk, 1.0)
        dim_minus = minus1_sk.shape[1]
        dim_plus = plus1_sk.shape[1]

        # 1. Intersection dimension (should be 0)
        inter_dim = subspace_intersection_dim(im, minus1_sk)

        # 2. Principal angles
        angles = principal_angles(im, minus1_sk)

        # 3. Identify support set S and problematic T
        tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
        support_set = []      # T ∈ S (satisfies row condition)
        prob_tabs = []        # T ∈ S with rho(T,k) = -1

        for i, T in enumerate(tabs):
            in_S = satisfies_support_condition(T, k)
            if in_S:
                support_set.append(i)
                if k < n:  # s_k is only defined when k+1 <= n
                    rho = axial_distance(T, k)
                    if abs(rho + 1) < 0.01:  # rho = -1
                        prob_tabs.append(i)

        # 4. For each problematic T: project |T> onto im(Pi_{k-1})
        prob_projections = []
        for i in prob_tabs:
            e_i = np.zeros(d)
            e_i[i] = 1.0
            comp, frac = project_onto_subspace(e_i, im)
            T = tabs[i]
            # Find T' (swap k and k+1 in T)
            T_sw = swap_entries(T, k)
            T_prime_in_S = None
            T_prime_idx = None
            if T_sw is not None:
                T_prime_idx = tab_idx.get(tab_key(T_sw))
                if T_prime_idx is not None:
                    T_prime_in_S = satisfies_support_condition(T_sw, k)

            # rho(T,k) = -1 means s_k|T> = -|T>, so T' doesn't exist
            # (swap would create invalid tableau)
            prob_projections.append({
                'tab_idx': i,
                'T': tab_key(T),
                'proj_frac': frac,
                'T_prime_exists': T_sw is not None,
                'T_prime_in_S': T_prime_in_S,
            })

        # 5. For each v ∈ im(Pi_{k-1}), decompose into s_k eigenbasis
        #    and measure +1 vs -1 components
        plus_fracs = []
        minus_fracs = []
        for j in range(rank_im):
            v = im[:, j]
            _, frac_plus = project_onto_subspace(v, plus1_sk)
            _, frac_minus = project_onto_subspace(v, minus1_sk)
            plus_fracs.append(frac_plus)
            minus_fracs.append(frac_minus)

        # 6. Check: for T ∈ S with rho(T,k) ≠ -1 (2x2 block)
        #    examine the paired T' and whether T' ∈ S
        block_analysis = []
        for i, T in enumerate(tabs):
            if not satisfies_support_condition(T, k):
                continue
            rho = axial_distance(T, k)
            if abs(rho + 1) < 0.01:
                # Pure -1 eigenvector, no T'
                continue
            # 2x2 block: T has s_k partner T'
            T_sw = swap_entries(T, k)
            if T_sw is None:
                continue
            T_prime_idx = tab_idx.get(tab_key(T_sw))
            if T_prime_idx is None:
                continue
            T_prime_in_S = satisfies_support_condition(T_sw, k)
            # Axial distance determines 2x2 block structure
            # s_k on {|T>, |T'>} = [[1/rho, sqrt(1-1/rho^2)],
            #                        [sqrt(1-1/rho^2), -1/rho]]
            # -1 eigenvector of [[a, b],[b,-a]] where a=1/rho, b=sqrt(1-1/rho^2):
            # eigenvalues +-1, -1 eigenvect proportional to (b, -(1+a))
            a = 1.0 / rho
            b = np.sqrt(max(0, 1 - 1/rho**2))
            norm_minus1_ev = np.sqrt(b**2 + (1+a)**2)
            alpha_minus = b / norm_minus1_ev      # coeff of T
            beta_minus = -(1+a) / norm_minus1_ev  # coeff of T'
            block_analysis.append({
                'T_idx': i, 'T_prime_idx': T_prime_idx,
                'rho': rho,
                'T_prime_in_S': T_prime_in_S,
                'alpha_minus_coeff': alpha_minus,
                'beta_minus_coeff': beta_minus,
            })

        # Summary for this irrep
        irrep_data = {
            'lam': lam,
            'd': d,
            'rank_im': rank_im,
            'dim_minus1': dim_minus,
            'dim_plus1': dim_plus,
            'inter_dim': inter_dim,
            'min_principal_angle_deg': np.degrees(angles[0]) if len(angles) > 0 else None,
            'all_angles_deg': np.degrees(angles).tolist() if len(angles) > 0 else [],
            'support_set_size': len(support_set),
            'prob_tabs_count': len(prob_tabs),
            'prob_projections': prob_projections,
            'plus_fracs': plus_fracs,
            'minus_fracs': minus_fracs,
            'block_analysis': block_analysis,
        }
        summary['irreps'].append(irrep_data)

        # Print summary
        print(f"\nλ={lam}: dim={d}, rank(im)={rank_im}, dim(-1(s_k))={dim_minus}")
        print(f"  Intersection dim: {inter_dim}  (should be 0)")
        if len(angles) > 0:
            print(f"  Principal angles (deg): {[f'{a:.2f}' for a in np.degrees(angles)]}")
            print(f"  Min angle: {np.degrees(angles[0]):.4f}°  (bounded away from 0: {angles[0] > 0.01})")
        print(f"  |S|={len(support_set)}, |S∩rho=-1|={len(prob_tabs)}")

        if prob_projections:
            print(f"  Problematic T (rho(T,k)=-1, T∈S) projections onto im:")
            for pp in prob_projections:
                print(f"    T={pp['T']}: proj_frac={pp['proj_frac']:.6f}, T'_exists={pp['T_prime_exists']}")

        if plus_fracs:
            max_minus = max(minus_fracs)
            min_plus = min(plus_fracs)
            print(f"  Im vector fracs: max(-1 comp)={max_minus:.6f}, min(+1 comp)={min_plus:.6f}")
            dominated = all(m < p for m, p in zip(minus_fracs, plus_fracs))
            print(f"  +1 dominates -1 for all im basis vectors: {dominated}")

        if block_analysis:
            t_prime_outside_S = [b for b in block_analysis if not b['T_prime_in_S']]
            print(f"  2x2 blocks with T∈S: {len(block_analysis)}, T'∉S: {len(t_prime_outside_S)}")
            for b in block_analysis[:3]:  # Show first few
                print(f"    rho={b['rho']:.2f}: α={b['alpha_minus_coeff']:.4f}(T) + β={b['beta_minus_coeff']:.4f}(T'), T'∈S={b['T_prime_in_S']}")

    return summary

def analyze_restriction_pattern(n_vals, k_vals):
    """
    Cross-cutting analysis: for T ∈ S with rho(T,k) = -1, check if
    these T are actually in the image. Since rho(T,k)=-1 means s_k|T>=-|T>,
    the key question is: what is the image of Pi on the subspace spanned by T?
    """
    print("\n" + "="*60)
    print("RESTRICTION PATTERN: T∈S with ρ(T,k)=-1 vs im(Pi_{k-1})")
    print("="*60)

    for n, k in zip(n_vals, k_vals):
        if k >= n:
            continue
        for lam in partitions(n):
            tabs = standard_tableaux(lam)
            d = len(tabs)
            if d <= 1:
                continue

            Pi, _ = staircase_product(lam, k - 1)
            Sk, _ = seminormal_matrix(lam, k)

            tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
            prob_tabs = []
            for i, T in enumerate(tabs):
                if satisfies_support_condition(T, k):
                    rho = axial_distance(T, k)
                    if abs(rho + 1) < 0.01:
                        prob_tabs.append((i, T))

            for (i, T) in prob_tabs:
                # Pi * e_i = image of |T> under Pi
                e_i = np.zeros(d)
                e_i[i] = 1.0
                v = Pi @ e_i
                norm_v = np.linalg.norm(v)

                # The key: since T∈S satisfies the row condition, we expect
                # |T> appears in the support. But the IMAGE Pi*|T> is NOT the
                # same as the indicator that T∈S!
                # Pi*e_i gives the column i of Pi, which shows how |T> maps under Pi.
                # If norm(Pi*e_i) ≈ 0, then |T> is killed by Pi!

                print(f"n={n}, k={k}, λ={lam}: T={tab_key(T)}, ρ(T,k)={axial_distance(T,k):.0f}")
                print(f"  ||Pi * |T>|| = {norm_v:.8f}")
                if norm_v > 1e-8:
                    # Decompose Pi*|T> into s_k eigenstates
                    minus1_sk = eigenspace(Sk, -1.0)
                    plus1_sk = eigenspace(Sk, 1.0)
                    _, frac_m = project_onto_subspace(v, minus1_sk)
                    _, frac_p = project_onto_subspace(v, plus1_sk)
                    print(f"  Pi*|T> frac in -1(s_k): {frac_m:.8f}, +1(s_k): {frac_p:.8f}")

def identify_critical_mechanism(n, k):
    """
    Identify the precise mechanism: what happens to |T> (T∈S, ρ=-1) under Pi?

    Key insight to check: Pi_{k-1} acting on a -1(s_k) eigenvector... does it
    necessarily kill it? Let's trace through the last few projections.
    """
    print(f"\n{'='*60}")
    print(f"CRITICAL MECHANISM TRACE: n={n}, k={k}")
    print(f"{'='*60}")

    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue

        Pi_full, _ = staircase_product(lam, k - 1)
        Sk, _ = seminormal_matrix(lam, k)
        Sk1, _ = seminormal_matrix(lam, k-1)  # s_{k-1}
        Sk2, _ = seminormal_matrix(lam, k-2)  # s_{k-2}

        tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}

        if k >= n:
            continue
        prob_tabs = []
        for i, T in enumerate(tabs):
            if satisfies_support_condition(T, k):
                rho = axial_distance(T, k)
                if abs(rho + 1) < 0.01:
                    prob_tabs.append((i, T))

        if not prob_tabs:
            continue

        print(f"\nλ={lam}:")
        for (i, T) in prob_tabs:
            e_i = np.zeros(d)
            e_i[i] = 1.0

            # Trace through each block
            state = e_i.copy()
            print(f"  T={tab_key(T)}, ρ(T,k)={axial_distance(T,k):.0f}")
            print(f"  Initial norm: {np.linalg.norm(state):.6f}")

            # Apply blocks one by one
            for block in range(1, k):
                prev_state = state.copy()
                for ki in range(block, 0, -1):
                    P, _ = projection_plus(lam, ki)
                    state = P @ state
                # Check -1(s_k) component after each block
                minus1_sk = eigenspace(Sk, -1.0)
                norm_state = np.linalg.norm(state)
                if norm_state > 1e-12 and minus1_sk.shape[1] > 0:
                    _, frac_m = project_onto_subspace(state, minus1_sk)
                    print(f"  After block {block}: norm={norm_state:.6f}, frac_in_-1(s_k)={frac_m:.6f}")
                elif norm_state > 1e-12:
                    print(f"  After block {block}: norm={norm_state:.6f}")
                else:
                    print(f"  After block {block}: KILLED (norm={norm_state:.2e})")
                    break


def analyze_P_k2_kills_minus1(n, k):
    """
    Test specific hypothesis: P_{k-2} kills all -1(s_k) components
    that might survive after P_{k-1}.

    The intuition: P_{k-2} = (I + s_{k-2})/2. s_{k-2} and s_k satisfy
    the braid relation s_{k-2} s_k s_{k-2} = s_k s_{k-2} s_k (far apart: they commute).
    Actually s_{k-2} and s_k commute since |k-2 - k| = 2 > 1.

    So P_{k-2} commutes with s_k. Therefore P_{k-2} preserves the -1(s_k) eigenspace!
    This means P_{k-2} CANNOT kill the -1(s_k) component.

    Instead, let's check P_{k-1}: s_{k-1} and s_k satisfy braid relation.
    P_{k-1} = (I + s_{k-1})/2 does NOT commute with s_k.
    """
    print(f"\n{'='*60}")
    print(f"COMMUTATIVITY ANALYSIS: n={n}, k={k}")
    print(f"{'='*60}")

    for lam in partitions(n):
        tabs = standard_tableaux(lam)
        d = len(tabs)
        if d <= 1:
            continue

        Sk, _ = seminormal_matrix(lam, k)
        Sk1, _ = seminormal_matrix(lam, k-1)
        Pk1 = (np.eye(d) + Sk1) / 2  # P_{k-1}
        minus1_sk = eigenspace(Sk, -1.0)

        if minus1_sk.shape[1] == 0:
            continue

        # Does P_{k-1} preserve -1(s_k) eigenspace?
        # Check: does P_{k-1} map -1(s_k) into itself?
        im_P_on_minus1 = image_basis(Pk1 @ minus1_sk)
        inter_dim = subspace_intersection_dim(im_P_on_minus1, minus1_sk)

        # Measure: what fraction of P_{k-1} * v is still in -1(s_k) for v in -1(s_k)?
        fracs = []
        for j in range(minus1_sk.shape[1]):
            v = minus1_sk[:, j]
            Pv = Pk1 @ v
            norm_Pv = np.linalg.norm(Pv)
            if norm_Pv > 1e-10:
                _, frac_m = project_onto_subspace(Pv, minus1_sk)
                fracs.append(frac_m)

        print(f"\nλ={lam}: dim(-1(s_k))={minus1_sk.shape[1]}")
        print(f"  P_{{k-1}} maps -1(s_k) into intersection with -1(s_k): dim={inter_dim}")
        if fracs:
            print(f"  Fracs of P_{{k-1}}*v in -1(s_k): {[f'{f:.4f}' for f in fracs]}")
            print(f"  Max frac: {max(fracs):.6f}  (0 means P_{{k-1}} kills all -1(s_k) components)")


def main():
    print("TRANSVERSALITY ANALYSIS: im(Pi_{k-1}) vs -1(s_k)")
    print("Checking WHY the image avoids -1(s_k) eigenvectors")
    print()

    # Test cases: (n, k) pairs with even k, k < n
    # s_k acts on positions k and k+1, so need k+1 <= n, i.e., k <= n-1
    test_cases = [
        (5, 4),
        (6, 4),
        (7, 4),
        (7, 6),
        (8, 4),
        (8, 6),
        (9, 4),
        (9, 6),
        (9, 8),
    ]

    # Main transversality analysis
    for n, k in test_cases:
        analyze_transversality(n, k)

    # Restriction pattern: what happens to problematic T under Pi?
    print("\n\n" + "#"*70)
    print("# RESTRICTION PATTERN ANALYSIS")
    print("#"*70)
    analyze_restriction_pattern(
        [5, 6, 7, 8],
        [4, 4, 4, 4]
    )

    # Critical mechanism trace for small cases
    print("\n\n" + "#"*70)
    print("# CRITICAL MECHANISM TRACE (n=5, k=4)")
    print("#"*70)
    identify_critical_mechanism(5, 4)
    identify_critical_mechanism(6, 4)

    # Commutativity analysis
    print("\n\n" + "#"*70)
    print("# COMMUTATIVITY ANALYSIS")
    print("#"*70)
    for n, k in [(5, 4), (6, 4), (7, 4), (7, 6)]:
        analyze_P_k2_kills_minus1(n, k)

    # Final summary: collect all intersection dims
    print("\n\n" + "#"*70)
    print("# GLOBAL SUMMARY: Intersection dims and min principal angles")
    print("#"*70)
    print(f"{'n':>3} {'k':>3} {'λ':>15} {'d':>4} {'r':>4} {'dim∩':>6} {'minAngle':>10}")
    for n, k in test_cases:
        if k >= n:
            continue
        for lam in partitions(n):
            tabs = standard_tableaux(lam)
            d = len(tabs)
            if d <= 1:
                continue

            Pi, _ = staircase_product(lam, k - 1)
            im = image_basis(Pi)
            rank_im = im.shape[1]
            if rank_im == 0:
                continue

            Sk, _ = seminormal_matrix(lam, k)
            minus1_sk = eigenspace(Sk, -1.0)
            if minus1_sk.shape[1] == 0:
                continue

            inter_dim = subspace_intersection_dim(im, minus1_sk)
            angles = principal_angles(im, minus1_sk)
            min_angle = np.degrees(angles[0]) if len(angles) > 0 else float('nan')

            print(f"{n:>3} {k:>3} {str(lam):>15} {d:>4} {rank_im:>4} {inter_dim:>6} {min_angle:>10.4f}")


if __name__ == '__main__':
    main()
