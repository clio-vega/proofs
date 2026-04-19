"""
Even-block gap analysis for H-invariant theorem.

The question: In the staircase symmetrizing product for S_n,
the gap occurs at even block starts k >= 4. By parabolic reduction,
even block k reduces to checking all irreps of S_{k+1}.

The operator: A = (1+s_3)(1+s_1)(1+s_2) in S_4 (relabeled).

Key question: Is ker(A) = ker(1+s_2) in every irrep of S_4?
If yes, the gap closes because the transported image lives in E_{k-1}^{(+1)},
not E_{k-1}^{(-1)} = ker(1+s_{k-1}).
"""

from sage.combinat.tableau import StandardTableaux
from sage.algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra
import numpy as np

print("=" * 70)
print("EVEN-BLOCK GAP ANALYSIS FOR H-INVARIANT THEOREM")
print("=" * 70)
print()

def get_seminormal_rep(partition, n):
    """
    Compute Young's seminormal representation of S_n for partition lambda.
    Returns matrices for each simple transposition s_i = (i, i+1).

    Uses the standard seminormal form where for each adjacent pair of SYT T, T'
    (T' = s_i(T)), the matrix entries are determined by the axial distance.
    """
    STab = StandardTableaux(partition)
    tabs = list(STab)
    d = len(tabs)

    # Index tablets
    tab_idx = {T: i for i, T in enumerate(tabs)}

    matrices = {}
    for i in range(1, n):  # s_i = transposition (i, i+1)
        M = matrix(QQ, d, d)
        for idx, T in enumerate(tabs):
            # Find positions of i and i+1 in T
            pos_i = None
            pos_i1 = None
            for r, row in enumerate(T):
                for c, val in enumerate(row):
                    if val == i:
                        pos_i = (r, c)
                    elif val == i + 1:
                        pos_i1 = (r, c)

            # Axial distance: d = (col(i+1) - col(i)) - (row(i+1) - row(i))
            # Content: c(box) = col - row
            content_i = pos_i[1] - pos_i[0]
            content_i1 = pos_i1[1] - pos_i1[0]
            axial_dist = content_i1 - content_i

            # Diagonal entry: 1/d
            M[idx, idx] = QQ(1) / axial_dist

            # Try to apply s_i to T: swap i and i+1
            # Build new tableau
            new_tab_list = [list(row) for row in T]
            new_tab_list[pos_i[0]][pos_i[1]] = i + 1
            new_tab_list[pos_i1[0]][pos_i1[1]] = i

            # Check if this is a valid SYT
            try:
                new_T = Tableau(new_tab_list)
                if new_T.is_standard():
                    jdx = tab_idx[new_T]
                    # Off-diagonal: sqrt(1 - 1/d^2)
                    off_diag_sq = 1 - QQ(1) / (axial_dist ** 2)
                    if off_diag_sq > 0:
                        # Use exact square root if possible
                        off_diag = sqrt(QQ(off_diag_sq))
                        M[idx, jdx] = off_diag
            except:
                pass

        matrices[i] = M
    return matrices, tabs, d

def get_seminormal_rep_exact(partition, n):
    """
    Compute Young's seminormal representation using the exact rational form
    (without square roots - works over QQ directly using the approach where
    the basis is the set of SYT and matrices are explicit).

    Actually use the approach of working in QQ[sqrt(various)] or just use
    the QSym/SymmetricGroupAlgebra approach in Sage.
    """
    pass

def analyze_operator_S4():
    """
    Full analysis for S_4.
    The operator is A = (1+s_3)(1+s_1)(1+s_2).
    H = <s_1, s_3> (parabolic subgroup S_2 x S_2).
    """
    print("=" * 60)
    print("S_4 ANALYSIS")
    print("Operator: A = (1+s_3)(1+s_1)(1+s_2)")
    print("H = <s_1, s_3> ≅ S_2 × S_2")
    print("=" * 60)
    print()

    n = 4
    # Use Sage's symmetric group algebra
    S4 = SymmetricGroup(4)
    QS4 = S4.algebra(QQ)

    # Generators: s_i = (i, i+1) as permutations
    s1 = S4([(1,2)])
    s2 = S4([(2,3)])
    s3 = S4([(3,4)])

    e1 = QS4.one()
    S1 = QS4(s1)
    S2 = QS4(s2)
    S3 = QS4(s3)

    # Operator A = (1+s3)(1+s1)(1+s2)
    A_elem = (e1 + S3) * (e1 + S1) * (e1 + S2)

    # Also compute individual operators
    op_s2 = e1 + S2
    op_s1s2 = (e1 + S1) * (e1 + S2)
    op_A = A_elem

    print("Operator A expanded:", A_elem)
    print()

    # Analyze each irrep
    partitions_4 = Partitions(4).list()
    print(f"Irreps of S_4: {partitions_4}")
    print()

    results = {}

    for lam in partitions_4:
        lam_tuple = tuple(lam)
        print(f"--- Irrep V_{lam_tuple} (dim = {SemistandardTableaux(lam).cardinality() if False else StandardTableaux(lam).cardinality()}) ---")

        # Get the irrep via Specht module
        try:
            # Use Sage's built-in representation
            rho = S4.representation(lam, QQ)

            # Get matrices for generators
            basis_size = rho.dimension()
            print(f"  Dimension: {basis_size}")

            # Build matrix for A
            A_mat = rho.representation_matrix(A_elem)
            S2_mat = rho.representation_matrix(op_s2)
            S1S2_mat = rho.representation_matrix(op_s1s2)

            ker_A = A_mat.right_kernel()
            ker_s2 = S2_mat.right_kernel()
            ker_s1s2 = S1S2_mat.right_kernel()

            print(f"  dim ker(1+s_2) = {ker_s2.dimension()}")
            print(f"  dim ker((1+s_1)(1+s_2)) = {ker_s1s2.dimension()}")
            print(f"  dim ker(A) = {ker_A.dimension()}")

            # Are they equal?
            eq_A_s1s2 = (ker_A == ker_s1s2)
            eq_s1s2_s2 = (ker_s1s2 == ker_s2)
            eq_A_s2 = (ker_A == ker_s2)

            print(f"  ker(A) = ker((1+s_1)(1+s_2))? {eq_A_s1s2}")
            print(f"  ker((1+s_1)(1+s_2)) = ker(1+s_2)? {eq_s1s2_s2}")
            print(f"  ker(A) = ker(1+s_2)? {eq_A_s2}")

            # Compute V^H where H = <s1, s3>
            # V^H = vectors fixed by all h in H
            # = ker(s1 - 1) ∩ ker(s3 - 1)
            S1_mat = rho.representation_matrix(S1)
            S3_mat = rho.representation_matrix(S3)
            I = matrix.identity(QQ, basis_size)

            fixed_s1 = (S1_mat - I).right_kernel()
            fixed_s3 = (S3_mat - I).right_kernel()
            V_H = fixed_s1.intersection(fixed_s3)

            print(f"  dim V^H (fixed by <s1,s3>) = {V_H.dimension()}")

            # Compute the staircase image transport
            # Block 1: apply P_1 = (1+s_1)/2... wait, P_k is the symmetrizer
            # P_k = sum_{w in S_{k+1}/S_k} w ... no.
            # Actually P_k = (1 + s_k)(1 + s_{k-1})...(1+s_1) / (k+1)!/(k)! ...
            # Let me think. The staircase symmetrizer:
            # Π = P_{n-1} P_{n-2} ... P_1 where P_k acts on block ending at position k.
            # P_k = (1 + s_k)(1 + s_{k-1})...(1 + s_1) [unnormalized]
            #
            # For S_4 (n=4), the staircase is:
            # Block 1: P_1 = 1 + s_1
            # Block 2: P_2 * P_1 style... let me use the correct definition.
            #
            # From the H-invariant theorem context:
            # The staircase product Π = (1+s_1)(1+s_2)(1+s_1)(1+s_3)(1+s_2)(1+s_1)...
            # No - let me use the actual definition.
            #
            # The staircase is: for S_n,
            # Π = ∏_{k=1}^{n-1} (1+s_k)(1+s_{k-1})...(1+s_1)
            # where the product goes from k=1 to n-1.
            #
            # Transport of V^H at each stage:
            # After block 1 (P_1 = 1+s_1): image of V^H under (1+s_1)
            # After block 2 (P_2 P_1 = (1+s_2)(1+s_1) applied to above): etc.
            #
            # The "image at start of block 3" means after blocks 1 and 2.
            # Block 1: apply (1+s_1)
            # Block 2: apply (1+s_2)(1+s_1)
            # Wait - re-reading the problem:
            # "block k" corresponds to P_{k-2} applied after P_{k-1}"
            # The relabeling is: k-2->1, k-1->2, k->3
            # So we're checking the image at the "start of block k-1"
            # which in the S_4 language is "after blocks 1 and 2 have been applied"

            # Let me compute the image of V^H under the staircase operator up to block 2.
            # P_1 operator: (1+s_1)
            # P_2 operator: (1+s_2)(1+s_1) -- but is this applied to the IMAGE of P_1?
            # Or is it the composition?
            #
            # The staircase product Π = P_{n-1}...P_2 P_1 where P_k = (1+s_k)...
            # For a vector v in V^H:
            # After block 1: w_1 = (1+s_1)v
            # After block 2: w_2 = (1+s_2)(1+s_1) w_1 ...
            # No. The staircase is one big product:
            # Π = [(1+s_3)(1+s_2)(1+s_1)] [(1+s_2)(1+s_1)] [(1+s_1)]
            # Applied right to left:
            # First (1+s_1), then (1+s_2)(1+s_1), then (1+s_3)(1+s_2)(1+s_1)
            #
            # The "transported image at start of block 3" = image after the first two groups:
            # [(1+s_2)(1+s_1)] [(1+s_1)] V^H

            # Compute staircase through blocks 1 and 2
            P1 = I + S1_mat
            P2_full = (I + S2_mat) * (I + S1_mat)  # block 2 factor

            # Transport: apply P1 then P2_full
            staircase_12 = P2_full * P1  # combined operator through blocks 1 and 2

            # Image of V^H under staircase_12
            if V_H.dimension() > 0:
                VH_basis = V_H.basis_matrix()
                transported = VH_basis * staircase_12.transpose()  # rows are images
                transported_space = transported.row_space()

                print(f"  dim(transported V^H after blocks 1,2) = {transported_space.dimension()}")

                # Check intersection with ker(A)
                intersection = transported_space.intersection(ker_A)
                print(f"  dim(transported ∩ ker(A)) = {intersection.dimension()}")

                # The key question: is the transported image in E_2^{(+1)} = +1 eigenspace of s_2?
                # s_2 acts on transported image
                # E_2^{(+1)} = ker(s_2 - 1) = +1 eigenspace
                fixed_s2 = (S2_mat - I).right_kernel()
                anti_s2 = (S2_mat + I).right_kernel()  # -1 eigenspace = ker(1+s_2)

                print(f"  dim E_2^{{+1}} (fixed by s_2) = {fixed_s2.dimension()}")
                print(f"  dim E_2^{{-1}} (anti by s_2) = {anti_s2.dimension()}")

                trans_in_pos = transported_space.intersection(fixed_s2)
                trans_in_neg = transported_space.intersection(anti_s2)
                print(f"  transported ∩ E_2^{{+1}} dim = {trans_in_pos.dimension()}")
                print(f"  transported ∩ E_2^{{-1}} dim = {trans_in_neg.dimension()}")

                # P_{k-1} = (1+s_2) projects onto E_2^{(+1)}
                # So the image of P_{k-1} lives in E_2^{(+1)} = ker(s_2 - 1)
                # Hence transported image ⊆ E_2^{(+1)} iff P2 maps into +1 eigenspace of s_2
                # But (1+s_2) projects ONTO ker(s_2+1) complement...
                # Actually (1+s_2)v has s_2(1+s_2)v = (s_2 + s_2^2)v = (s_2 + 1)v = (1+s_2)v
                # So indeed (1+s_2)v is in the +1 eigenspace of s_2!
                # So the transported image DEFINITELY lives in E_2^{(+1)}, not E_2^{(-1)}.

                # But the transported space at "start of block 3" includes MORE than just P2 image.
                # Let me reconsider: the image is (1+s_2)(1+s_1)[(1+s_1)V^H]
                # The outermost factor is (1+s_2), so the image IS in the +1 eigenspace of s_2.

                print(f"  [CHECK] Is transported ⊆ E_2^{{+1}}? (should be yes due to (1+s2) factor)")
                is_subset = transported_space.is_subspace(fixed_s2)
                print(f"  transported ⊆ E_2^{{+1}}: {is_subset}")
            else:
                print(f"  V^H = 0 (trivial)")

            results[lam_tuple] = {
                'dim': basis_size,
                'ker_A': ker_A.dimension(),
                'ker_s2': ker_s2.dimension(),
                'ker_s1s2': ker_s1s2.dimension(),
                'ker_A_eq_ker_s2': eq_A_s2,
                'V_H_dim': V_H.dimension() if 'V_H' in dir() else 0
            }

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

        print()

    return results

def analyze_operator_Sn(n, even_blocks):
    """
    For S_n, analyze the even block gaps.
    even_blocks: list of k values (even block starts) to check.
    For each even block k, the operator in the relabeled S_{k+1} is:
      A = (1+s_k)(1+s_{k-2})(1+s_{k-1})
    But by parabolic reduction to S_{k+1}, this becomes A = (1+s_3)(1+s_1)(1+s_2).

    Actually we analyze directly in S_n for each even block.
    """
    print(f"\n{'='*60}")
    print(f"S_{n} ANALYSIS")
    print(f"{'='*60}")

    Sn = SymmetricGroup(n)
    QSn = Sn.algebra(QQ)

    gens = {}
    for i in range(1, n):
        gens[i] = Sn([(i, i+1)])

    e = QSn.one()
    S = {i: QSn(gens[i]) for i in range(1, n)}
    I_mat_cache = {}

    partitions_n = Partitions(n).list()

    for k in even_blocks:
        print(f"\n  Even block k={k} (gap at P_{{k-2}} step after P_{{k-1}} in block k-1)")
        print(f"  Operator: A = (1+s_{k})(1+s_{k-2})(1+s_{k-1})")

        # A = (1+s_k)(1+s_{k-2})(1+s_{k-1})
        A_elem = (e + S[k]) * (e + S[k-2]) * (e + S[k-1])
        op_sk = e + S[k-1]  # (1+s_{k-1})
        op_sk_prev = (e + S[k-2]) * (e + S[k-1])  # (1+s_{k-2})(1+s_{k-1})

        print(f"  Checking all {len(partitions_n)} irreps:")

        all_equal = True

        for lam in partitions_n:
            lam_tuple = tuple(lam)
            try:
                rho = Sn.representation(lam, QQ)
                basis_size = rho.dimension()

                A_mat = rho.representation_matrix(A_elem)
                Sk_mat = rho.representation_matrix(op_sk)
                Sk_prev_mat = rho.representation_matrix(op_sk_prev)

                ker_A = A_mat.right_kernel()
                ker_sk = Sk_mat.right_kernel()
                ker_sk_prev = Sk_prev_mat.right_kernel()

                eq_A_sk = (ker_A == ker_sk)
                eq_prev_sk = (ker_sk_prev == ker_sk)

                if not eq_A_sk:
                    all_equal = False

                status = "OK" if eq_A_sk else "GAP!"
                print(f"    λ={lam_tuple}: dim={basis_size}, "
                      f"ker(A)={ker_A.dimension()}, ker(1+s_{{k-1}})={ker_sk.dimension()}, "
                      f"ker(A)=ker(1+s_{{k-1}})? {eq_A_sk} [{status}]")

            except Exception as ex:
                print(f"    λ={lam_tuple}: ERROR - {ex}")

        if all_equal:
            print(f"  >>> ALL IRREPS: ker(A) = ker(1+s_{{k-1}}) for k={k}. Gap CLOSES!")
        else:
            print(f"  >>> SOME IRREPS: ker(A) ≠ ker(1+s_{{k-1}}) for k={k}. Gap remains!")

def analyze_transported_image_Sn(n):
    """
    For S_n, check the transported V^H image at each even block start.
    Verify that it lands in E_{k-1}^{(+1)} (not the -1 eigenspace).
    """
    print(f"\n{'='*60}")
    print(f"S_{n} TRANSPORTED IMAGE ANALYSIS")
    print(f"{'='*60}")

    Sn = SymmetricGroup(n)
    QSn = Sn.algebra(QQ)

    gens = {}
    for i in range(1, n):
        gens[i] = Sn([(i, i+1)])

    e = QSn.one()
    S = {i: QSn(gens[i]) for i in range(1, n)}

    # H = <s_1, s_3, s_5, ...> (non-adjacent generators)
    # For H-invariant theorem, H = parabolic subgroup generated by non-adjacent s_i
    # Standard choice: H = <s_1, s_3, s_5, ...>
    H_gens = list(range(1, n, 2))  # 1, 3, 5, ...
    print(f"  H = <s_{{{', s_'.join(str(i) for i in H_gens)}}}>")

    partitions_n = Partitions(n).list()

    for lam in partitions_n:
        lam_tuple = tuple(lam)
        try:
            rho = Sn.representation(lam, QQ)
            basis_size = rho.dimension()

            I = matrix.identity(QQ, basis_size)

            # Get all generator matrices
            S_mats = {i: rho.representation_matrix(QSn(gens[i])) for i in range(1, n)}

            # Compute V^H = intersection of +1 eigenspaces of all H generators
            V_H = matrix.identity(QQ, basis_size).row_space()  # Start with full space
            for i in H_gens:
                fixed_i = (S_mats[i] - I).right_kernel()
                V_H = V_H.intersection(fixed_i)

            if V_H.dimension() == 0:
                print(f"  λ={lam_tuple}: V^H = 0, skip")
                continue

            print(f"  λ={lam_tuple}: dim V^H = {V_H.dimension()}")

            # Transport V^H through the staircase
            # Staircase: Π = [B_{n-1}][B_{n-2}]...[B_1]
            # where B_k = (1+s_k)(1+s_{k-1})...(1+s_1)
            # Applied right to left: B_1 first, then B_2, etc.

            # Build the staircase block by block
            current_space = V_H

            for k in range(1, n):
                # Apply block k: (1+s_k)(1+s_{k-1})...(1+s_1)
                # Build the block operator matrix
                block_mat = I
                for j in range(1, k+1):
                    block_mat = (I + S_mats[j]) * block_mat  # apply s_j from left...
                    # Wait: (1+s_k)...(1+s_1) means s_1 applied first
                    # So we want (I + S_mats[k]) * ... * (I + S_mats[1])
                    # Let me rebuild:

                # Correct: block_k = (1+s_k)(1+s_{k-1})...(1+s_1)
                block_mat = I
                for j in range(k, 0, -1):
                    block_mat = (I + S_mats[j]) * block_mat

                if current_space.dimension() > 0:
                    basis = current_space.basis_matrix()
                    transported = basis * block_mat.transpose()
                    current_space = transported.row_space()

                # At even block starts (k+1 is even, i.e., k is odd), check eigenspace
                if k >= 2 and (k % 2 == 0):  # k even means next block k+1 is odd...
                    # Let me reconsider: "even block start k >= 4"
                    # After applying blocks 1 through k-1, we're at start of block k
                    pass

            # Simpler: just check the final transported image
            print(f"  Final transported dim: {current_space.dimension()}")

        except Exception as ex:
            print(f"  λ={lam_tuple}: ERROR - {ex}")
            import traceback
            traceback.print_exc()

def check_kernel_equality_S4_direct():
    """
    Direct, clean check: for each S_4 irrep, compute ker(A) and ker(1+s_2)
    and verify equality. Also verify the transported image lands in E_2^{(+1)}.
    """
    print("\n" + "="*60)
    print("DIRECT KERNEL EQUALITY CHECK: S_4")
    print("A = (1+s_3)(1+s_1)(1+s_2)")
    print("="*60 + "\n")

    S4 = SymmetricGroup(4)
    QS4 = S4.algebra(QQ)

    s1 = S4([(1,2)])
    s2 = S4([(2,3)])
    s3 = S4([(3,4)])

    e = QS4.one()
    S1, S2, S3 = QS4(s1), QS4(s2), QS4(s3)

    # The three operators
    op1 = e + S2                    # (1+s_2)
    op2 = (e + S1) * (e + S2)      # (1+s_1)(1+s_2)
    opA = (e + S3) * (e + S1) * (e + S2)  # A = (1+s_3)(1+s_1)(1+s_2)

    # H = <s_1, s_3> for the H-invariant

    summary = []

    for lam in Partitions(4):
        lam_t = tuple(lam)
        rho = S4.representation(lam, QQ)
        d = rho.dimension()
        I = matrix.identity(QQ, d)

        M1 = rho.representation_matrix(op1)
        M2 = rho.representation_matrix(op2)
        MA = rho.representation_matrix(opA)
        MS1 = rho.representation_matrix(S1)
        MS2 = rho.representation_matrix(S2)
        MS3 = rho.representation_matrix(S3)

        k1 = M1.right_kernel()
        k2 = M2.right_kernel()
        kA = MA.right_kernel()

        # V^H
        VH = (MS1 - I).right_kernel().intersection((MS3 - I).right_kernel())

        # Transported image: apply staircase blocks 1 and 2 to V^H
        # Block 1: (1+s_1)
        # Block 2: (1+s_2)(1+s_1)
        # Combined: (1+s_2)(1+s_1) * (1+s_1) * V^H
        # Wait: the staircase up to block 2 (excluding block 3) is:
        # [(1+s_2)(1+s_1)] then [(1+s_1)]
        # Applied to V^H: first (1+s_1)V^H, then (1+s_2)(1+s_1) applied to result

        P1_mat = I + MS1
        P2_mat = (I + MS2) * (I + MS1)  # block 2 operator

        # Staircase through blocks 1 and 2:
        staircase_mat = P2_mat * P1_mat

        if VH.dimension() > 0:
            VH_basis = VH.basis_matrix()
            transported_basis = VH_basis * staircase_mat.transpose()
            transported = transported_basis.row_space()
        else:
            transported = VH  # 0-dimensional

        # +1 eigenspace of s_2
        E2_plus = (MS2 - I).right_kernel()
        # -1 eigenspace of s_2
        E2_minus = (MS2 + I).right_kernel()

        trans_in_plus = transported.intersection(E2_plus)
        trans_in_minus = transported.intersection(E2_minus)
        trans_in_kerA = transported.intersection(kA)

        # Check: does (1+s_2) kill the transported image?
        # (1+s_2)w = 0 iff w in E2_minus
        # We need to check if transported ⊆ E2_plus (complement of E2_minus)

        row = {
            'lambda': lam_t,
            'dim': d,
            'dim_VH': VH.dimension(),
            'dim_transported': transported.dimension(),
            'ker(1+s2)_dim': k1.dimension(),
            'ker(op2)_dim': k2.dimension(),
            'ker(A)_dim': kA.dimension(),
            'ker(A)=ker(1+s2)': kA == k1,
            'ker(op2)=ker(1+s2)': k2 == k1,
            'transported⊆E2+': transported.is_subspace(E2_plus) if transported.dimension() > 0 else True,
            'transported∩ker(A)_dim': trans_in_kerA.dimension(),
        }
        summary.append(row)

        print(f"λ = {lam_t}  (dim={d})")
        print(f"  dim V^H = {VH.dimension()}")
        print(f"  dim transported image = {transported.dimension()}")
        print(f"  ker(1+s_2) dim = {k1.dimension()}")
        print(f"  ker((1+s_1)(1+s_2)) dim = {k2.dimension()}")
        print(f"  ker(A) dim = {kA.dimension()}")
        print(f"  ker(A) = ker(1+s_2)? {kA == k1}")
        print(f"  ker((1+s_1)(1+s_2)) = ker(1+s_2)? {k2 == k1}")
        print(f"  transported ⊆ E_2^{{+1}}? {transported.is_subspace(E2_plus) if transported.dimension() > 0 else 'vacuous'}")
        print(f"  transported ∩ ker(A) dim = {trans_in_kerA.dimension()}")
        if trans_in_kerA.dimension() > 0:
            print(f"  *** POTENTIAL GAP: transported image meets ker(A)! ***")
        else:
            print(f"  Gap CLOSED for this irrep (no intersection with ker(A))")
        print()

    return summary

# Run the main S_4 analysis
print("\n### PART 1: S_4 Direct Kernel Analysis ###\n")
summary = check_kernel_equality_S4_direct()

# Run S_n for n = 5, 6, 7
print("\n### PART 2: ker(A) = ker(1+s_{k-1}) for S_5, S_6, S_7 ###\n")
analyze_operator_Sn(5, [4])
analyze_operator_Sn(6, [4, 6])
analyze_operator_Sn(7, [4, 6])

print("\n### SUMMARY ###\n")
print("Key result: If ker(A) = ker(1+s_{k-1}) in every irrep,")
print("then the even-block gap closes because the transported image")
print("lies in E_{k-1}^{(+1)} (image of (1+s_{k-1})), not E_{k-1}^{(-1)}.")
print()
print("If also transported ⊆ E_{k-1}^{(+1)} and ker(A) ∩ E_{k-1}^{(+1)} = {0},")
print("then the within-block injectivity argument applies and the gap is closed.")
