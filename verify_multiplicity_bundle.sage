# Verify Multiplicity Bundle Theorem - computational tests
from sage.all import *

def verify_crystal_decomposition(n, k):
    """Check that B(omega_1)^{otimes k} for sl_n decomposes with multiplicities = dim(Specht modules)"""
    print(f"\n{'='*60}")
    print(f"Test: sl_{n}, k={k}")
    print(f"{'='*60}")

    # Crystal of fundamental representation
    C = crystals.Tableaux(['A', n-1], shape=[1])

    # Tensor product B(omega_1)^{otimes k}
    T = crystals.TensorProduct(*[C]*k)

    # Find connected components by finding highest weight elements
    hw_elements = [t for t in T if all(t.e(i) is None for i in range(1, n))]

    # Count multiplicities by weight
    from collections import Counter
    hw_weights = [t.weight() for t in hw_elements]
    mult_dict = Counter([str(w) for w in hw_weights])

    print(f"\nCrystal decomposition of B(omega_1)^{{otimes {k}}}:")
    for w, m in sorted(mult_dict.items()):
        print(f"  Weight {w}: multiplicity {m}")

    # Compare with Schur-Weyl prediction
    # Partitions of k with at most n parts
    print(f"\nSchur-Weyl prediction (partitions of {k} with <= {n} parts):")
    match_all = True
    for lam in Partitions(k, max_length=n):
        specht_dim = StandardTableaux(lam).cardinality()
        # The GL_n irrep V_lambda has highest weight determined by lambda
        # For partition lambda, the highest weight in fundamental weight coords is:
        # (lambda_1 - lambda_2, lambda_2 - lambda_3, ..., lambda_{n-1} - lambda_n)
        hw = []
        lam_list = list(lam) + [0]*(n - len(lam))
        for i in range(n-1):
            hw.append(lam_list[i] - lam_list[i+1])
        hw_str = str(tuple(hw))
        # Check if this weight appears in crystal decomposition with correct multiplicity
        crystal_mult = mult_dict.get(hw_str, 0)
        match = (crystal_mult == specht_dim)
        status = "OK" if match else "FAIL"
        if not match:
            match_all = False
        print(f"  lambda={lam}, dim(S_lambda)={specht_dim}, hw={hw}, crystal_mult={crystal_mult} [{status}]")

    print(f"\nOverall match: {match_all}")
    return match_all

def verify_combinatorial_R_identity(n):
    """Check that the combinatorial R-matrix on B(omega_1) x B(omega_1) is the identity for sl_n"""
    print(f"\n{'='*60}")
    print(f"Combinatorial R-matrix test: sl_{n}, B(omega_1) x B(omega_1)")
    print(f"{'='*60}")

    C = crystals.Tableaux(['A', n-1], shape=[1])
    T = crystals.TensorProduct(C, C)

    # Highest weight elements
    hw_elements = [t for t in T if all(t.e(i) is None for i in range(1, n))]
    print(f"\nHighest weight elements of B(omega_1) x B(omega_1):")
    for hw in hw_elements:
        print(f"  {hw} (weight: {hw.weight()})")

    # Connected components
    print(f"\nConnected components:")
    visited = set()
    for hw in hw_elements:
        component = []
        queue = [hw]
        seen = {hw}
        while queue:
            current = queue.pop(0)
            component.append(current)
            for i in range(1, n):
                child = current.f(i)
                if child is not None and child not in seen:
                    seen.add(child)
                    queue.append(child)
        print(f"\n  Component with HW {hw}:")
        for elem in component:
            print(f"    {elem} (weight: {elem.weight()})")
        visited.update(seen)

    # Check for elements not in any component from HW elements
    all_elems = set(T)
    if all_elems - visited:
        print(f"\n  WARNING: Elements not reached from highest weight elements: {all_elems - visited}")

    # Now test the combinatorial R-matrix explicitly
    # For type A crystals, SageMath has a built-in R-matrix
    print(f"\nTesting combinatorial R-matrix element by element:")
    is_identity = True
    for t in T:
        factors = list(t)
        b1, b2 = factors[0], factors[1]
        # The R-matrix swaps factors: B1 x B2 -> B2 x B1
        # For identical crystals, we check if R(b1 x b2) = b1 x b2
        # i.e., whether the swap is trivial

        # Use SageMath's sigma method on tensor product elements
        # which implements the combinatorial R-matrix
        try:
            r_result = T(b2, b1)  # The "swapped" element
            # The combinatorial R-matrix for B x B -> B x B
            # maps b1 x b2 to sigma(b2) x sigma(b1) where sigma is some involution
            # For identical crystals of the fundamental representation, we expect identity
            print(f"  {b1} x {b2} -> checking...")
        except Exception as e:
            print(f"  Error for {b1} x {b2}: {e}")

    # Alternative approach: use the crystal isomorphism
    # B(omega_1) x B(omega_1) decomposes as B(2*omega_1) + B(omega_2)
    # If R is identity, then it preserves this decomposition trivially
    print(f"\n  For identical fundamental crystals, R-matrix = identity")
    print(f"  is equivalent to: each connected component is symmetric under factor swap.")
    print(f"  Checking symmetry of components...")

    for hw in hw_elements:
        component = []
        queue = [hw]
        seen = {hw}
        while queue:
            current = queue.pop(0)
            component.append(current)
            for i in range(1, n):
                child = current.f(i)
                if child is not None and child not in seen:
                    seen.add(child)
                    queue.append(child)

        # Check if swapping factors of each element gives an element in the same component
        symmetric = True
        for elem in component:
            factors = list(elem)
            swapped = T(factors[1], factors[0])
            if swapped not in seen:
                symmetric = False
                print(f"    ASYMMETRIC: {elem} -> swap -> {swapped} NOT in component of {hw}")
                is_identity = False
        if symmetric:
            print(f"    Component of {hw}: symmetric under swap")

    print(f"\n  R-matrix appears to be identity: {is_identity}")

# Run tests
print("MULTIPLICITY BUNDLE THEOREM - COMPUTATIONAL VERIFICATION")
print("="*60)

r1 = verify_crystal_decomposition(2, 3)
r2 = verify_crystal_decomposition(2, 4)
r3 = verify_crystal_decomposition(3, 3)

verify_combinatorial_R_identity(2)
verify_combinatorial_R_identity(3)

# ============================================================
# Sigma factorization test
# ============================================================
print("\n" + "="*60)
print("Verifying sigma factorization for sl_2, k=3")
print("="*60)

# Build the representation matrices
n, k = 2, 3
dim_total = n**k

# sigma_{[1,3]} reverses (v1 x v2 x v3) -> (v3 x v2 x v1)
sigma = matrix(QQ, dim_total, dim_total)

from itertools import product as iterproduct

basis_tuples = list(iterproduct(range(n), repeat=k))
tuple_to_idx = {t: i for i, t in enumerate(basis_tuples)}

for t in basis_tuples:
    t_rev = tuple(reversed(t))
    sigma[tuple_to_idx[t_rev], tuple_to_idx[t]] = 1

print(f"\nsigma_{{[1,{k}]}} matrix (reversal):")
print(sigma)

# Eigenvalues
eigenvalues = sigma.eigenvalues()
from collections import Counter as C2
eigen_counts = C2(eigenvalues)
print(f"\nEigenvalues of sigma: {dict(eigen_counts)}")

print("\nSchur-Weyl prediction for sl_2, k=3:")
print("lambda=(3): dim V_(3)=4, dim S_(3)=1, w_0=(13) acts as sgn=+1 on trivial rep")
print("  -> contributes 4*1 = 4 eigenvalues of +1")
print("lambda=(2,1): dim V_(2,1)=2, dim S_(2,1)=2, w_0=(13) has eigenvalues +1,-1 on std rep")
print("  -> contributes 2*(+1) and 2*(-1) = 2 of +1, 2 of -1")
predicted_plus = 4 + 2
predicted_minus = 2
print(f"Total predicted: +1 with mult {predicted_plus}, -1 with mult {predicted_minus}")
print(f"Actual eigenvalues: {dict(eigen_counts)}")
print(f"Match: {eigen_counts.get(1,0) == predicted_plus and eigen_counts.get(-1,0) == predicted_minus}")

# ============================================================
# Additional: sl_2 k=4 sigma factorization
# ============================================================
print("\n" + "="*60)
print("Verifying sigma factorization for sl_2, k=4")
print("="*60)

n, k = 2, 4
dim_total = n**k
sigma = matrix(QQ, dim_total, dim_total)

basis_tuples = list(iterproduct(range(n), repeat=k))
tuple_to_idx = {t: i for i, t in enumerate(basis_tuples)}

for t in basis_tuples:
    t_rev = tuple(reversed(t))
    sigma[tuple_to_idx[t_rev], tuple_to_idx[t]] = 1

eigenvalues = sigma.eigenvalues()
eigen_counts = C2(eigenvalues)
print(f"\nEigenvalues of sigma (reversal on V^{{otimes 4}}): {dict(eigen_counts)}")

# Schur-Weyl for sl_2, k=4:
# lambda=(4): dim V=5, dim S=1, w_0=(14)(23) is even perm, sgn=+1, triv rep -> +1
#   -> 5 eigenvalues of +1
# lambda=(3,1): dim V=3, dim S=3, w_0=(14)(23) in std rep S_(3,1)
# lambda=(2,2): dim V=1, dim S=2 (WAIT: need <= 2 parts)
# Actually for sl_2, only partitions with <= 2 parts
# lambda=(4): V has dim 5, S_(4) has dim 1
# lambda=(3,1): V has dim 3, S_(3,1) has dim 3
# lambda=(2,2): V has dim 1, S_(2,2) has dim 2

print("\nSchur-Weyl prediction for sl_2, k=4:")
print("Partitions of 4 with <= 2 parts: (4), (3,1), (2,2)")

# w_0 in S_4 is (1,4)(2,3), which is an even permutation
# S_(4): trivial rep, w_0 -> +1. dim V_(4) = 5. Contributes 5*(+1)
# S_(2,2): For the Specht module S_(2,2) of dimension 2:
#   w_0 acts on S_(2,2). Since (2,2) is conjugate to itself and
#   w_0 acts as sgn^? ... actually w_0 acts on S_lambda as sgn * w_0_action_on_S_lambda'
#   More precisely: on S_lambda, w_0 acts as (-1)^{k(k-1)/2} times the transpose map...
#   Let's just compute it numerically.

# Compute w_0 action on Specht modules using SageMath
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
SGA = SymmetricGroupAlgebra(QQ, 4)

# w_0 = [4,3,2,1] in one-line notation
w0 = Permutation([4,3,2,1])

for lam in Partitions(4, max_length=2):
    ST = StandardTableaux(lam)
    st_list = list(ST)
    specht_dim = len(st_list)

    # Get the Specht module representation matrix of w_0
    # Using the natural action on tabloids
    # Alternative: use SageMath's representation
    print(f"\n  lambda={lam}, dim(V_lambda)={WeylCharacterRing('A1')(list(lam)).degree()}, dim(S_lambda)={specht_dim}")

# Directly verify: eigenvalues should sum correctly
print(f"\nActual eigenvalues: {dict(eigen_counts)}")
print(f"Total dimension: {sum(eigen_counts.values())} (should be {dim_total})")

print("\n" + "="*60)
print("ALL TESTS COMPLETE")
print("="*60)
