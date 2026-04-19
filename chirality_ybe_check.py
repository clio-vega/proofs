"""
Chirality Obstruction Theorem — Computational verification.

We work with the six-vertex R-matrix R(u) on C^2 ⊗ C^2, then set b2=0
and examine which YBE equations fail.

Convention: R(u) in basis {|00>, |01>, |10>, |11>}:
  R(u) = diag-block with:
    |00> -> a1(u)|00>
    |01> -> b2(u)|01> + c2(u)|10>
    |10> -> c1(u)|10> + b1(u)|10>   ... wait, let me be careful.

Vertex weights: R^{kl}_{ij} = weight of vertex with left-in=i, bottom-in=j, right-out=k, top-out=l.
Ice rule: i+j = k+l.

As a matrix acting on V⊗V where the first tensor factor is "horizontal" and second is "vertical":
  R|ij> = sum_{k,l} R^{kl}_{ij} |kl>

The six vertices (i,j) -> (k,l):
  (0,0) -> (0,0): a1
  (1,1) -> (1,1): a2
  (1,0) -> (1,0): b1  (1 passes horizontally, 0 passes vertically)
  (0,1) -> (0,1): b2  (0 passes horizontally, 1 passes vertically)
  (1,0) -> (0,1): c1  (1 turns up, 0 turns right)
  (0,1) -> (1,0): c2  (0 turns up, 1 turns right)

Matrix in basis |00>, |01>, |10>, |11>:
  R = [[a1, 0,  0,  0 ],
       [0,  b2, c2, 0 ],
       [0,  c1, b1, 0 ],
       [0,  0,  0,  a2]]
"""

from sympy import symbols, expand, factor, Matrix, eye, zeros, simplify
from sympy import kronecker_product

u, v = symbols('u v')
q = symbols('q')

# ============================================================
# PART 1: Parametric six-vertex R-matrix
# ============================================================

# Use the standard trigonometric/rational parametrization.
# For the Hecke-type R-matrix: R(u) = u*T + I where T is the Hecke generator.
#
# From kl_ybe_test.py, the Hecke generator on V⊗V is:
#   T = [[q, 0, 0, 0],
#        [0, q-1, 1, 0],   <-- note: T|01> = (q-1)|01> + 1|10>
#        [0, q, 0, 0],     <-- note: T|10> = q|01> + 0|10>
#        [0, 0, 0, q]]
#
# Wait, from the KL test: the 2x2 block on {|01>,|10>} is [[q-1, q], [1, 0]].
# This means T|01> = (q-1)|01> + |10> and T|10> = q|01>.
# Hmm, let me re-read the block ordering.

# From kl_ybe_test.py Part 3:
#   R_hecke = [[q,   0,   0,   0],
#              [0,   q-1, q,   0],
#              [0,   1,   0,   0],
#              [0,   0,   0,   q]]
#
# So: T|01> = (q-1)|01> + q|10>, T|10> = |01>.
# Wait no: R_hecke[1,1]=q-1, R_hecke[1,2]=q means |01> -> (q-1)|01> + q|10>?
# No! Matrix acts as R|v> where v is a column vector.
# |01> is the second basis vector (index 1).
# R * |01> = column 1 of R = [0, q-1, 1, 0]^T = (q-1)|01> + 1|10>.
# R * |10> = column 2 of R = [0, q, 0, 0]^T = q|01>.
#
# So the 2x2 block acting on (|01>, |10>) with the COLUMN convention is:
#   [[q-1, q],    <- coefficients of |01>
#    [1,   0]]    <- coefficients of |10>
# Applied to |01>: (q-1)|01> + 1|10>.
# Applied to |10>: q|01> + 0|10>.

# This gives vertex weights:
#   T|01> = (q-1)|01> + 1|10>  =>  b2 = q-1, c2 = 1
#   T|10> = q|01> + 0|10>      =>  c1 = q, b1 = 0
# And a1 = a2 = q.
#
# Interesting: in the HECKE form, it's b1 = 0 that vanishes, not b2.
# The chirality is reversed compared to PROVE.md's convention.
#
# For the spectral-parameter version R(u) = u*T + I:
#   a1(u) = uq + 1
#   a2(u) = uq + 1
#   b2(u) = u(q-1) + 1
#   c2(u) = u
#   c1(u) = uq
#   b1(u) = 1
#
# ALL six weights are nonzero for generic u,q. This is the FULL six-vertex model.

# Now, for the FIVE-vertex (chiral) model, we want one of the b-weights to vanish.
# Let's work with general symbolic weights to keep things clean.

# General six-vertex R-matrix with spectral parameter
def R_general(a1, a2, b1, b2, c1, c2):
    """R-matrix with given weights."""
    return Matrix([
        [a1, 0,  0,  0 ],
        [0,  b2, c2, 0 ],
        [0,  c1, b1, 0 ],
        [0,  0,  0,  a2]
    ])

# ============================================================
# PART 2: YBE computation — general framework
# ============================================================

I2 = eye(2)
I4 = eye(4)

def compute_R12(R):
    """R acting on spaces 1,2 in V^⊗3."""
    return kronecker_product(R, I2)

def compute_R23(R):
    """R acting on spaces 2,3 in V^⊗3."""
    return kronecker_product(I2, R)

def compute_R13(R):
    """R acting on spaces 1,3 in V^⊗3. Uses P23 conjugation."""
    # P23 swaps spaces 2 and 3
    P = Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])  # swap on C^2⊗C^2
    P23 = kronecker_product(I2, P)
    R12 = kronecker_product(R, I2)
    return P23 * R12 * P23

# ============================================================
# PART 3: Five-vertex model — set b2 = 0
# ============================================================

# Use symbolic weights with spectral parameters u, v
# Standard parametrization for five-vertex:
#   a1(u), a2(u), b1(u), c1(u), c2(u) are nonzero
#   b2(u) = 0

# Let's use the most general form: independent symbolic weights for each spectral parameter.
a1u, a2u, b1u, c1u, c2u = symbols('a1u a2u b1u c1u c2u')
a1v, a2v, b1v, c1v, c2v = symbols('a1v a2v b1v c1v c2v')
a1w, a2w, b1w, c1w, c2w = symbols('a1w a2w b1w c1w c2w')

# The YBE in braid form: R12(u) R23(u+v) R12(v) = R23(v) R12(u+v) R23(u)
# But with general weights, we write it as:
# R12(p) R13(r) R23(s) = R23(s) R13(r) R12(p)
# where p, r, s are three potentially different sets of weights.
# (In the difference form: p=u-v, r=u, s=v.)

# FIVE-VERTEX: b2 = 0 for all spectral parameters
R_p = R_general(a1u, a2u, b1u, 0, c1u, c2u)  # R12 with parameter p=u-v
R_r = R_general(a1w, a2w, b1w, 0, c1w, c2w)  # R13 with parameter r=u
R_s = R_general(a1v, a2v, b1v, 0, c1v, c2v)  # R23 with parameter s=v

print("=" * 70)
print("FIVE-VERTEX YBE COMPUTATION (b2 = 0)")
print("=" * 70)

print("\nR-matrix (five-vertex, b2=0):")
print(R_p)

R12_p = compute_R12(R_p)
R13_r = compute_R13(R_r)
R23_s = compute_R23(R_s)

LHS = R12_p * R13_r * R23_s
RHS = R23_s * R13_r * R12_p

diff = LHS - RHS

# Expand each entry
print("\nComputing LHS - RHS...")
diff_expanded = diff.applyfunc(expand)

# Find nonzero entries
print("\nNonzero entries of LHS - RHS:")
basis_labels = ['000', '001', '010', '011', '100', '101', '110', '111']
nonzero_count = 0
nonzero_positions = []
for i in range(8):
    for j in range(8):
        if diff_expanded[i, j] != 0:
            nonzero_count += 1
            nonzero_positions.append((i, j))
            print(f"  ({basis_labels[i]}, {basis_labels[j]}): {factor(diff_expanded[i, j])}")

print(f"\nTotal nonzero entries: {nonzero_count}")
print(f"Positions: {[(basis_labels[i], basis_labels[j]) for i,j in nonzero_positions]}")

# ============================================================
# PART 4: Understand which boundary conditions fail
# ============================================================

print("\n" + "=" * 70)
print("BOUNDARY CONDITION ANALYSIS")
print("=" * 70)

# The YBE says: for all boundary (i,j,k) in and (l,m,n) out,
# sum over internal states of LHS = sum over internal states of RHS.
# In the matrix form, (LHS - RHS)[out_state, in_state] = 0.
#
# The basis is |ijk> where i,j,k ∈ {0,1}.
# A nonzero entry at position (out, in) means the YBE fails for
# that particular (in-state, out-state) pair.

print("\nFailing boundary configurations (in -> out):")
for i, j in nonzero_positions:
    in_state = basis_labels[j]
    out_state = basis_labels[i]
    print(f"  |{in_state}> -> |{out_state}>: residual = {factor(diff_expanded[i,j])}")

# ============================================================
# PART 5: Now check the FULL six-vertex (b2 ≠ 0) to confirm it CAN satisfy YBE
# ============================================================

print("\n" + "=" * 70)
print("COMPARISON: FULL SIX-VERTEX (b2 ≠ 0)")
print("=" * 70)

# Use the Hecke parametrization: R(u) = u*T + I
# Weights: a1=a2=uq+1, b1=1, b2=u(q-1)+1, c1=uq, c2=u

def hecke_weights(param):
    """Return (a1, a2, b1, b2, c1, c2) for R(param) = param*T + I."""
    return (param*q + 1, param*q + 1, 1, param*(q-1) + 1, param*q, param)

a1p, a2p, b1p, b2p, c1p, c2p = hecke_weights(u-v)
a1r, a2r, b1r, b2r, c1r, c2r = hecke_weights(u)
a1s, a2s, b1s, b2s, c1s, c2s = hecke_weights(v)

R_6v_p = R_general(a1p, a2p, b1p, b2p, c1p, c2p)
R_6v_r = R_general(a1r, a2r, b1r, b2r, c1r, c2r)
R_6v_s = R_general(a1s, a2s, b1s, b2s, c1s, c2s)

R12_6v = compute_R12(R_6v_p)
R13_6v = compute_R13(R_6v_r)
R23_6v = compute_R23(R_6v_s)

LHS_6v = R12_6v * R13_6v * R23_6v
RHS_6v = R23_6v * R13_6v * R12_6v
diff_6v = (LHS_6v - RHS_6v).applyfunc(expand)

is_zero_6v = all(diff_6v[i,j] == 0 for i in range(8) for j in range(8))
print(f"\nFull six-vertex YBE (Hecke parametrization): {'VERIFIED' if is_zero_6v else 'FAILED'}")

if not is_zero_6v:
    print("Nonzero entries:")
    for i in range(8):
        for j in range(8):
            if diff_6v[i,j] != 0:
                print(f"  ({basis_labels[i]}, {basis_labels[j]}): {diff_6v[i,j]}")

# ============================================================
# PART 6: Set b2=0 in the Hecke parametrization specifically
# ============================================================

print("\n" + "=" * 70)
print("HECKE WITH b2 KILLED")
print("=" * 70)

# In the Hecke parametrization, b2(u) = u(q-1)+1.
# Setting this to 0 means u = 1/(1-q).
# But we want b2=0 for ALL u, which means we need a DIFFERENT parametrization.
#
# The five-vertex model is NOT a specialization of the Hecke R-matrix.
# It's a separate object. Let me use general polynomial weights.

# For the five-vertex model, the standard parametrization from the literature:
# a1(u) = 1, a2(u) = 1-qu, b1(u) = 1-u, b2(u) = 0, c1(u) = u(1-q), c2(u) = u
# (This is the "stochastic" five-vertex model, up to normalization.)
# Actually let me try: a1=1, a2=qu, b1=u, b2=0, c1=1-qu, c2=...
# There are various conventions. Let me just use fully symbolic and see what the
# residuals look like.

print("\nThe five-vertex model with fully symbolic weights:")
print("  R(u) with b2=0: a1u, a2u, b1u, c1u, c2u free")
print("  R(v) with b2=0: a1v, a2v, b1v, c1v, c2v free")
print("  R(w) with b2=0: a1w, a2w, b1w, c1w, c2w free")
print()
print("The nonzero residuals from Part 3 show exactly which configs fail.")
print("Now let me analyze each one in detail.")

# ============================================================
# PART 7: Detailed analysis of each failing equation
# ============================================================

print("\n" + "=" * 70)
print("DETAILED RESIDUAL ANALYSIS")
print("=" * 70)

for idx, (i, j) in enumerate(nonzero_positions):
    in_state = basis_labels[j]
    out_state = basis_labels[i]
    residual = diff_expanded[i, j]

    print(f"\n--- Failing equation {idx+1}: |{in_state}> -> |{out_state}> ---")
    print(f"Residual = {residual}")
    print(f"Factored = {factor(residual)}")

    # Also compute LHS and RHS separately
    lhs_val = expand(LHS[i, j])
    rhs_val = expand(RHS[i, j])
    print(f"LHS = {lhs_val}")
    print(f"RHS = {rhs_val}")

# ============================================================
# PART 8: Check if ANY choice of weights can make the five-vertex YBE hold
# ============================================================

print("\n" + "=" * 70)
print("CAN ANY WEIGHT CHOICE SAVE THE FIVE-VERTEX YBE?")
print("=" * 70)

# The residuals are polynomial expressions in the 15 weight variables
# (5 weights × 3 spectral parameters). For the YBE to hold, ALL residuals
# must vanish simultaneously. Let's check if this forces triviality.

print("\nNumber of independent residual equations:", nonzero_count)
print("\nChecking if residuals can simultaneously vanish for nontrivial weights...")

# Check each residual for what it implies
for idx, (i, j) in enumerate(nonzero_positions):
    residual = diff_expanded[i, j]
    print(f"\nResidual {idx+1}: {factor(residual)}")
    # Get the variables that appear
    vars_in_residual = residual.free_symbols
    print(f"  Variables: {vars_in_residual}")
