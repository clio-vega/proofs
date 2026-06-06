"""
Goertzen-Williamson Gram test for V_(3,1) of H_q(S_4) in the POLYTABLOID
(Specht / Murphy) basis instead of the seminormal basis.

The previous test (gw_gram_n4.py) used the seminormal Hoefsmit basis with the
Hecke-invariant diagonal Gram and got a "partial pass":
   <E_4 v, E_4 v>_G / <v, v>_G  =  4 q^3 (q^2+4q+1) / (1+q)^8.
The factor q^2+4q+1 is the conjectured P_{(3,1),1}(q); the prefactor
4 q^3 / (1+q)^8 is artefactual.

The Goertzen-Williamson framework lives on the polytabloid integral basis where
the Gram matrix has a clean integer expression at q=1, q-deforming to the
"standard form" on the Murphy basis. This script:

  (a) Builds the symmetric-group action of S_4 on V_(3,1) in the polytabloid
      basis, with rational integer entries.
  (b) Computes the standard Specht Gram (number of common columnar fillings).
  (c) Diagonalises (T_1, T_3) jointly to find the H-invariant vector v_poly at q=1.
  (d) Performs the "shrunk" GW test  <pi v, pi v>_G / <v,v>_G  where pi is the
      symmetriser (1+s_1)(1+s_2)(1+s_1)(1+s_3)(1+s_2)(1+s_1)/(2^6) at q=1.
  (e) Then lifts to Hecke at generic q by switching to the MURPHY basis (a
      q-deformation of the polytabloid basis): we identify a triangular
      change-of-basis from seminormal to "Murphy" via the Gram-Schmidt-of-JM
      decomposition, transport v and the Hecke-invariant Gram, and recompute.
"""

import sympy as sp
from sympy import Matrix, Rational, simplify, factor, expand, symbols, eye, zeros, sqrt, S
from itertools import permutations

q = symbols('q')

# -----------------------------------------------------------------------------
# Reuse the seminormal construction from gw_gram_n4.py
# -----------------------------------------------------------------------------
print("="*70)
print("Step 0: Rebuild seminormal H_q(S_4) action on V_(3,1)")
print("="*70)

T1 = sp.diag(q, q, -1)

A_val = -sp.Rational(1)/(q**2 + q + 1)
D_val = q**3/(q**2 + q + 1)
BC = sp.simplify(A_val*D_val + q)
T3 = sp.Matrix([
    [A_val,         1,     0],
    [BC,            D_val, 0],
    [0,             0,     q]
])

A2_val = -sp.Rational(1)/(q + 1)
D2_val = q**2/(q + 1)
BC2 = sp.simplify(A2_val*D2_val + q)
T2 = sp.Matrix([
    [q,      0,      0],
    [0,      A2_val, 1],
    [0,      BC2,    D2_val]
])

I3 = sp.eye(3)
# All Hecke relations checked in gw_gram_n4.py.

# Hecke-invariant diagonal Gram (from previous test)
ga = (q**4 + q**2)/(q**2 + q + 1)
gb = (q**3 + q**2 + q)/(q + 1)**2
gc = sp.Integer(1)
G_seminormal = sp.diag(ga, gb, gc)

# H-invariant vector in seminormal basis
M1 = T1 - q*I3; M3 = T3 - q*I3
combined = M1.col_join(M3)
ker = combined.nullspace()
assert len(ker) == 1
v_sn = ker[0]
denoms = [sp.fraction(sp.together(x))[1] for x in v_sn]
common = denoms[0]
for d in denoms[1:]:
    common = sp.lcm(common, d)
v_sn = sp.simplify(v_sn * common)
print(f"v in seminormal basis: {v_sn.T}")

Pi = (I3+T1)*(I3+T2)*(I3+T1)*(I3+T3)*(I3+T2)*(I3+T1)
Pi = sp.simplify(Pi)

# -----------------------------------------------------------------------------
# Step 1: Build POLYTABLOID basis of V_(3,1).
#
# SYTs of shape (3,1):
#   T_a = [[1,2,3],[4]]    column structure: col1=(1,4), col2=(2),  col3=(3)
#   T_b = [[1,2,4],[3]]    col1=(1,3), col2=(2), col3=(4)
#   T_c = [[1,3,4],[2]]    col1=(1,2), col2=(3), col3=(4)
#
# Polytabloid e_T = sum_{sigma in C(T)} sgn(sigma) sigma . {T}, where C(T) is
# the column stabilizer and {T} is the tabloid (row equivalence class of T).
#
# At q=1, S_4 acts on the span of polytabloids of standard tableaux as V_(3,1).
# The Specht inner product <e_T, e_T'> is given by the Young's inner product:
# orthonormalize tabloids -> compute polytabloid pairings.
#
# We build everything concretely using *tabloids* (a tabloid of shape (3,1) is
# a partition of {1,2,3,4} into a 3-set (row 1) and a 1-set (row 2)). There are
# 4 tabloids, indexed by which element is in row 2.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 1: Build polytabloid basis of V_(3,1) at q=1 (S_4)")
print("="*70)

# Tabloids: indexed by element in row 2. Tabloid t_k: row2 = {k}, row1 = {1,2,3,4}\{k}.
# So tabloid index k in {1,2,3,4} -> tabloid_index = k-1 in 0..3.
# S_4 acts on tabloids by permuting elements (move element k to position pi(k)).
def perm_action_on_tabloid(pi, k):
    """pi is a tuple of length 4 representing image of (1,2,3,4); k is the row-2 element.
       After permuting, the new row-2 element is pi(k)."""
    return pi[k-1]

# Polytabloids:
#   T_a = [[1,2,3],[4]] : col1 = {1,4}.  e_a = {T_a} - (14).{T_a}
#         {T_a} = tabloid with 4 in row 2  -> index 4
#         (14).{T_a} = tabloid with 1 in row 2 -> index 1
#         => e_a = t_4 - t_1
#   T_b = [[1,2,4],[3]] : col1 = {1,3}.  e_b = {T_b} - (13).{T_b}
#         {T_b} -> 3 in row 2 -> t_3
#         (13).{T_b} -> 1 in row 2 -> t_1
#         => e_b = t_3 - t_1
#   T_c = [[1,3,4],[2]] : col1 = {1,2}.  e_c = {T_c} - (12).{T_c}
#         {T_c} -> 2 in row 2 -> t_2
#         (12).{T_c} -> 1 in row 2 -> t_1
#         => e_c = t_2 - t_1
#
# So (e_a, e_b, e_c) in the basis (t_1, t_2, t_3, t_4) is:
#   e_a = (-1, 0, 0, 1)
#   e_b = (-1, 0, 1, 0)
#   e_c = (-1, 1, 0, 0)
#
# Polytabloid embedding into tabloid space (as columns):
P_embed = sp.Matrix([
    [-1, -1, -1],
    [ 0,  0,  1],
    [ 0,  1,  0],
    [ 1,  0,  0],
])  # 4x3, columns are e_a, e_b, e_c in tabloid basis (t_1..t_4)
print("Polytabloid embedding (columns e_a, e_b, e_c in tabloid basis t_1..t_4):")
sp.pprint(P_embed)

# Tabloid Young inner product: <t_i, t_j> = delta_{ij} (orthonormal).
# Then Specht Gram in polytabloid basis: G_poly = P_embed^T * I_4 * P_embed.
G_poly_q1 = P_embed.T * P_embed
print("\nPolytabloid Gram matrix at q=1 (G_poly = P^T P):")
sp.pprint(G_poly_q1)

# Build S_4 action on tabloids (4x4 permutation matrices for s_1, s_2, s_3).
def perm_matrix_on_tabloids(pi):
    """pi tuple of length 4. Returns 4x4 matrix M such that M[j,i]=1 iff pi sends
       tabloid t_i (row-2 element i+1) to tabloid t_{j+1} (row-2 element pi(i+1))."""
    M = sp.zeros(4, 4)
    for i in range(4):  # tabloid t_{i+1}: row-2 element is i+1
        new = perm_action_on_tabloid(pi, i+1)  # new row-2 element
        M[new-1, i] = 1
    return M

s1 = (2, 1, 3, 4)  # swap 1,2
s2 = (1, 3, 2, 4)  # swap 2,3
s3 = (1, 2, 4, 3)  # swap 3,4

S1_tab = perm_matrix_on_tabloids(s1)
S2_tab = perm_matrix_on_tabloids(s2)
S3_tab = perm_matrix_on_tabloids(s3)

# Project onto polytabloid subspace: find S_i_poly such that
#   s_i . P_embed = P_embed . S_i_poly   (action restricted to Specht subspace)
# Solve P_embed * S_i_poly = S_i_tab * P_embed
# S_i_poly = (P_embed^T P_embed)^{-1} P_embed^T S_i_tab P_embed  (least-squares)
G_poly_q1_inv = G_poly_q1.inv()

def restrict_to_poly(M_tab):
    return G_poly_q1_inv * P_embed.T * M_tab * P_embed

S1_poly = restrict_to_poly(S1_tab)
S2_poly = restrict_to_poly(S2_tab)
S3_poly = restrict_to_poly(S3_tab)

print("\ns_1 in polytabloid basis (q=1):")
sp.pprint(S1_poly)
print("\ns_2 in polytabloid basis (q=1):")
sp.pprint(S2_poly)
print("\ns_3 in polytabloid basis (q=1):")
sp.pprint(S3_poly)

# Verify these are involutions and satisfy braid relations.
for name, S in [("s1", S1_poly), ("s2", S2_poly), ("s3", S3_poly)]:
    assert sp.simplify(S*S - sp.eye(3)) == sp.zeros(3,3), f"{name}^2 != 1"
assert sp.simplify(S1_poly*S3_poly - S3_poly*S1_poly) == sp.zeros(3,3)
assert sp.simplify(S1_poly*S2_poly*S1_poly - S2_poly*S1_poly*S2_poly) == sp.zeros(3,3)
assert sp.simplify(S2_poly*S3_poly*S2_poly - S3_poly*S2_poly*S3_poly) == sp.zeros(3,3)
print("\nAll S_4 relations verified in polytabloid basis.")

# Sanity: at q=1, T_i = s_i (Hecke degenerates). The polytabloid basis matrices
# should be conjugate to the seminormal matrices at q=1.
T1_q1 = T1.subs(q, 1)
T2_q1 = T2.subs(q, 1)
T3_q1 = T3.subs(q, 1)
print("\nT_1 at q=1 (seminormal):")
sp.pprint(T1_q1)
print("T_1 at q=1 (polytabloid):")
sp.pprint(S1_poly)

# -----------------------------------------------------------------------------
# Step 2: Find change-of-basis P : seminormal -> polytabloid at q=1, by
# simultaneous similarity. We need invertible M with
#    S_i_poly = M^{-1} T_i_q1 M  for i=1,2,3.
# Solve for M by demanding it intertwines.
# Equivalently: find M with T_i_q1 M = M S_i_poly.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 2: Find change-of-basis seminormal -> polytabloid at q=1")
print("="*70)

m11, m12, m13, m21, m22, m23, m31, m32, m33 = symbols(
    'm11 m12 m13 m21 m22 m23 m31 m32 m33')
M_unk = sp.Matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

eqs = []
for T_sn, S_pl in [(T1_q1, S1_poly), (T2_q1, S2_poly), (T3_q1, S3_poly)]:
    R = T_sn * M_unk - M_unk * S_pl
    for i in range(3):
        for j in range(3):
            eqs.append(R[i,j])

sol = sp.solve(eqs, [m11,m12,m13,m21,m22,m23,m31,m32,m33], dict=True)
print(f"# free-parameter solutions: {len(sol)}")
print(f"Solution (with free params): {sol[0]}")

# Pick a specific normalization. The intertwiner is determined up to scalar
# (since V_(3,1) is irreducible). Substitute by setting one free parameter to 1.
free = [sym for sym in [m11,m12,m13,m21,m22,m23,m31,m32,m33] if sym not in sol[0]]
print(f"Free parameters: {free}")
subs_dict = {free[0]: 1}
for f in free[1:]:
    subs_dict[f] = 0
M_q1 = M_unk.subs(sol[0]).subs(subs_dict)
print("Change-of-basis M (seminormal -> polytabloid) at q=1:")
sp.pprint(M_q1)
print(f"det M = {sp.simplify(M_q1.det())}")

# Translate v_sn at q=1 into the polytabloid basis.
v_sn_q1 = v_sn.subs(q, 1)
print(f"\nv (seminormal) at q=1 = {v_sn_q1.T}")
v_poly_q1 = M_q1.inv() * v_sn_q1
v_poly_q1 = sp.simplify(v_poly_q1)
print(f"v (polytabloid) at q=1 = {v_poly_q1.T}")

# -----------------------------------------------------------------------------
# Step 3: At q=1, compute the polytabloid Gram and the GW ratio
#
#    <pi v, pi v>_{G_poly} / <v, v>_{G_poly}
#
# where pi = (1+s_1)(1+s_2)(1+s_1)(1+s_3)(1+s_2)(1+s_1) acting on the
# polytabloid-basis matrices, and G_poly = P_embed^T P_embed.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 3: q=1 Gram ratio in polytabloid basis")
print("="*70)
Pi_poly_q1 = ((sp.eye(3)+S1_poly)*(sp.eye(3)+S2_poly)*(sp.eye(3)+S1_poly)
              *(sp.eye(3)+S3_poly)*(sp.eye(3)+S2_poly)*(sp.eye(3)+S1_poly))
Pi_poly_q1 = sp.simplify(Pi_poly_q1)

E4_poly_q1 = Pi_poly_q1 / (1+1)**6  # (1+q)^6 at q=1 = 64
E4v_poly_q1 = sp.simplify(E4_poly_q1 * v_poly_q1)
print(f"E_4 v (polytabloid, q=1) = {E4v_poly_q1.T}")

def inner_poly(u, w, G):
    return sp.simplify((u.T * G * w)[0,0])

vv_q1 = inner_poly(v_poly_q1, v_poly_q1, G_poly_q1)
ee_q1 = inner_poly(E4v_poly_q1, E4v_poly_q1, G_poly_q1)
ratio_q1 = sp.simplify(ee_q1 / vv_q1)
print(f"\nAt q=1:")
print(f"  <v,v>_poly = {vv_q1}")
print(f"  <E_4 v, E_4 v>_poly = {ee_q1}")
print(f"  ratio = {ratio_q1}")
print(f"  P_{{(3,1),1}}(1) = 1 + 4 + 1 = 6")
print(f"  ratio / 6 = {sp.Rational(ratio_q1, 6)}")

# -----------------------------------------------------------------------------
# Step 4: Lift to generic q. We construct a "Murphy/polytabloid-like" basis at
# generic q by demanding the same change-of-basis structure. The cleanest q-lift:
# transport the polytabloid Gram via the change-of-basis M(q), where M(q) is
# the q-analogue intertwiner. We obtain M(q) by re-solving the intertwining
# equations T_i M = M S_i^{(q)}, where S_i^{(q)} is the q-deformed action that
# specialises to the polytabloid action at q=1.
#
# The natural q-deformation of the SYMMETRIC GROUP polytabloid action is via
# the MURPHY BASIS of the Specht module: the Murphy basis matrices are
# upper-triangular q-analogues, and the change-of-basis from seminormal to
# Murphy is a triangular matrix in the dominance order on tableaux.
#
# For shape (3,1) and SYTs ordered T_a > T_b > T_c by dominance (a is
# super-standard, c is sub-dominant), the Murphy basis is obtained from the
# seminormal basis by an upper-triangular change of basis.
#
# Rather than write out the explicit Murphy theory, we exploit a shortcut:
# we know the polytabloid basis at q=1, and the Hecke-invariant Gram
# G_seminormal at generic q. We TRANSPORT G_seminormal back to the polytabloid
# basis at q=1 (i.e. compute G_poly_seminormal = M(q=1)^T G_seminormal M(q=1)).
# This gives a *q-dependent Gram on the polytabloid basis*. Then the ratio is
# basis-invariant -- the polytabloid GW ratio equals the seminormal GW ratio,
# which we already computed!
#
# So the right test is: use the *integer* polytabloid Gram G_poly_q1 itself
# (NOT the Hecke-invariant one) and the seminormal v transported via M(q=1)
# (with q-dependence carried in v only). This is what GW propose.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 4: Generic q -- transport v(q) into polytabloid basis via M(q=1)")
print("="*70)

# Use M_q1 to transport v(q) (which is symbolic in q in the seminormal basis)
# into the polytabloid basis. CAUTION: this is NOT a Hecke-equivariant
# transport at generic q -- it's an ANSATZ that the polytabloid basis is the
# correct integral basis to measure the GW Gram in.
v_poly_q = sp.simplify(M_q1.inv() * v_sn)
print(f"v(q) in polytabloid basis (via M(q=1)):")
sp.pprint(v_poly_q)

# Symmetrize using polytabloid-basis matrices at q=1 (Specht/symmetric-group
# structure -- GW's framework is at q=1 in their main theorem)
E4_specht = Pi_poly_q1 / 64  # rank-1 idempotent of S_4 on V_(3,1)
E4v_specht = sp.simplify(E4_specht * v_poly_q)
print(f"\nE_4(q=1) v(q) in polytabloid basis:")
sp.pprint(E4v_specht)

vv_specht = inner_poly(v_poly_q, v_poly_q, G_poly_q1)
ee_specht = inner_poly(E4v_specht, E4v_specht, G_poly_q1)
ratio_specht = sp.factor(sp.simplify(ee_specht / vv_specht))
print(f"\n<E_4 v, E_4 v>_poly / <v,v>_poly (Specht Gram, q=1 symmetrizer):")
print(f"  = {ratio_specht}")

# -----------------------------------------------------------------------------
# Step 5: True q-deformed Gram. The right object: q-symmetrize using the
# *Hecke* operators (which depend on q), but measure with the polytabloid
# (Specht-integral) Gram. The Hecke operators are in the seminormal basis,
# so we conjugate them to the polytabloid basis at q=1 (the only place we
# have an explicit identification).
#
# The "honest" GW prescription:
#   G_poly is the q=1 integer Gram on Specht polytabloids.
#   The Hecke action acts via its own q-deformed matrices T_i^Hecke.
#   To compare, we use M(q=1) to identify the bases at q=1, and propagate.
#
# Let T_i_hecke_in_poly := M(q=1)^{-1} T_i M(q=1)  (seminormal->poly at q=1).
# This gives Hecke matrices in the poly basis (at q=1 they are S_i_poly).
# Build Pi_q in this basis, apply to v_poly_q, measure with G_poly_q1.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 5: Full Hecke symmetrizer in polytabloid basis, polytabloid Gram")
print("="*70)
Minv = M_q1.inv()
T1_poly = sp.simplify(Minv * T1 * M_q1)
T2_poly = sp.simplify(Minv * T2 * M_q1)
T3_poly = sp.simplify(Minv * T3 * M_q1)
print("T_1 in polytabloid basis (q-dependent):")
sp.pprint(T1_poly)

# Sanity: at q=1 these reduce to S_i_poly
for name, Ti, Si in [("T1", T1_poly, S1_poly), ("T2", T2_poly, S2_poly),
                     ("T3", T3_poly, S3_poly)]:
    diff = sp.simplify(Ti.subs(q, 1) - Si)
    assert diff == sp.zeros(3,3), f"{name}(q=1) != {name[0]}_poly"
print("Verified T_i(q=1) = S_i_poly.")

Pi_hecke_poly = ((sp.eye(3)+T1_poly)*(sp.eye(3)+T2_poly)*(sp.eye(3)+T1_poly)
                 *(sp.eye(3)+T3_poly)*(sp.eye(3)+T2_poly)*(sp.eye(3)+T1_poly))
Pi_hecke_poly = sp.simplify(Pi_hecke_poly)

E4_hecke_poly = sp.simplify(Pi_hecke_poly / (1+q)**6)
E4v_hecke_poly = sp.simplify(E4_hecke_poly * v_poly_q)

# v in poly basis using same M_q1
v_in_poly = sp.simplify(Minv * v_sn)
E4v_in_poly = sp.simplify(E4_hecke_poly * v_in_poly)

vv_full = inner_poly(v_in_poly, v_in_poly, G_poly_q1)
ee_full = inner_poly(E4v_in_poly, E4v_in_poly, G_poly_q1)
ratio_full = sp.factor(sp.simplify(ee_full / vv_full))
print(f"\nFull Hecke + polytabloid Gram ratio:")
print(f"  <v,v>_poly = {sp.factor(vv_full)}")
print(f"  <E_4 v, E_4 v>_poly = {sp.factor(ee_full)}")
print(f"  ratio = {ratio_full}")

# Note: this is just the seminormal ratio in disguise IF we use the
# Hecke-invariant Gram. We're using the POLYTABLOID gram instead, so this
# differs from the seminormal answer.

# -----------------------------------------------------------------------------
# Step 6: Compare with target P_{(3,1),1}(q) = q^2 + 4q + 1
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 6: Compare candidates with P_{(3,1),1}(q) = q^2 + 4q + 1")
print("="*70)
target = q**2 + 4*q + 1
print(f"Target: {target}")
print(f"  P(0) = {target.subs(q, 0)}")
print(f"  P(1) = {target.subs(q, 1)}")
print()

candidates = {
    "Specht-Gram + q=1 symmetrizer + v(q)": ratio_specht,
    "Specht-Gram + Hecke symmetrizer + v(q)": ratio_full,
}

for name, r in candidates.items():
    print(f"\n--- {name} ---")
    r = sp.simplify(r)
    rf = sp.factor(r)
    print(f"  raw = {rf}")
    print(f"  at q=0: {sp.simplify(r.subs(q, 0))}")
    print(f"  at q=1: {sp.simplify(r.subs(q, 1))}")
    # Try monomial rescalings
    matched = False
    for a in range(-12, 13):
        for b in range(-12, 13):
            cand = sp.simplify(r * (1+q)**a * q**b)
            if sp.simplify(cand - target) == 0:
                print(f"  *** EXACT after * (1+q)^{a} * q^{b}: matches q^2+4q+1 ***")
                matched = True
                break
        if matched:
            break
    if not matched:
        # try multiplicative scaling by a constant
        if sp.simplify(r) != 0:
            const = sp.simplify(target.subs(q, 0) / r.subs(q, 0)) if r.subs(q,0) != 0 else None
            print(f"  q=0 normalisation factor: {const}")
        # ratio with target
        ratio_target = sp.factor(sp.simplify(r / target))
        print(f"  r / target = {ratio_target}")

# -----------------------------------------------------------------------------
# Step 7: Try yet another normalisation: pick v with a different scale.
# v_sn was scaled to have polynomial entries. Try the natural Specht-norm
# normalisation: rescale v_poly so its first nonzero entry is 1.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("Step 7: Try Specht-normalised v (first nonzero entry = 1)")
print("="*70)
# Normalize v_in_poly so leading coefficient = 1
nz = None
for i in range(3):
    if sp.simplify(v_in_poly[i]) != 0:
        nz = v_in_poly[i]
        break
print(f"Leading nonzero entry of v_in_poly: {sp.factor(nz)}")
v_norm = sp.simplify(v_in_poly / nz)
print(f"Normalised v: {v_norm.T}")

E4v_norm = sp.simplify(E4_hecke_poly * v_norm)
vv_n = inner_poly(v_norm, v_norm, G_poly_q1)
ee_n = inner_poly(E4v_norm, E4v_norm, G_poly_q1)
ratio_n = sp.factor(sp.simplify(ee_n / vv_n))
print(f"\nNormalised ratio: {ratio_n}")
# basis-invariance: should equal ratio_full
print(f"Equal to ratio_full? {sp.simplify(ratio_n - ratio_full) == 0}")

# -----------------------------------------------------------------------------
# Step 8: Final summary -- record what works and what doesn't.
# -----------------------------------------------------------------------------
print()
print("="*70)
print("SUMMARY")
print("="*70)
print(f"Target P_{{(3,1),1}}(q) = {target}")
print(f"Polytabloid Gram matrix (q=1):")
sp.pprint(G_poly_q1)
print(f"Polytabloid Gram det = {G_poly_q1.det()}")
print()
print(f"v in polytabloid basis (q-dependent):")
sp.pprint(sp.simplify(v_in_poly))
print()
print(f"<E_4 v, E_4 v>_poly / <v,v>_poly = {ratio_full}")
print(f"At q=0: {sp.simplify(ratio_full.subs(q, 0))}")
print(f"At q=1: {sp.simplify(ratio_full.subs(q, 1))}")
