"""
Chirality Obstruction — Functional equation analysis.

The 10 YBE residuals for the five-vertex model (b2=0) with general weights
give functional equations when we require the weights to come from functions
of a spectral parameter:
  R12 has parameter p (= u-v in difference form)
  R13 has parameter r (= u)
  R23 has parameter s (= v)

I derived by hand:
  c2 = c1 (constant ratio)
  a1 = c1 (constant ratio)
  a2 = B*c1 (constant ratio)
  b1(t) = alpha*t*c1(t) (linear ratio)
  With B=1 all residuals vanish.

This would mean R~(t) = [[1,0,0,0],[0,0,1,0],[0,1,alpha*t,0],[0,0,0,1]]
satisfies YBE — a NONTRIVIAL five-vertex model!

Must verify computationally before declaring the theorem false.
"""

from sympy import symbols, expand, Matrix, eye, zeros, Rational
from sympy import kronecker_product

u, v, alpha = symbols('u v alpha')

I2 = eye(2)

def R_five(t):
    """Five-vertex R-matrix: a1=a2=c1=c2=1, b1=alpha*t, b2=0."""
    return Matrix([
        [1,       0, 0,       0],
        [0,       0, 1,       0],
        [0,       1, alpha*t, 0],
        [0,       0, 0,       1]
    ])

def compute_R13(R):
    P = Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    P23 = kronecker_product(I2, P)
    R12 = kronecker_product(R, I2)
    return P23 * R12 * P23

# Standard YBE: R12(u-v) R13(u) R23(v) = R23(v) R13(u) R12(u-v)
R12 = kronecker_product(R_five(u - v), I2)
R13 = compute_R13(R_five(u))
R23 = kronecker_product(I2, R_five(v))

LHS = R12 * R13 * R23
RHS = R23 * R13 * R12
diff = (LHS - RHS).applyfunc(expand)

basis = ['000', '001', '010', '011', '100', '101', '110', '111']

print("Testing R~(t) with a1=a2=c1=c2=1, b1=alpha*t, b2=0")
print("YBE: R12(u-v) R13(u) R23(v) = R23(v) R13(u) R12(u-v)")
print()

nonzero = 0
for i in range(8):
    for j in range(8):
        if diff[i,j] != 0:
            nonzero += 1
            print(f"  ({basis[i]},{basis[j]}): {diff[i,j]}")

if nonzero == 0:
    print("ALL ZERO — YBE SATISFIED!")
    print(f"\nThis is a counterexample to the chirality obstruction theorem!")
    print(f"The five-vertex model with b1 = alpha*t, b2 = 0, a1=a2=c1=c2=1")
    print(f"satisfies the Yang-Baxter equation.")
else:
    print(f"\n{nonzero} nonzero entries — YBE FAILS.")
    print("The hand analysis had an error. Let me check more carefully.")

# ============================================================
# Also try the permutation-based five-vertex:
# R(t) = P + t * |10><10|  (where P swaps)
# ============================================================
print("\n" + "="*70)
print("VARIANT: R(t) = P + t * E_{10,10}")
print("="*70)

def R_perm_variant(t):
    """R(t) = permutation + t*(projection onto |10>)."""
    P = Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    E = zeros(4)
    E[2,2] = 1  # |10><10|
    return P + t * E

R12b = kronecker_product(R_perm_variant(u - v), I2)
R13b = compute_R13(R_perm_variant(u))
R23b = kronecker_product(I2, R_perm_variant(v))

LHSb = R12b * R13b * R23b
RHSb = R23b * R13b * R12b
diffb = (LHSb - RHSb).applyfunc(expand)

nonzero2 = 0
for i in range(8):
    for j in range(8):
        if diffb[i,j] != 0:
            nonzero2 += 1
            print(f"  ({basis[i]},{basis[j]}): {diffb[i,j]}")

if nonzero2 == 0:
    print("ALL ZERO — YBE SATISFIED!")
else:
    print(f"\n{nonzero2} nonzero entries — YBE FAILS.")

# ============================================================
# Let's try the most general approach: parameterize by spectral parameter
# and check what the residuals force
# ============================================================
print("\n" + "="*70)
print("GENERAL FIVE-VERTEX WITH SPECTRAL PARAMETER")
print("="*70)

# Use polynomial weights of degree ≤ 1 in the spectral parameter
# (most common for rational vertex models)
a, b, c, d, e, f, g, h, k, l = symbols('a b c d e f g h k l')

# a1(t) = a + b*t
# a2(t) = c + d*t
# b1(t) = e + f*t
# c1(t) = g + h*t
# c2(t) = k + l*t
# b2(t) = 0

def R_poly(t):
    a1t = a + b*t
    a2t = c + d*t
    b1t = e + f*t
    c1t = g + h*t
    c2t = k + l*t
    return Matrix([
        [a1t, 0,   0,   0  ],
        [0,   0,   c2t, 0  ],
        [0,   c1t, b1t, 0  ],
        [0,   0,   0,   a2t]
    ])

R12c = kronecker_product(R_poly(u - v), I2)
R13c = compute_R13(R_poly(u))
R23c = kronecker_product(I2, R_poly(v))

LHSc = R12c * R13c * R23c
RHSc = R23c * R13c * R12c

print("Computing YBE residual for degree-1 polynomial weights...")
diffc = (LHSc - RHSc).applyfunc(expand)

# Collect residual equations
print("Nonzero residuals:")
residuals = []
for i in range(8):
    for j in range(8):
        if diffc[i,j] != 0:
            residuals.append((i, j, diffc[i,j]))
            print(f"  ({basis[i]},{basis[j]}): {diffc[i,j]}")

print(f"\nTotal: {len(residuals)} nonzero entries")

# Each residual is a polynomial in u, v with coefficients that are polynomials
# in a,b,c,d,e,f,g,h,k,l. For YBE to hold for all u,v, each coefficient of
# each power of u,v must vanish.

print("\n--- Extracting coefficient equations ---")
from sympy import Poly, collect
from sympy import degree

# For each residual, expand as polynomial in u,v and extract coefficients
all_eqs = set()
for idx, (i, j, res) in enumerate(residuals):
    p = Poly(res, u, v)
    for monom, coeff in p.as_dict().items():
        if expand(coeff) != 0:
            all_eqs.add(expand(coeff))

print(f"Total independent coefficient equations: {len(all_eqs)}")
for eq in sorted(all_eqs, key=str):
    print(f"  {eq} = 0")
