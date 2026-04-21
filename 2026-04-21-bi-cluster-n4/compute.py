"""
Staircase n=4 principal-minor data for comparison with Bi's cluster algebra seed.

Outputs to /home/clio/projects/scratch/2026-04-21-staircase-n4-data.md:
  * D_4 ordering (6 permutations, sorted lexicographically)
  * Staircase word (1)(2,1)(3,2,1) and factor-by-factor decomposition
  * M^(1)[D_4], M^(2)[D_4], M^(3)[D_4] as 6x6 matrices over Z[q]
  * d_k(4) = det M^(k)[D_4] for k=1,2,3, factored over Z[q]
  * Eigenvalue decompositions of each M^(k)[D_4]
  * First-pair identity check d_3 * d_1 == q^6 * d_2^2

Re-uses compute_column_of_MR / build_MR_k / descent_paired from xmx_qdef_compute.py.
"""
import sys
sys.path.insert(0, '/home/clio/projects/scratch')

import sympy as sp
from sympy import symbols, factor, expand, Matrix, eye
from xmx_qdef_compute import (
    descent_paired,
    build_MR_k,
    staircase_word_partial,
    q,
)

OUT = "/home/clio/projects/scratch/2026-04-21-staircase-n4-data.md"
lines = []
def w(s=""): lines.append(s)

n = 4
D = descent_paired(n)
w(f"# Staircase n=4 principal-minor data (compiled 2026-04-21)")
w()
w(f"Hecke conventions: T_i^2 = (q-1)T_i + q.")
w()
w(f"Staircase element: Pi_q = (1+T_1)(1+T_2)(1+T_1)(1+T_3)(1+T_2)(1+T_1)")
w(f"  = B_1 B_2 B_3  where B_j = (1+T_j)(1+T_{{j-1}})...(1+T_1).")
w()
w(f"Truncations: Pi_q^(k) = B_1 B_2 ... B_k.  For n=4: k in {{1,2,3}}, Pi_q^(3) = Pi_q.")
w()
w(f"Matrix entry: M^(k)[u,v] = [T_u] (T_v * Pi_q^(k))  for u,v in D_n.")
w()
w(f"## D_4 ordering (lexicographic on one-line notation)")
w()
w(f"|D_4| = {len(D)}.  Entries:")
w()
w("| index | permutation |")
w("|-------|-------------|")
for i, perm in enumerate(D):
    w(f"| {i} | {perm} |")
w()

# Factor orderings / staircase words
w(f"## Staircase factor orderings")
w()
for k in (1, 2, 3):
    word = staircase_word_partial(k)
    factors_str = " · ".join(f"(1+T_{i})" for i in word)
    w(f"- `Pi_q^({k})` = {factors_str}  (word = {word})")
w()

# Build M^(k) and record
Ms = {}
dets = {}
for k in (1, 2, 3):
    M, _ = build_MR_k(n, k, q, D)
    Ms[k] = M
    dets[k] = expand(M.det())

# Print matrices
for k in (1, 2, 3):
    M = Ms[k]
    w(f"## M^({k})[D_4, D_4]  (6x6)")
    w()
    w("Rows indexed by u (output), columns by v (input).  M[u,v] = [T_u](T_v Pi_q^(k)).")
    w()
    w("| u \\\\ v | " + " | ".join(str(v) for v in D) + " |")
    w("|" + "---|" * (len(D) + 1))
    for i, u in enumerate(D):
        row = [factor(expand(M[i, j])) for j in range(len(D))]
        w(f"| {u} | " + " | ".join(str(x) for x in row) + " |")
    w()

# Determinants & factorizations
w(f"## d_k(4) = det M^(k)[D_4]")
w()
w("| k | det (expanded) | det (factored) |")
w("|---|----------------|----------------|")
for k in (1, 2, 3):
    det = dets[k]
    w(f"| {k} | `{det}` | `{factor(det)}` |")
w()

# First-pair identity
lhs = expand(dets[3] * dets[1])
rhs = expand(q**6 * dets[2]**2)
identity_holds = (expand(lhs - rhs) == 0)
w(f"## First-pair identity  d_3 * d_1 ?= q^6 * d_2^2")
w()
w(f"- LHS (expanded): `{factor(lhs)}`")
w(f"- RHS (expanded): `{factor(rhs)}`")
w(f"- **LHS == RHS?**  `{identity_holds}`")
w()

# Eigenvalues / char polys
w(f"## Eigenvalue decomposition of M^(k)[D_4]")
w()
for k in (1, 2, 3):
    M = Ms[k]
    w(f"### k = {k}")
    w()
    cp = expand(M.charpoly().as_expr())
    w(f"- Characteristic polynomial in t: `{factor(cp)}`")
    try:
        ev = M.eigenvals()
        ev_list = []
        for e, mult in ev.items():
            ev_list.append(f"{factor(sp.together(e))} (mult {mult})")
        w(f"- Eigenvalues:")
        for e in ev_list:
            w(f"    - {e}")
    except Exception as exc:
        w(f"- Eigenvalue computation failed: {exc}")
    w()

# Known closed-form cross-check
w(f"## Known closed forms (from proof 2026-04-18-first-q-shifted-pair.tex)")
w()
w(f"Expected:")
w(f"- d_1(4) = q^6")
w(f"- d_2(4) = q^9 (q+1)^3 (2q+1)")
w(f"- d_3(4) = q^18 (q+1)^6 (2q+1)^2  =  (d_2(4))^2")
w()
expected = {
    1: q**6,
    2: q**9 * (q + 1)**3 * (2*q + 1),
    3: q**18 * (q + 1)**6 * (2*q + 1)**2,
}
w("| k | matches closed form? |")
w("|---|----------------------|")
for k in (1, 2, 3):
    match = (expand(dets[k] - expected[k]) == 0)
    w(f"| {k} | {match} |")
w()
# Also check d_3 == d_2^2 (the n=4 squaring identity)
sq_ident = (expand(dets[3] - dets[2]**2) == 0)
w(f"- n=4 squaring identity  d_3(4) == d_2(4)^2 ?  **{sq_ident}**")
w(f"  (The corollary of locality: det M^(3) = (det M^(2))^2 at n=4.)")
w()

# write file
with open(OUT, "w") as fh:
    fh.write("\n".join(lines) + "\n")

# print concise summary to stdout
print(f"D_4 ordering: {D}")
print()
for k in (1, 2, 3):
    print(f"d_{k}(4) = {factor(dets[k])}")
print()
print(f"First-pair identity  d_3*d_1 == q^6 * d_2^2 :  {identity_holds}")
print(f"n=4 squaring  d_3 == d_2^2 :  {sq_ident}")
print()
print(f"Output written to {OUT}")
