# Theorem B — the level-0 coefficient: VERIFIED, and PROVED for τ>0

**Date:** 2026-05-31 (prove session)
**Status:** PARTIAL, with the level-0 piece on solid footing. The level-0 coefficient of
Theorem B, `#{j:d_j=0}`, is reduced to a single-matrix invariant; the identity
`#{d_j=0}=#{s(T)=0}=K_{λ,β(E)}` is **verified n≤7** (project `fastp` survey + an independent
self-checking script) and **PROVED for all τ>0** (corollary of the established `Ω=0 on V^λ`).
The τ=0 case is reduced to a clean finite linear-algebra statement. Higher levels remain the
Puiseux frontier.

## Target (Theorem B)

`{d_j}(λ)={s(T):T∈SYT(λ)}` as multisets, `s(T)=Σ_{i∈Des(T)}w_i`, `w_i=2i−1` if `n−i` odd
else `0`; `{d_j}` = vanishing orders at `x=q²` of the eigenvalue branches of the Baxterised
w₀ monodromy `M(x)=∏_{a=1}^{N} Řᵢₐ(x)` along the staircase reduced word
`staircase_word(n)=[1, 2,1, 3,2,1, …]` (blocks `s_{k-1}…s_1`, k=2..n), with
`Řᵢ=P_q^{(i)}+r(x)P_{-1}^{(i)}`, `r(x)=(q²−x)/(q(x−1))`, `N=C(n,2)`, `P_q=(T_i+1)/(q+1)`.

Established earlier: min `=τ(τ+1)/2`; sum `Σd_j=ord det M=Σ_T s(T)`; conjugation
`{d_j}(λ')={C(n,2)−d_j}`; min-level multiplicity `K_{λ,ν}`, `ν=(2^{n−ℓ},1^{τ+1})`.

## Parity sets and level-0 combinatorics

`w_i≠0 ⟺ n−i odd ⟺ i≡n−1 (mod 2)`. Paid set `O` (one parity class), free `E={1..n−1}\O`;
`s(T)=0 ⟺ Des(T)⊆E`, and `#{T:Des(T)⊆E}=K_{λ,β(E)}` (Stanley EC2 Thm 7.19.7),
`β(E)=(2^{n/2})` (n even) / `(1,2^{(n−1)/2})` (n odd). Hand-check: `K_{(4,2),(2,2,2)}=3`.

## The level-0 coefficient — exact reduction

At fusion `r(q²)=0`, each factor collapses: `M(q²)=∏_a P_q^{(i_a)}` (staircase order). By
continuity of the roots of `det(tI−M(x))`,

    #{j:d_j=0} = #(nonzero eigenvalues of M(q²)) = n_λ − mult₀(charpoly M(q²)).      (L0)

Moreover `M(q²)` has eigenvalue 0 **semisimple** (alg = geo mult — verified by the project,
RESULTS.md), so equivalently `#{d_j=0} = rank M(q²) = rank ∏_a P_q^{(i_a)}`.

Theorem B at level 0 is then

    #{d_j=0} = K_{λ,β(E)}.        (★)

### (★) is VERIFIED for n≤7

- Project engine `fastp.py` (exact GF(p) ε-series, two primes, cross-checked) computes the full
  `{d_j}` for all 43 shapes n≤7 (`survey7.out`). Reading the multiplicity of `0`: e.g.
  `(4,2): {d_j}=[0,0,0,1,5,6,9,10,14]` → 3 zeros `= #{s=0}=K_{(4,2),(2,2,2)}=3`;
  `(4,1,1)`→1; `(3,2,1)`→2; `(2,2,2)`→1 — all equal `#{s(T)=0}`. Holds for every shape n≤7.
- Independent self-checking rep (`verify_B2.py`, Hecke relations verified): `#nonzero eig
  M(q²) == #{s(T)=0}` for ALL shapes n≤6, and — a bonus — for FIVE different reduced words
  for w₀ (`down,up,rowdec,coldec,down_rev`). So the level-0 count is, empirically,
  **word-independent** (the operator `∏P_q` is word-dependent, but its rank is not, n≤6).

### (★) is PROVED for τ>0

`τ>0 ⟹ |S|=0<τ(τ+1)/2 ⟹ Ω=(q+1)^N M(q²)=0` on `V^λ` (Strong Pillar 1 / the established
`Ω=0 on V^λ`). Then charpoly `=t^{n_λ}`, `#{d_j=0}=0`; and RHS`=0` (`τ>0 ⟹` no `T` with
`s(T)=0`). **Level-0 of Theorem B is closed for every τ>0 shape**, as a corollary of work
already done.

### (★) for τ=0 — the residual, a clean finite statement

For τ=0, `M(q²)≠0` (singular, since `det M(x)=r(x)^{Σ dim E_i^-}→0`). Residual:

    rank ∏_a P_q^{(i_a)} = K_{λ,(2^{n/2})}    (τ=0, n even; (1,2^{…}) odd).      (★₀)

No Puiseux content — linear algebra over `Q(q)` at `x=q²`. This is the cleanest remaining
closeable piece of B. `rank ∏P_q = K` reconfirmed for all τ=0 shapes n≤6 a third way
(`image_probe.py`, exact column-space rank).

**Two naive routes RULED OUT (image_probe.py).** (i) The `s(T)=0` seminormal coordinates do
NOT form a section of the image — the K×rank submatrix is singular for (2,1),(4,1),(3,2),
(3,1,1),(2,2,1). (ii) The image is NOT the coordinate subspace `span{e_T : s(T)=0}` — false for
almost every shape (true only for the trivial 1-dim cases and (2,2),(2,2,2)). So `∏P_q` lands
in a *mixed* subspace with no seminormal-coordinate labeling; any proof of (★₀) must use the
operator structure (rank identity / fusion), not a basis relabeling.

**Lead (still open):** `K_{λ,(2^{n/2})}` = dim of the `(2^{n/2})`-weight space of the GL-irrep
`V_λ`. Conjecture: the (word-independent) image of the staircase `∏P_q` is isomorphic to that
weight space — "fusion at q² pairs the strands" — but NOT via the seminormal coordinates (ruled
out above); the isomorphism must be intrinsic. If so, (★₀) proves by inevitability and templates
the higher levels. Untested.

### Bug caught (recorded so it isn't repeated)
A first standalone script (`verify_B.py`) applied `sp.nsimplify` to already-exact rational
matrix entries before taking the char poly; this corrupted entries and produced a SPURIOUS
mismatch at `(4,2)` (4 vs 3) and `(4,1,1)` (2 vs 1). Removing `nsimplify` (`verify_B2.py`)
restores agreement with the authoritative `fastp` survey. Lesson: never `nsimplify` exact input.

### GUARDRAIL: (★) is NOT the trivial symmetrizer x_{w0}
`∏(T_{i_a}+1)` along a reduced word is word-dependent (n=3:
`(T_1+1)(T_2+1)(T_1+1)−(T_2+1)(T_1+1)(T_2+1)=q(T_1−T_2)≠0`), and `x_{w0}=Σ_w T_w` acts as 0 on
every `V^λ≠(n)`. So `M(q²)` is genuinely the staircase projector product, not `x_{w0}` — even
though, interestingly, its *rank* turns out word-independent.

## Higher levels — the polygon ladder (the real frontier; NOT hand-waved)

`det(tI−M(x))=Σ_m(−1)^m a_m(x)t^{n_λ−m}`; when the Newton polygon in `(x−q²)` is the lower
hull with no leading cancellation, `ord a_m=Σ(m smallest d_j)`; (★) pins the first off-axis
vertex at `m=K_{λ,β(E)}`; `m=n_λ` is the proved det/sum. `a_m=Σ(m×m principal minors)`.
Two gaps, both needing the ε-series engine on examples:
1. no leading cancellation / Newton = eigenvalue-valuation multiset;
2. eigenvalue-valuation multiset vs Smith exponents (agree only on the determinant).
Per-eigenvalue↔SYT is REFUTED (fractional individual slopes; confluent Puiseux clusters —
RESULTS.md); only the multiset is integral.

## Next steps
1. Prove (★₀): `rank ∏P_q = K_{λ,(2^{n/2})}` for τ=0. Try the `(2^{n/2})`-weight-space /
   pairwise-fusion identification of the image; exploit the observed word-independence of rank.
2. ε-series Newton polygon (n≤5) for the polygon ladder (†) and the Smith-vs-Newton question.

## Verification
- `~/projects/scratch/2026-05-31-branch-multiset/` (project, authoritative): `fastp.py`,
  `survey7.out`, `RESULTS.md` — full `{d_j}` n≤7; min, sum, conjugation all hold.
- `~/projects/scratch/2026-05-31-thmB/`: `verify_B2.py` (independent level-0 + word test),
  `verify_B.py` (buggy `nsimplify` — kept as cautionary), `prove-journal.md`.
Convention-independent facts (`M(q²)=∏P_q`, (L0), τ>0 closure via `Ω=0`, n=3 word-dependence,
`K_{(4,2),(2,2,2)}=3`) and the n≤7 verification of (★) stand.
