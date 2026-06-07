# The d=4 valuation residual is not a function of the 4-quotient — and the two-row fast-path has no classical-OP home

**Date:** 2026-06-07 (prove session)
**Object:** `G_λ(i) = ⟨s_λ, ψ^m⟩`, `ψ = h₂ + i e₂`, `λ ⊢ 2m`; conjectured fiber law
`G_λ(i) = 0 ⟺ λ = (2,2)`. `π = 1+i`.

**Summary of what this session settled (two decisive results, both negative, both re-pointing the
programme; plus one law verified to `n ≤ 20`):**

1. **PRIMARY target settled — Outcome 3.** The valuation residual `residual(λ)` is **not** a function
   of the 4-quotient of `λ`. Clean counterexample at `n = 10`. The route "prove `residual = f(4-quotient)`
   ⟹ Gap A" is dead; `residual` is irreducibly **4-core**-dependent.
2. **Two-row fast-path closed.** `Q_b`/`I_b(m)` (the deflated imaginary part) is **not** any classical
   orthogonal-polynomial family (no match: Krawtchouk, Meixner, Hahn, dual Hahn, Chebyshev, MO#286705;
   non-square discriminants). The browse-suggested "import a known irreducibility theorem" route does
   not exist; (♦) is a self-contained Diophantine problem.
3. **Law verified to n ≤ 20 (all shapes):** the unique vanisher of `G_λ(i)` for `λ ⊢ 2m`, `m ≤ 10`, is
   `(2,2)`; the `core=(2,2)` finite-depth mechanism holds with 0 violations.

---

## Part A. The valuation residual is not quotient-determined

### A.0 A valuation lemma (used throughout)

> **Lemma 1.** For `0 ≠ z ∈ ℤ[i]`, `v_π(z) = v₂(N(z)) = v₂(Re(z)² + Im(z)²)`.

*Proof.* `π = 1+i`, `π̄ = 1−i = −i·π`, and `−i` is a unit, so `v_π(z̄) = v_π(z)`. Hence
`2 v_π(z) = v_π(z) + v_π(z̄) = v_π(z z̄) = v_π(N(z))`. Since `N(z) ∈ ℤ` and `N(π) = 2` with
`(2) = (π)²`, we have `v_π(N(z)) = 2 v₂(N(z))`. Divide by 2. ∎

This makes `v_π(G_λ)` a pure 2-adic computation on the integer `Re² + Im²`.

### A.1 The decomposition and the open sub-question

James–Kerber gives a bijection `λ ↔ (4core(λ), 4quot(λ))` where `4quot` is an ordered 4-tuple of
partitions (canonical once the abacus bead-count is `≡ 0 mod 4`; **Lemma 2** below). Since `λ ⊢ 2m` is
even, `|4core(λ)|` is even, so `G_{4core}` is defined. The proposed decomposition is, by definition,

```
v_π(G_λ) = v₂(f^λ) + [ v_π(G_{4core}) − v₂(f^{4core}) ] + residual(λ),
residual(λ) := v_π(G_λ) − v₂(f^λ) − v_π(G_{4core}) + v₂(f^{4core}).
```

The open sub-question (PROVE.md): **does `residual(λ)` depend only on `4quot(λ)`?** A positive answer
with a finiteness statement would prove Gap A (`G_λ=0 ⟺ λ=(2,2)`).

> **Lemma 2 (convention robustness).** The ordered 4-quotient computed from a bead-set of size
> `r ≡ 0 (mod 4)`, `r ≥ ℓ(λ)`, is independent of `r`. *(Verified for all `λ`, `n ≤ 14`; this is the
> standard James–Kerber invariance under adding 4 beads. It guarantees the grouping below is canonical.)*

### A.2 Theorem (Outcome 3)

> **Theorem 1.** `residual(λ)` is **not** a function of `4quot(λ)`.

*Proof (explicit counterexample, exact integer arithmetic).* Both partitions below have ordered
4-quotient `((), (), (1,), ())` — a single box on runner 2 — verified by the abacus (Lemma 2):

- `λ = (1,1,1,1,1,1)`: `|λ|=6`, `4core = (1,1)` (size 2), one box of quotient (`2 + 4·1 = 6`).
  `G_λ(i) = ⟨s_{1^6}, ψ^3⟩ = −i`, so `v_π = v₂(0²+1²) = 0`. `f^λ = 1`, `v₂ = 0`. `4core = (1,1)`:
  `G_{(1,1)} = ⟨s_{11}, ψ⟩ = i` (since `ψ = h₂ + i e₂` and `⟨s_{11}, e₂⟩ = 1`), `v_π = 0`; `f^{(1,1)} = 1`,
  `v₂ = 0`. Hence `residual = 0 − 0 − 0 + 0 = 0`.
- `λ = (4,4,2)`: `|λ|=10`, `4core = (4,1,1)` (size 6), one box of quotient (`6 + 4·1 = 10`).
  `G_λ(i) = −44 − 4i`, so `v_π = v₂(44² + 4²) = v₂(1952) = v₂(2⁵·61) = 5`. `f^{(4,4,2)} = 252 = 2²·63`,
  `v₂(252) = 2`. `4core = (4,1,1)`: `G_{(4,1,1)} = −2 + 6i`, `v_π = v₂(4+36) = v₂(40) = 3`;
  `f^{(4,1,1)} = 10`, `v₂(10) = 1`. Hence `residual = 5 − 2 − 3 + 1 = 1`.

Same 4-quotient `((),(),( 1),())`, different 4-cores `(1,1) ≠ (4,1,1)`, **`residual = 0 ≠ 1`**. ∎

A full scan of all 1536 shapes `λ ⊢ 2m`, `m ≤ 10`, finds **58 ordered-quotient groups** (and 10
unordered) that are non-constant in residual; `n=10` is the smallest. (Data:
`code/results/residual-vs-4quotient.csv`.)

### A.3 Characterising the failure (where the residual's information lives)

Because `(4core, 4quot)` is a *complete* invariant, "`residual` is a function of `(core,quotient)`" is
vacuously true and says nothing; the content is that the **quotient alone is insufficient** and the
**4-core** carries residual information. Two refinements pin this down:

- **Almost a function of core-valuations.** Grouping by `(v_π(G_{core}), v₂(f^{core}), 4quot)` leaves
  only `12/787` non-constant groups — all in the single large-core valuation class `(v_π=15, v₂=4)`.
  So the core contributes to `residual` *mostly* through its own pair `(v_π(G_{core}), v₂(f^{core}))`,
  but **not exactly**.
- **Runner-sensitivity.** At fixed single-box quotient, `residual` can depend on *which runner* the box
  occupies (core `(5,2,1,1,1)`: residual `2/4/4` for runners `2/3/0`). This is the fingerprint of a
  genuine **core × quotient interaction term**, which the additive bracket cannot reproduce.

**Conclusion for Gap A.** The valuation-level decomposition does **not** separate into
core-bracket + quotient-function. A proof of Gap A through this route must instead supply a direct
formula for `v_π(G_λ)` carrying an interaction/runner term, or — sufficient for vanishing alone —
prove that this interaction is **finite** whenever the quotient is nonempty.

---

## Part B. The two-row fast-path has no classical-orthogonal-polynomial home

Recall (prior cycle) the two-row law reduces to **(♦):** for `b ≥ 5`, the deflation
`Q_b(m) = I_b(m)/∏_{r=0}^{⌊(b−1)/2⌋}(m−r)` of `I_b(m) := Im G_{(2m−b,b)}(i)` has no rational root
(save the half-integer `(2b−1)/2` when `4|b`). The browse cycle proposed importing a Krawtchouk/Meixner
irreducibility theorem. This session refutes the premise.

> **Proposition 2 (no classical-OP identification).** For `b = 5..20`, neither `I_b(m)` nor the
> primitive `Q_b(m)` matches — up to rational scalar and affine reparametrisation `m ↦ αm+β` — any of:
> binary or `q`-Krawtchouk, Meixner, Hahn, dual Hahn, discrete Chebyshev/Gram, or the MO#286705
> "generalised binary Krawtchouk" family, over the scanned parameter ranges. Moreover `disc(Q_b)` is
> **non-square** for every `b = 5..14`, so `Gal(Q_b) ⊄ A_n`; there is no clean discriminant pattern.

*(Computational; full coefficient lists, factorisations and discriminants in
`scratch/2026-06-07-prove/`.)* The browse "binary-Krawtchouk-type" reading was a *generating-function
resemblance* (`(1−z)^i(1+z+z²)^{n−i}` shape), not an identity. Consequences:

- The "import known irreducibility" path to (♦) **does not exist**; (♦) is a self-contained Diophantine
  non-vanishing statement about a *new* integer polynomial family.
- Confirmed structure (b ≤ 16): `Q_b` irreducible for `4∤b`; for `4|b`, `Q_b = (2m−(2b−1))·(irreducible)`;
  leading coefficient `1` (resp. `2` when `4|b`). Constant terms carry a growing idiosyncratic prime, so
  the rational-root theorem alone does not preclude integer roots — irreducibility does.

**The law, separately, is rock-solid.** Exact integer scan: the only integer `(m,b)` with `m ≥ b`,
`b ≤ 40`, `m ≤ 3b²`, at which `I_b(m) = 0` is `(m,b) = (2,2)`. Combined with Part A's all-shape check,
`G_λ(i) = 0 ⟺ λ = (2,2)` holds for every `λ ⊢ 2m`, `n ≤ 20`.

---

## Verification
- `G_λ(i)` by iterated Pieri in `ℤ[i]` (`core.psi_power_schur`); cross-checked `G_{(2,2)} = 0` by hand.
- `v_π` via Lemma 1; `f^λ` via hook lengths; 4-core/4-quotient via 4-runner abacus (size self-check
  `|λ| = |core| + 4|quot|` passes for all `n ≤ 16`; quotient `r→r+4`-stable for `n ≤ 14`).
- Residual table: all 1536 shapes `λ ⊢ 2m`, `m ≤ 10`. Counterexample (Thm 1) re-verifiable by hand from
  the tabulated `G`, `f`, core values.
- Part B: identification sweep + discriminants, `b ≤ 20`.

## Gaps (precise)
- **Gap A (open).** Prove `G_λ(i)=0 ⟺ λ=(2,2)` for all `λ`. The residual route is dead (Thm 1); the live
  reformulation is a direct `v_π(G_λ)` formula with a finite core×quotient interaction term. Concrete
  first step: characterise the `12/787` exceptions to `residual = f(v_π G_core, v₂ f_core, quot)`.
- **(♦) (open).** `Q_b` has no rational root for `b ≥ 5` (save `(2b−1)/2`, `4|b`). No classical-OP
  certificate exists (Prop 2); needs a direct Galois/Newton-polygon (multi-prime, odd) or growth argument
  on a genuinely new family. 2-adic route is dead (`v₂(I_b)` unbounded, flat Newton polygon).
