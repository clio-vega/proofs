# PROVE 2026-10-19 c2 — the assignment's own grouping was wrong, and that was the first result

**Artifact:** `proofs/2026-10-19-c2-sa-divisors-full-q.tex` (13 pp., compiled, page-1 gate passed)
**On GitHub:** https://github.com/clio-vega/proofs/blob/main/2026-10-19-c2-sa-divisors-full-q.tex
**PDF:** https://github.com/clio-vega/proofs/blob/main/2026-10-19-c2-sa-divisors-full-q.pdf
**Scripts:** https://github.com/clio-vega/proofs/tree/main/code-1019c2-sa-converse
**Commit:** `6fad30c` (source + scripts), `7cd78ba` (page-1 block), pushed and verified as a set
difference against `origin/main`: 0 commits ahead, 0 files differing, 18/18 files present.

## The thing worth telling you

Rick assigned me, in writing twice, *"the boundary/A-type (S/A) divisors at full q"*. The brief
made me determine scope before attempting anything. **The scope determination is what produced the
best result: "S/A" is not one family, and the two halves have different answers.**

- **The A half carries no converse at all.** Theorem B of my 10-18 paper holds for *every* unit of
  the A-type GKM ring, so no renormalisation whatever produces an A-type pole. A statement
  "pole ⟺ combinatorial condition" there holds with the condition identically false.
- **The S half has one**, and I proved it.

The mechanism is structural and I want to record it because it says *why* the halves differ. On the
A-type hyperplane the two colliding fixed points carry the **same multiset of `w`-values**, so the
twist `Ψ_χ` has collision multiplier 1 for *every* `χ` — no free parameter is left to break. On the
boundary divisor the multisets differ by the single element `w_m` and the multiplier is `χ(ε)` —
**exactly one** free parameter. The assignment grouped a family with one degree of freedom together
with a family that has none.

## The theorem

Writing `γ^ν_M` for the residue at `{w_m = ε}` of the `(ν,M)` entry of `Φ^{-1}`, the boundary
collapse multiplier being *exactly 1* forces

    γ^ν_{R ∪ {m}} = − γ^ν_R        for every R ∌ m.

The residue functional is a **plain difference across the collision pairing**. That one identity
collapses a 2^N-term functional to

    Res_ν(Ψ_χ φ_λ φ_μ) = (1 − χ(ε)) · Σ_{R ∈ E_S} γ^ν_R · bar(Ψ_χ(R) φ_λ(R) φ_μ(R)),
    E_S(λ,μ,ν;m) = { R ⊆ [N]∖{m} : R ∪ {m} ≤ ν, R ≥ λ, R ≥ μ },

every term individually nonzero. Unconditionally: no pole if `χ(ε)=1` or `E_S = ∅`; a pole forces
`E_S ≠ ∅` and is simple. And **the converse, with the quantifier where it belongs**:

    E_S ≠ ∅  ⟺  some admissible χ with χ(ε) ≠ 1 makes ĉ^ν_{λμ} singular at {w_m = ε}.

The (⇒) half is a Vandermonde argument: the `N−1` linear forms `χ(w_r)` in the coefficients of
`χ` are independent, so the **squarefree monomials** `Π_{r∈R} χ(w_r)` over distinct `R` are linearly
independent, and the residue sum cannot vanish identically. *The quantifier is not a weakening* —
on this family the untwisted constants have no poles at all, so a converse must quantify over `χ`.
Rick's pair-divisor theorem needs no quantifier because there the untwisted constants already have
poles.

## Two gaps, named, and they are distinguishable by one experiment

1. **Hypothesis L^S**: `A ≥ B`, `m ∉ A` ⟹ `Ã_{A,B}|_{w_m=ε} ≠ 0`. The boundary analogue of Rick's
   Lemma L, whose degeneration `w_c = T^{-e_c}` cannot be run with `w_m` pinned at `ε`. Measured
   699/699 with a 1179/1179 plant. **This is the first thing to attempt next.**
2. **Non-cancellation at fixed χ** when `|E_S| ≥ 2`. Bypassed by the quantifier in the theorem.
   0 cancellations in 185 (N=3) and 1886 (N=4) opportunities, four to five twists each.

The pre-registered discriminator: **an FN that survives a change of `χ` refutes L^S; an FN that
moves with `χ` is cancellation.** One experiment separates them.

## By-products

- **Motegi's order is the tail-counting order**: `B ≤ A` ⟺ `#{c∈B : c≥k} ≤ #{c∈A : c≥k}` ∀k.
  Transitivity, `A∖{m} ≤ A`, and **complementation as an order anti-automorphism** all fall out of
  one line. 21844/21844, exhaustive, N ≤ 7. The anti-automorphism is why the support of Motegi's
  inverse formula agrees with the triangularity of `Φ^{-1}` — two mechanisms, same support.
- **Gate 0b closed at the definition.** `Δ^X_A := WX̃_{N,N,|A|}(w_A|w|A)` (L4257–4262) **is** my
  diagonal entry `Ã_{A,A}`. So Motegi's all-cardinality biorthogonality hands me `Φ^{-1}` in my own
  notation. I verified the parent theorem **first-hand against my own symmetrisation engine**,
  **336/336**, both types, N ≤ 3 — a different mechanism from his branching+RTV proof, so
  corroboration rather than a second reading of one engine. **This discharges Gap 1 of Rick's
  D∞ at N ≤ 3 from my side**, and it is worth telling him: that one unverified dependency sits
  under *three* results (his D∞, his Theorem W det lemma, my inverse lemma).

## Pre-registration outcomes, including the third arm

- **P1 → SPLIT**: arm 3 for A (no converse exists), arm 1 for S. **The presupposition failed, not
  the claim.** The brief's arm-3 detector fired, but on the S half and in type B only: `{w_m = −1}`
  is not a root hyperplane of `B_N`.
- **P2 → arm 3.** `Φ̂` *is* a different object (the type-A, unit-diagonal factor). The gate's
  purpose was served by the parent. The brief's instruction — *specialise the defining equation,
  not the statement* — is what separated them, and the `E^B`-error shape did not recur.
- **P3 → arm 3.** "Canonically" has no content: the boundary colength is one count by one index
  set, so there is no second index set for the first to be canonically isomorphic to. What
  biorthogonality actually bought was not the colength but the **non**-vanishing of `γ` on its
  support, which no colength argument gives.

## Two instrument faults, both of which produced confident legible wrong numbers

1. **A non-generic sample read as a phenomenon.** `ord_S det Φ = 8` against a predicted 4. The
   random point had `w_3 = 7/3` with `q = 3/7`, i.e. `w_3 = 1/q`, so `(1−q w_m w_3)` and
   `(1−q w_3/w_m)` vanish at `w_m = +1` — and the extra order appeared for exactly the four `I`
   that diagnosis predicts. ***The hypothesis "at generic q no other root factor vanishes on S" is
   a condition on the SAMPLE POINT, not only on q.*** Caught only because I had written the
   predicted number down before the run. Fixed with a genericity gate + assert.
2. **A simplifier that guesses.** `sympy.nsimplify` on a rational returned algebraic nonsense
   (`2**(41/58)`) and manufactured **20 false mismatches** in the biorthogonality check.
   ***`nsimplify` is a guesser, not a simplifier.*** `Rational(cancel(·))` gives 336/336.

One read as a catastrophe, one as a phenomenon. **The sign of the error said nothing about its
cause** — same shape as the four self-matching detectors of 10-19 c1.

## Trust boundary

`tracecheck`: 24 events, schema OK, **0 violations**. Last session 3 of 3 of my verify events were
flagged circular because I read `inputs` as *what this is about*; this session I emitted the
evidence as its own event and pointed `inputs` at dependencies only. The one verification that *is*
conditional — the non-cancellation census, which uses the residue formula — is **flagged in the
event itself** as `cross-verified`, not `independent`.

`trustcheck --files-dir .`: **`OK: ... is valid`**, 0 file-not-found. The 9 warnings are the
pre-existing `arXiv: names unknown peer` format noise (7 of them predate today).
