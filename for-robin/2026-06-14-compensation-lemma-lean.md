# Lean: Compensation Lemma C (`c=1`) + Number Lemma (`c=2`) — both done, target exceeded

**Clio, lean session 2026-06-14.**

## TL;DR

The lemma `LEAN.md` asked for — the three-row **Compensation Lemma C** — is formalised, and so is
the fast follow-on, the **`c=2` Number Lemma**. No `sorry`, no warnings, only the three standard
axioms. Two new files in `tworow_d4_kernel/TworowD4Kernel/`:
`CompensationLemma.lean` and `NumberLemmaC2.lean`. Snapshot in the proofs repo at
`2026-06-14-lean-compensation-and-number-lemma/`.

The headline theorems (namespace `TworowD4Kernel`):

```lean
theorem compensation_lemma_C (a b j : ℕ) (hj1 : 1 ≤ j) (hja : j ≤ a + 1)
    (hD : ((b * (a + 1) : ℕ) : ℤ) - ((j * (j - 1) : ℕ) : ℤ) ≠ 0) :
    padicValNat 2 (a + 1) ≤
      padicValNat 2 ((a + 2).choose j)
        + padicValInt 2 (((b * (a + 1) : ℕ) : ℤ) - ((j * (j - 1) : ℕ) : ℤ))

theorem number_lemma_C2 (F j : ℕ) (hF2 : 2 ∣ F) (hF : 2 ≤ F) (hj4 : 4 ≤ j) (hjF : j ≤ F + 3) :
    padicValNat 2 F + 1 ≤
      padicValNat 2 ((F + 3).choose j) + padicValNat 2 (j * (j - 1) * (j - 2) * (j - 3))
```

Plus the reusable helper you flagged for `c ≥ 2`:

```lean
theorem vz_choose_ge {n k : ℕ} (hk : 1 ≤ k) (hkn : k ≤ n) :
    vz n - vz k ≤ vz (n.choose k)          -- v₂(C(n,k)) ≥ v₂(n) − v₂(k)
```

(`vz n = (padicValNat 2 n : ℤ)`, reused from `D0ClosedForms.lean`.)

## Two corrections to the informal statements (please sanity-check me)

Both lemmas are written in the papers "for all `j`", but each proof silently uses a restriction
that is **genuinely necessary** — without it the literal statement is false under the standard
`padicVal … 0 = 0` convention. I formalised the *true* statement with the restriction explicit,
after a numerical sweep to pin the exact hypothesis. Neither correction touches the actual
application (the restriction always holds there), so the closures stand — but the paper statements
of Lemma C and the Number Lemma should be amended to match.

1. **Lemma C needs `D := b(a+1) − j(j−1) ≠ 0`.** The §4 proof's "strict-min rule" assumes it. At
   `D = 0` the bare inequality fails: smallest witness `(a,b,j) = (3,3,4)` gives
   `v₂(C(5,4)) + v₂(0) = 0 < 2 = v₂(4)`. But in the application range `1 ≤ j ≤ b−1`, `a ≥ b ≥ 1`,
   we get `b(a+1) ≥ b(b+1) > (b−1)(b−2) ≥ j(j−1)`, so `D > 0` always — the corner is unreachable
   there. **Sweep** (`a ≤ 400`, all `b`, `1 ≤ j ≤ a+1`): 630 violations of the bare form, *every
   one* with `D = 0`; with `D ≠ 0`, zero violations.

2. **Number Lemma needs `j ≤ F + 3`.** The subset identity
   `C(F+3,j)·C(j,4) = C(F+3,4)·C(F−1,j−4)` is only a *valuation* identity while all four binomials
   are nonzero. The unbounded "all `j ≥ 4`" form fails, e.g. `F=8, j=12`:
   `v₂(8)+1 = 4 > 0 + 3`. In the application `j ≤ b ≤ F+3` (with `F ∈ {a, b−1}`, so
   `F+3 ∈ {a+3, b+2} ≥ b`), so the used range is exactly `j ≤ F+3`. **Sweep** (even `F ≤ 1000`,
   `4 ≤ j ≤ F+3`): 250 500 cases, zero violations.

These are the same kind of convention-corner I'd want flagged in a referee report — the math is
right, the quantifier is one notch too generous.

## What carried the proofs

**Lemma C.** Over `ℤ` (the quadratic `b(a+1)−j(j−1)` is signed, hence `padicValInt`). The key
Mathlib import is `padicValInt_dvd_iff` (`(p:ℤ)^n ∣ a ↔ a = 0 ∨ n ≤ padicValInt p a`), which turns
divisibility of the difference into the valuation bound — and the `a = 0` branch is exactly where
the `D ≠ 0` hypothesis is spent. I unified the paper's case split into a single split on
`2^A ∣ Q` (`A = v₂(a+1)`, `Q = j(j−1)`): if it holds, `2^A ∣ D` directly gives
`A ≤ padicValInt 2 D`; otherwise `e := v₂(Q) < A`, and the absorption
`(a+2)·C(a+1,j−1) = C(a+2,j)·j` (via `vz_choose_ge`) gives `v₂(C(a+2,j)) ≥ A − e` while
`2^e ∣ D` gives `padicValInt 2 D ≥ e`. Sum `≥ A`.

**Number Lemma.** A clean `ℕ` inequality (no negatives in the conclusion). Two ingredients:
`j(j−1)(j−2)(j−3) = 24·C(j,4)` (from `Nat.descFactorial_eq_factorial_mul_choose`; I unfold
`descFactorial 4` by hand into the product, `descFactorial_four`) with `v₂(24)=3`; and the subset
identity `Nat.choose_mul` plus the exact value `3 + v₂(C(F+3,4)) = v₂(F)+v₂(F+2)` (only `F`, `F+2`
even among four consecutive integers), giving `v₂(C(F+3,4)) ≥ v₂(F) − 2`.

Every step is a named identity or `omega`/`linarith`/`positivity`. No `change`, no non-terminal
`simp`, no `native_decide` (would have added an axiom).

## Scope

I deliberately did **not** touch the symmetric-function assembly (the `M_j` determinant, `G_λ`,
`g(j)>0`, the Prop-2 collapse) — that needs machinery Mathlib doesn't have. These two files are the
self-contained 2-adic kernels: the parts a skeptic most wants machine-checked, and the parts that
generalise (the `vz_choose_ge` helper and the subset-identity pattern are wanted again for the open
`c=3` two-generator Number Lemma). The `c=3` analogue (`v₂C(b+3,j)+v₂Q₃(j)≥v₂(b+3)+1` and
Compensation Lemma B) is the natural next lean target once its prove session lands.

— Clio
