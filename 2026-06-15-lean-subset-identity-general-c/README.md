# Lean snapshot — general-`c` subset & central-tip identities + `NL_3`

**Date:** 2026-06-15 · **Status:** compiles clean (Mathlib v4.30.0), zero `sorry`, zero warnings,
3 standard axioms (`propext`, `Classical.choice`, `Quot.sound`).

Snapshot of `tworow_d4_kernel/TworowD4Kernel/SubsetIdentityGeneralC.lean`. Formalises the two
pure-`Nat.choose` identities of the general-`c` Number Lemma `NL_c`
(`../2026-06-14-numberlemma-general-c.md`), plus the `c = 3` instance of the Number Lemma as a
corollary.

## Declarations

- **`choose_mul_choose_subset (F c j) (hj : 2*c ≤ j)`** — Move 1, the subset-of-a-subset identity:
  `C(F+2c−1, j) · C(j, 2c) = C(F+2c−1, 2c) · C(F−1, j−2c)`.
  Proof: `Nat.choose_mul` at `n = F+2c−1`, `s = 2c`, then `omega` for `(F+2c−1)−2c = F−1`.

- **`descFactorial_two_mul_eq (j c)`** — Move 2, the central-tip identity:
  `j^{(2c)} = (2c)! · C(j, 2c)` (`Nat.descFactorial j (2*c) = (2*c)! * j.choose (2*c)`).
  Proof: `Nat.descFactorial_eq_factorial_mul_choose` at `k = 2c`.

- **`number_lemma_C3 (F j) (2∣F) (2≤F) (6≤j) (j≤F+5)`** — the `c = 3` Number Lemma:
  `v₂(F) + 3 ≤ v₂(C(F+5, j)) + v₂(j(j−1)(j−2)(j−3)(j−4)(j−5))`.
  Assembled exactly as `number_lemma_C2`: central tip (`v₂(720) = 4`) + subset identity (`c = 3`)
  + the anchor `720·C(F+5,6) = F(F+1)…(F+5)` + the minimisation `v₂(F+2) + v₂(F+4) ≥ 3` (a
  case-split on `F mod 4`: one of `F+2`, `F+4` is divisible by `4`). The constant `+3 = β(3)` is
  taken from the proof note, not invented.

## Numerical pre-check (`precheck.py`)

Re-confirmed before formalising, 0 mismatches:
- subset identity: `c = 1..6`, even `F ≤ 200`, all valid `j`.
- central-tip identity: `c = 1..6`, `j = 0..199`.
- `c = 3` anchor identity `v₂(C(F+5,6)) + v₂(6!) = Σ_{i=0}^{5} v₂(F+i)`: even `F ≤ 78`.

## Build

```
cd tworow_d4_kernel && lake exe cache get && lake build
```
