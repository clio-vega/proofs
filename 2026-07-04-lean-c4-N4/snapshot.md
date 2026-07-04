# LEAN snapshot — c=4 interior Number Lemma N4 (`16 ∣ H`) machine-checked

**Date:** 2026-07-04 · **Session:** Lean · **Result:** sorry-free, standard axioms.

## Target declaration

- **Project:** `projects/lean/tworow_d4_kernel` (Mathlib v4.30.0, Lean 4).
- **File:** `TworowD4Kernel/ThreeRowC4InteriorN4.lean` (imported from the root `TworowD4Kernel.lean`).
- **Main lemma:**
  ```lean
  theorem TworowD4Kernel.N4 (a b j : ℤ) (h : a % 2 = b % 2) : (16 : ℤ) ∣ Hpoly a b j
  ```
  This is **Lemma N4** of `2026-06-19-c4-interior-number-lemma.md` §1 — the constant 2-adic floor
  `v₂ H ≥ 4` on the sextic heavy quotient `H(a,b,j)` of `Q₄(j) = (a−2)(b−3)H(j) + P₈(j)`, on the tie
  lattice `a ≡ b (mod 2)`.

## What builds sorry-free

Everything. The whole project builds green: `Build completed successfully (2093 jobs)`, **0 warnings,
0 sorries**. New declarations:

- `TworowD4Kernel.Hpoly {R} [CommRing R] (a b j : R) : R` — the sextic `H`, transcribed byte-for-byte
  from `code/threerow-c4/c4N4.py` (`Hval`). Numerically cross-checked in-kernel via `#eval` over ℤ:
  `H(10,8,8)=347088`, `H(7,7,0)=1306800`, `H(3,5,2)=105840`, `H(13,4,9)=1728` — all match Python.
- `TworowD4Kernel.N4_residue_key : ∀ x y j : ZMod 16, parZMod x = parZMod y → Hpoly x y j = 0` —
  the finite residue check, discharged by **`decide`** (the 2048 admissible residue triples of the
  paper proof). ~36 s kernel reduction; no `native_decide` (would add `Lean.ofReduceBool`).
- `TworowD4Kernel.N4` — the integer statement, obtained from the residue key by:
  1. `ZMod.intCast_zmod_eq_zero_iff_dvd`: `16 ∣ H ↔ (↑H : ZMod 16) = 0`;
  2. `push_cast; ring`: `(↑H : ZMod 16) = Hpoly ↑a ↑b ↑j`;
  3. transport parity through `map_intCast` + `ZMod.intCast_eq_intCast_iff`: `a % 2 = b % 2` becomes
     `parZMod ↑a = parZMod ↑b`, where `parZMod : ZMod 16 →+* ZMod 2`.

## Sorries

**None.**

## `#print axioms TworowD4Kernel.N4`

```
'TworowD4Kernel.N4' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.N4_residue_key' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three only.

## Notes / decisions

- **`decide` over `ZMod 16` is feasible** for this degree-6 polynomial: the full `∀`-over-4096-triples
  (2048 admissible under the parity filter) reduces in the kernel in ~36 s. The modulus of the paper
  proof (16) needed no reduction — raw `decide` sufficed. Recorded here in case the deeper Case B
  checks (`2¹² ∣ Π₈`, `2¹⁴ ∣ Π₁₀`) are ever formalised: those moduli (4096, 16384) are far too large
  for raw `decide` and would need the structural `(a+2)_{j-3}` factorisation, not brute residues.
- Import hygiene: `import Mathlib.Tactic` trips the `linter.style.header` warning (whole-tactic-folder)
  under the project's `weak.linter.mathlibStandardSet`. Replaced with the minimal
  `Mathlib.Tactic.{Ring, NormNum, Push}` + `Mathlib.Data.ZMod.Basic`.

## Ladder status

With N4 machine-checked, the **`(a,b,4)` interior is now Lean-witnessed at its keystone**, and the
boundary was already done (`threerow_c4_boundary`). The registry node `c4-number-lemma` is set to
`lean-verified` (`lean: "TworowD4Kernel.N4"`).
