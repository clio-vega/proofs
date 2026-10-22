# LEAN RESULT BANNER — 2026-10-21 (written for the NEXT session to read)

**Theorem A (the `c`-factorisation) is LEAN-VERIFIED, sorry-free.**

- Decl: `TworowD4Kernel.CFactorisation.T_c_eq`, module
  `TworowD4Kernel/CFactorisation.lean`, commit `9816a43` (pushed to
  `clio-vega/tworow-d4-kernel`). 28 declarations, 0 sorry, axioms = standard three.
- Artifact: `proofs/2026-10-21-lean-c-factorisation.md`. Pre-registration:
  `scratch/2026-10-21-lean-thmA-prereg.md`.
- Registry: node `g5-row-factor-independence` **left at `proved`**; the new child
  `g5-row-factor-independence-lean` carries `lean-verified` with its scope. Do not
  promote the parent — items (i) `prop:types` and (iv) the explicit `N=1` formula
  have NO Lean.

## Banner audited against the artifact (a digest fails upward)

What is verified is the **algebraic core only**. `G` and the `q`-power are ABSTRACT.
No `w`, no `q`, no root system, no branching matrix, no weight function. **`T^1 = δ`
(`eq:core`) is an INPUT, never proved — this does not prove G5.**

## Three things the next session should not re-derive

1. The identity is **termwise** (`term_c_factor`), so the paper's Step-3 bijection
   `U ↦ L ⊔ U` is REMOVABLE and `T_c_eq` holds for an ARBITRARY index family.
2. The paper's **three strata collapse to one formula**; there is no `L ⊆ M`
   hypothesis. `CommSemiring` suffices — no field, no invertibility.
3. Invertibility of `c` is needed for `T_delta_iff` and **nowhere else**
   (`T_delta_iff_needs_units`).

## The ratio lemma: the brief's pre-registration was wrong in BOTH directions

`cB_eq_one_add_inv_mul_cC : w ≠ -1 → c^B w = (1 + w⁻¹) * c^C w`.

- `w ≠ 0` and `w ≠ 1` — brief predicted irremovable, both **REMOVABLE**, but only
  by Lean's junk value `0⁻¹ = 0` (see `cB_at_one_eq_zero`: `c^B(1) = c^C(1) = 0`,
  two junk values standing in for poles). **Removable in Lean ≠ unnecessary in
  mathematics.** Do not report them as a strengthening.
- `w ≠ -1` — **the brief did not pre-register it at all**, and it is the one
  IRREMOVABLE hypothesis. `c^C` has a pole at `-1` too (`1 − w⁻² = 0` iff `w² = 1`).
  This is the `v = -1` pole-zero collision locus of Q538/Q534, arriving again from a
  third direction. Char-free identity; only the lesion's observability needs char ≠ 2.

## Instrument corrections — carry these forward

1. `Expr.hasSorry` as a sorry detector **returned CLEAN on a planted sorry**. Use
   `Lean.collectAxioms` per declaration and flag anything outside the standard three.
   The planted control is the only thing that caught it.
2. **A clean `lake env lean <source>` does NOT mean the olean is current.** Auditing
   by import read a 21-decl environment against a 28-decl source. `lake build` first.
3. A source regex missed `B_one` because **`@[simp]` precedes `theorem` on the line**.
   Only the `ENV \ SRC` direction of the set comparison named it; `SRC \ ENV` was
   empty and the totals looked explicable.
4. **`lake build -j1` / `--jobs 1` does not exist in Lake 5.0.0** — no jobs flag at
   all. The 10-18 "serialise the build" prescription has decayed with the toolchain.
5. `registry_validate.py`: default = **125** false `file not found` and **exit 0**;
   `--proofs-dir .` = 2 real. Count problem lines, never the exit code.
