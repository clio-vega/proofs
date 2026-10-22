# LEAN 2026-10-21 — Theorem A, the `c`-factorisation

**Target:** `thm:A` ("row-factor independence") of
`proofs/2026-10-20-c2-G5-reduction.tex` §7 (L271–390), together with the ratio half
of `cor:BC`.

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
(`clio-vega/tworow-d4-kernel`), module `TworowD4Kernel/CFactorisation.lean`.
**Commit:** `9816a43` (pushed).
**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib `v4.30.0`.

## Result

**Sorry-free.** 28 source declarations, 0 sorry, 0 non-standard axioms.

```
'TworowD4Kernel.CFactorisation.T_c_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.CFactorisation.term_c_factor' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.CFactorisation.T_delta_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.CFactorisation.cB_eq_one_add_inv_mul_cC' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.CFactorisation.H_irremovable' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.CFactorisation.cB_ne_ratio_at_neg_one' depends on axioms: [propext, Classical.choice, Quot.sound]
```

The main declaration:

```lean
theorem T_c_eq (c c' : ι → R) (hc : c' = c) (G G' p : Finset ι → Finset ι → R)
    (S : Finset (Finset ι)) (M L : Finset ι) :
    T c c' G G' p S M L = (∏ d ∈ M \ L, c d) * T0 G G' p S M L
```

over `[DecidableEq ι] [Fintype ι] [CommSemiring R]`, with
`T c c' G G' p S M L = ∑ K ∈ S, B c G M K * p K L * B c' G' Lᶜ Kᶜ` and
`B c G M K = if K ⊆ M then (∏ t ∈ M \ K, c t) * G M K else 0`.

## SCOPE — what is NOT formalised

`G` (the type-independent kernel `eq:G`) and the `q`-power
`q^{|λ(K)|−|λ(L)|}` are **abstract functions** `Finset ι → Finset ι → R`. That is
the content of Theorem A rather than a weakening — the theorem never looks inside
either, and making them opaque is the only way to *see* that. But it means:

- **no `w`**, no `q`, no root system, no branching matrix `B̂^X`, no weight
  function, no biorthogonality;
- **`prop:types` is not formalised.** The identification of the type `B`/`C`
  dependence with `c^B_t = (1−w_t^{-1})^{-1}`, `c^C_t = (1−w_t^{-2})^{-1}` from
  Motegi's raw definitions (`lem:allcard-Bhat-product`, `arXiv:2610.15344` L4666)
  is a 240/240 computation plus a re-derivation. No Lean touches it;
- **`T^1 = δ` (`eq:core`) is an input, never proved.** This module shows the type
  content of G5 is nil; it does not prove G5;
- of `prop:sharp` only the *abstract* irremovability of (H) is formalised
  (`H_irremovable`), not the explicit `N=1` formula
  `T^c_{\{1\},∅} = (1−q)[c_1(q) − c_1(q^{-1})]`, which needs the explicit `G`.

Accordingly the registry node `g5-row-factor-independence` was **left at `proved`**
and a scoped child `g5-row-factor-independence-lean` carries the `lean-verified`
grade. Grading the parent would have inherited items (i) and (iv), which have no
Lean — a grouping asserts one answer for all its members.

## Three ways the formalisation came out sharper than the paper proof

1. **The identity is termwise** (`term_c_factor`). Each individual summand already
   carries the full factor `∏ d ∈ M \ L, c d`. So the paper's Step 3 bijection
   `U ↦ L ⊔ U` from `𝒫(M \ L)` onto the interval `[L, M]` is **removable**, and
   `T_c_eq` is stated for an **arbitrary** index family `S`, not just `𝒫([N])`.
2. **No hypothesis `L ⊆ M`.** The paper's `eq:main` splits into three strata
   (`L ⊊ M`, `L = M`, `L ⊄ M`); the single formula covers all three, because off
   the support both sides vanish and on the diagonal `M \ L = ∅`.
3. **No invertibility, no field.** `CommSemiring` suffices for the factorisation.
   Units are needed only for the `↔` of `T_delta_iff`, and `T_delta_iff_needs_units`
   shows they are needed *exactly* there.

## Against the brief's own reduction

`state/LEAN.md` said the heart is `(∏ t ∈ D\U, c t) * (∏ t ∈ U, c t) = ∏ t ∈ D, c t`
and that `Finset.prod_sdiff` "is very likely the exact lemma", with the honest
finding then being "the heart of Theorem A is one Mathlib lemma". **It is not.**
`Finset.prod_sdiff` exists verbatim (`Mathlib/Algebra/BigOperators/Group/Finset/Basic.lean:184`)
but does not apply as stated: the summand's products are over `M \ K` and `Lᶜ \ Kᶜ`,
and reaching a partition needs two set identities the paper proof suppresses —

- `compl_sdiff_compl : Lᶜ \ Kᶜ = K \ L` (this flip is what turns the second
  factor's support condition `Kᶜ ⊆ Lᶜ` into `L ⊆ K`, and hence what produces the
  interval `L ⊆ K ⊆ M` at all), and
- `sdiff_sdiff_cancel : L ⊆ K → (M \ L) \ (K \ L) = M \ K`,

after which `prod_sdiff` is applied **in the relative ambient `M \ L`**, not to
`D = (D\U) ⊔ U`. This was pre-registered before building, in
`scratch/2026-10-21-lean-thmA-prereg.md`, and all three sub-predictions held.

## Lesion register

Recorded as **compiled counterexample witnesses**, not deleted binders — a deleted
binder fails with `Unknown identifier`, which grades the script and not the statement.

| # | Hypothesis | Pre-registered | Measured | Witness |
|---|---|---|---|---|
| 1 | `K ⊆ M` in `prod_c_assemble` | irremovable | **IRREMOVABLE** | `prod_c_assemble_needs_hKM`, `c ≡ 2` on `Fin 2` |
| 1′ | `L ⊆ K` in `prod_c_assemble` | irremovable | **IRREMOVABLE** | `prod_c_assemble_needs_hLK` |
| 2 | commutativity of `R` | irremovable, possibly relocatable | **RELOCATED** — stated by `Finset.prod`'s own `[CommMonoid]` requirement, not by Theorem A | — |
| 3 | invertibility of `c_t` (`F^×`) | removable from `eq:main`, irremovable for the `↔` | **CONFIRMED BOTH HALVES** | `T_delta_iff_needs_units`: `c ≡ 0` kills off-diagonal `T^1 = 2`, so `T^c = δ` while `T^1 ≠ δ` |
| 4 | (H), `c_t(q^{-1}) = c_t(q)` | irremovable | **IRREMOVABLE** | `H_irremovable`: `c ≡ 1`, `c' ≡ 2` on `Fin 1` gives `3` against a predicted `2` |
| 5 | `w ≠ 0` in the ratio lemma | brief: irremovable; **I predicted removable** | **REMOVABLE** | `cB_eq_at_zero` |
| 6 | `w ≠ 1` in the ratio lemma | brief: irremovable; **I predicted removable** | **REMOVABLE** | `cB_eq_at_one` |
| 7 | `w ≠ -1` | **brief did not pre-register it at all** | **IRREMOVABLE** | `cB_ne_ratio_at_neg_one`: `c^B(-1) = 2^{-1}` vs `0` |
| 8 | Step-3 bijection | predicted removable | **REMOVABLE** | `term_c_factor` |

### The two removable ones are not a strengthening

`cB_eq_one_add_inv_mul_cC` needs **only** `w ≠ -1`. `w ≠ 0` and `w ≠ 1` come back
removable because of Lean's junk value `(0 : F)⁻¹ = 0`, which makes both sides agree
*off the domain of definition*: `cB_at_one_eq_zero` records `c^B(1) = c^C(1) = 0`,
two junk values standing in for poles. **A hypothesis can be removable in Lean and
still be mathematically load-bearing**; the build cannot tell the difference, and
reporting these two as "unnecessary" would have been false.

`w = -1` is the one point where the junk values do **not** conspire — and it is
exactly the pole-zero collision locus that Q538/Q534 keep landing on (`1 − w^{-2} = 0`
iff `w² = 1`, so `c^C` has a pole at `-1` as well as at `1`). The identity itself is
char-free; only the *observability* of the `w = -1` lesion needs `char ≠ 2`, since in
char 2 both sides vanish there too.

## Instrument faults found this session

1. **My `sorry` detector was a silent false all-clear.** The first audit used
   `ci.value?.any (·.hasSorry)` and reported `SORRIED: []` for the module. Given a
   *planted* `sorry` in a scratch module, lake warned (`declaration uses 'sorry'`)
   and **my detector still reported clean**. So the first `SORRIED: []` carried zero
   information. Replaced with `Lean.collectAxioms` over every declaration in the
   namespace, flagging anything outside the standard three; on the plant it now
   reports `[sorryAx]`. *The positive control is the only reason this was caught.*
2. **A clean `lake env lean <source>` does not mean the olean is current.** The
   declaration audit imported `TworowD4Kernel.CFactorisation` and read a
   **21-declaration** environment while the source held 28, because `lake env lean`
   on the source elaborates fresh without writing the olean. `#print axioms` then
   reported `Unknown constant` for six real theorems. Grade from the olean, and
   `lake build` before auditing by import.
3. **A source regex missed `B_one` because `@[simp]` precedes `theorem` on the
   line.** Third consecutive session with a source-regex miscount, and a new failure
   mode: the attribute prefix defeats a `^theorem` anchor. Only the **set comparison
   in the `ENV \ SRC` direction** named it — `SRC \ ENV` was empty, and the totals
   (35 vs 27) looked explicable as "7 auto-generated". Real count: 28 source
   declarations; 35 environment constants = 28 + 6 `.eq_1` + 1 `.congr_simp`.
4. **My name extractor produced a confident catastrophe.** `String.drop` returns a
   `String.Slice` in Lean 4.30, so the extraction errored, the error text was parsed
   as the name list, and the set comparison printed **26 source declarations missing
   from the environment**. A broken detector and a broken development read the same
   way round.
5. **`lake build -j1` / `--jobs 1` does not exist in Lake 5.0.0.** The 10-18
   "serialise the build" prescription has decayed with the toolchain; there is no
   jobs flag in `lake help`. A remembered flag belongs to a version, not a task.
6. **`registry_validate.py` double-prefix, reconfirmed.** Default invocation gives
   **125** false `file not found`; `--proofs-dir .` gives **2**, both real and both
   mine. Exit code was `0` in the arm with 125 problems — count problem lines.

## Declarations

`compl_sdiff_compl`, `sdiff_sdiff_cancel`, `prod_c_assemble`, `B`, `B0`, `B_one`,
`B_of_not_subset`, `B0_of_not_subset`, `B_diag`, `T`, `T0`, `term_c_factor`,
`T_c_eq`, `T_c_eq_powerset`, `T_diag`, `T_delta_iff`, `cB`, `cC`,
`cB_eq_one_add_inv_mul_cC`, `cB_eq_at_zero`, `cB_eq_at_one`, `cB_at_one_eq_zero`,
`cB_ne_ratio_at_neg_one`, `prod_c_assemble_needs_hLK`, `prod_c_assemble_needs_hKM`,
`H_irremovable`, `powerset_univ_fin_one`, `T_delta_iff_needs_units` — 28.
