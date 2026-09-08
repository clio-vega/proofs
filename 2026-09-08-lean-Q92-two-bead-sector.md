# LEAN 2026-09-08 — the two-bead sector of `[R_e(t), R_f(s)]` on the locus `ts = 1`

**Project:** `clio-vega/tworow-d4-kernel`, local `/home/clio/projects/lean/tworow_d4_kernel`
**New module:** `TworowD4Kernel/CrossRankTwoBead.lean` (imports `TworowD4Kernel.CrossRankOneBead`)
**Paper proof:** `proofs/2026-09-07-c2-Q96-order-over-sublattice.tex` §4, Thm 4.1(2) and Cor 4.2(iii)
**Registry:** `proofs/registry/fock-ribbon-sign-operator.json`, new node
`Q96-two-bead-sector-ts-equals-one-lean` (`lean-verified`)

## Target declaration

```lean
theorem twoBeadSector_commutator_eq_zero_of_mul_eq_one {R : Type*} [CommRing R] (t s : R)
    (hts : t * s = 1) (he : 0 < e) (hf : 0 < f) :
    twoBeadSector f e M M' s t - twoBeadSector e f M M' t s = 0
```

`twoBeadSector e f M M' t s` is the two-bead part of `⟨M'|R_f(s) R_e(t)|M⟩`: the sum, over
routes whose four sites `x, y, x+e, y+f` are pairwise distinct, of
`t ^ ribbonHeight e M x * s ^ ribbonHeight f (addRibbon e M x) y`. Two parameters throughout;
`s = t` is the one-parameter case.

## Status: sorry-free

**0 sorries.** 14 theorems, 2 defs, 1 `Decidable` instance. `lake build` exit 0 (2980 jobs),
`lake test` exit 0 (2980/2980).

The three statements the brief asked for all landed:

1. **`ribbonHeight_addRibbon_shift`** (brief's `heightShift`) —
   `(ribbonHeight f (addRibbon e M x) y : ℤ) = ribbonHeight f M y + crossIndex e f x y`.
   Built from two smaller lemmas, `ribbonHeight_erase` and `ribbonHeight_insert`, rather than
   one insert/erase-and-count argument; each is stated in `ℤ` so no truncated subtraction
   appears anywhere.
2. **`twoBead_height_sum`** — `P' + Q' = P + Q` **in `ℕ`**, where `P'` and `Q'` are the heights
   of the moves taken *second*. This replaces the brief's `twoBead_contribution`: I took the
   brief's own warning seriously and never wrote `t^(P-k)`. `P'` and `Q'` are defined as the
   cardinalities they are; `Q' = Q + k` and `P' = P - k` are then *theorems*
   (`ribbonHeight_addRibbon_shift`, `ribbonHeight_addRibbon_shift_swap`), and the `k` cancels
   out of the sum. Nonnegativity is automatic rather than a hypothesis.
3. **`twoBeadSector_commutator_eq_zero_of_mul_eq_one`** — Cor 4.2(iii), forward direction.

## The thing I did not expect

The paper proves Cor 4.2(iii) by *reading off the factor* `1 - (ts)^k` from each summand, which
requires knowing what the summands are — hence the classification (one legal assignment when
`e ≠ f`, two when `e = f`) that the brief listed as stretch goal 4.

The Lean proof does not need it. Reversing the order of the two moves is an involution
`Prod.swap` on routes (`mem_routes_swap`: `(x,y) ∈ routes e f M M' → (y,x) ∈ routes f e M M'`),
and on `ts = 1` the two orderings carry the *same monomial*, because `P' + Q' = P + Q` and
`t, s` are inverse units. So the difference of the two sums cancels **termwise**, under a
bijection, without anyone counting the terms.

That means the `e = f` case — the one `Q92 Thm 2.3(2)` got wrong — is covered with **no case
split at all**. The correction the brief wanted handed to a type checker turns out to be a
place where the case split itself was the mistake: a symmetry between two assignments is one
lemma at two instantiations, not two lemmas. Stretch goal 4 is therefore not merely unfinished;
it is not needed for 1–3, which is a better outcome than finishing it would have been.

## Scope — what is NOT proved

- **The converse of Cor 4.2(iii).** The corollary says the sector vanishes identically *if and
  only if* `ts = 1`. Only the **if** direction is a Lean theorem. The **only if** is witnessed
  solely by `#guard`s exhibiting nonvanishing at `ts = 2`. For that reason
  `Q96-ts-equals-one-locus`, which states the iff, is deliberately **left at `proved`** rather
  than promoted.
- **The classification of legal assignments** (Thm 4.1(2), last sentence). Not formalised — see
  above; the main theorem does not use it.
- **The one-bead sector of Thm 4.1(1)** (the two-parameter version). Not formalised. The
  *one-parameter* one-bead sector was done on 2026-09-07 (`matrixElem_one_bead`).
- **The `e = f` correction** is kernel-checked as an evaluated counterexample, not as a named
  theorem: `e = f = 3`, `M = {0,1,2}`, `M' = {1,3,5}` gives `s²t - st² - s + t = (t-s)(1-st)`,
  `= 28` at `t=3, s=5` and `= 0` at `s=t`. This reproduces the paper's Remark 4.3 example
  independently (the paper reached it from `b=6, c=8` in a 9-bead window).
- `Q92-structural-divisibility` is **not** touched: that node is the one-bead `(1+t)`, a
  different statement.

## Guards — 33 `#guard`s, chosen to move

| instance | `e,f` | `M'` | `k` | commutator two-bead part |
|---|---|---|---|---|
| 1 | 2,3 | `{0,3,5}` | `+1` | `1 - st` |
| 2 | 2,3 | `{1,3,4}` | `-1` | `s²t - s` |
| 3 | 1,3 | `{0,3,4}` | `0` | `0` — **labelled in the file as a kernel** |
| 4 | 3,3 | `{1,3,5}` | `+1, -1` | `(t-s)(1-st)` — the correction |

`M = {0,1,2}` throughout. Instance 3 is identically zero for *all* `t, s`, so it cannot
distinguish `ts = 1` from anything; it is there to pin the third value of `k` and says so in
the file rather than being counted as evidence.

Each `ts = 1` vanishing is checked over `ℚ` at `t = 2, s = 1/2` — so `t ≠ s`, off the
`t = s = -1` anchor — and **paired with the same instance at `s = 1`** (`ts = 2`), where it does
*not* vanish.

Both hypotheses of `crossIndex_swap` are shown load-bearing by `#guard`ed counterexamples:
dropping `x ≠ y` (at `x=y=0, e=1, f=2`: `k=1`, `k'=0`) and dropping `x+e ≠ y+f` (at
`x=1, y=0, e=1, f=2`: `k=-1`, `k'=0`). The paper states exactly these two and no more.

**Negative control run live:** flipping the instance-4 guard from `28` to `0` turns `lake test`
red (exit 1); restored, exit 0.

## `#print axioms`

Every one of the 14 theorems lies inside `[propext, Classical.choice, Quot.sound]`, with the
main declaration exactly equal to those three:

```
'TworowD4Kernel.twoBeadSector_commutator_eq_zero_of_mul_eq_one'
    depends on axioms: [propext, Classical.choice, Quot.sound]
```

Two are strictly smaller: `TwoBead.swap` depends on *no* axioms, and `pow_swap_of_mul_eq_one`
on `[propext, Quot.sound]`. A grep sweep of the full output for `sorryAx`, `ofReduceBool`,
`_native`, `nativeDecide` returns nothing. The module is imported from the root
`TworowD4Kernel.lean`, so it is inside the CI axiom-audit import closure.

Full output for all 14: `ribbonHeight_erase`, `ribbonHeight_insert`,
`ribbonHeight_addRibbon_shift`, `crossIndex_swap`, `addRibbon_comm`, `TwoBead.swap`,
`routes_snd_legal`, `mem_routes_swap`, `twoBeadRoutes_swap`,
`ribbonHeight_addRibbon_shift_swap`, `twoBead_height_sum`, `pow_swap_of_mul_eq_one`,
`twoBead_contribution_of_mul_eq_one`, `twoBeadSector_commutator_eq_zero_of_mul_eq_one`.

## A note on the validator

`python3 code/registry_validate.py proofs/registry/fock-ribbon-sign-operator.json` on its bare
default reported **189** `file ... not found` violations, including
`proofs/2026-09-06-Q84-order-of-N-e.tex`, which exists. The default `--proofs-dir` is the
parent of the registry's directory (`proofs/`), and the node paths are themselves `proofs/...`,
so it looks under `proofs/proofs/`. With `--proofs-dir /home/clio/projects` the count is **2**,
both of them mine and both real. A wall of not-founds on files that exist is a root mismatch,
whichever binary printed it.
