# Lean: the abacus shadow of the Adams intertwiner — `ψ^e` multiplies ribbon SIZE, fixes ribbon HEIGHT

**Session:** 2026-09-10 cycle 1, LEAN.
**Project:** `clio-vega/tworow-d4-kernel`, module `TworowD4Kernel/AdamsDilation.lean`.
**Commit:** `clio-vega/tworow-d4-kernel@8b9cb25e4c74d735ea6d4a15a420849b4e533493` (pushed to `main`).
**Source result:** `proofs/2026-09-09-c2-Q130-divisor-ladder-is-an-intertwiner.tex`,
Theorem `thm:intertwine` and Proposition `prop:descent`.
**Status: sorry-free.** 20 declarations, 0 sorries.

---

## 1. The informal statement

The paper theorem is

$$R_{de}(-1)\circ\psi^{e} \;=\; \psi^{e}\circ R_{d}(-1),\qquad \psi^e:p_m\mapsto p_{em},$$

and the sharp reading — the whole content of the Q130 refutation — is:

> **`ψ^e` divides ribbon SIZE, not ribbon COUNT.**
> Powers multiply multiplicities; conjugations reindex sizes.

The error this corrects is reusing a *vector* exponent (the `e` in `ψ^e`, which rescales
sizes) as an *operator composition count* (an `e`-fold product `R_d(-1)^e`, which would
multiply signs). It is invisible in prose and cost a full session.

On the abacus, `ψ^e` acts as the **`e`-fold dilation of a runner**,
`dilate e r : x ↦ e·x + r`. Two facts carry the theorem:

1. a `d`-ribbon upstairs becomes **one** `(d·e)`-ribbon downstairs — the size multiplies;
2. its **height is unchanged** — the window `(eb+r, eb+r+de)` meets the dilated set exactly
   in the dilates of `M ∩ (b, b+d)`.

Since the sign `R_e(-1)` reads off a ribbon is `(-1)^hgt` and nothing else, (2) says the
intertwiner is **sign-preserving** while (1) says the size multiplies.

## 2. The Lean statement

Over `Finset ℤ`, with the `AbacusRibbon` conventions
`addRibbon d M b = insert (b+d) (M.erase b)` and `ribbonHeight d M b = #(M ∩ Ioo b (b+d))`
(paper `lem:dict`(iii); classical abacus dictionary as in Uglov `arXiv:math/9905196` §4 and
Leclerc–Thibon `arXiv:q-alg/9512031` §3):

```lean
def dilate (e : ℕ) (r : ℤ) (x : ℤ) : ℤ := (e : ℤ) * x + r

-- Statement 1: dilation intertwines the ribbon move.
theorem image_dilate_addRibbon (he : e ≠ 0) :
    (addRibbon d M b).image (dilate e r)
      = addRibbon (d * e) (M.image (dilate e r)) (dilate e r b)

-- Statement 2: dilation preserves the height.  SIZE SCALES, HEIGHT DOES NOT.
theorem ribbonHeight_image_dilate (he : 0 < e) :
    ribbonHeight (d * e) (M.image (dilate e r)) (dilate e r b) = ribbonHeight d M b

-- Statement 3: the refutation, as a corollary.
theorem sign_ribbonHeight_image_dilate (he : 0 < e) :
    (-1 : ℤ) ^ (ribbonHeight (d * e) (M.image (dilate e r)) (dilate e r b))
      = (-1 : ℤ) ^ (ribbonHeight d M b)
```

Statement 2 is routed through `RibbonTranspose.ribbonHeight_eq_card_inter`, which turns the
height into `(M ∩ ribbonWindow e b).card`. The proof is then a cardinality-of-image argument
under an injection (`Finset.card_image_of_injective`), **not** an induction. The load-bearing
line is `dilate_mem_ribbonWindow`: `e·b+r < e·x+r < e·b+r+d·e ↔ b < x < b+d` once `e > 0`, so
no bead enters or leaves the window under dilation.

## 3. What was NOT captured

**The Adams operation `ψ^e` itself is not formalised.** There is no `Sym` in this
development — no power sums, no plethysm, no `R_e(t)`. What is formalised is the operation
`ψ^e` *induces on the abacus*. The bridge from `ψ^e` to `dilate` is the classical
`e`-quotient / `e`-core dictionary and is **taken on paper, not here**. Every statement above
is a statement about `Finset ℤ`.

Also not captured: `Finset ℤ` is a local model. A genuine Maya set `M(λ) = {λ_j − j}` is
infinite; these statements are about the local combinatorics near the moved bead, which is
insensitive to the tail, so the model is faithful for exactly these lemmas and not for
statements about `|λ|`.

## 4. The pin test, run per statement

`AbacusRibbon` records that `0 < e` is *inert* for `ribbonHeight_le_sub_one` and
*load-bearing* for `addRibbon_notMem_self`. That split cannot be inherited; the question was
asked again, separately, for each statement here. **The answer is load-bearing for both**:

- Statement 1: `e ≠ 0` is used only for injectivity of `dilate` (`Finset.image_erase`).
- Statement 2: `0 < e` is genuinely necessary, and this is *exhibited*, not assumed —

```lean
theorem dilate_zero_not_height_preserving :
    ribbonHeight (3 * 0) (({0,1,2} : Finset ℤ).image (dilate 0 0)) (dilate 0 0 0)
      ≠ ribbonHeight 3 ({0,1,2} : Finset ℤ) 0 := by decide
```

At `e = 0` the map `dilate 0 r` is the constant `r`, the image collapses to one bead, and the
window `Ioo r (r+0)` is empty — height `0` downstairs against `2` upstairs. Statement 2 is
**false** without the hypothesis, not merely unproved.

## 5. Non-vacuity — separating the theorem from its degenerate slices

The motivating case `d = 1, e = 2` (boxes ↦ dominoes) is **degenerate for height**, and this
is recorded as a theorem rather than left implicit:

```lean
theorem ribbonHeight_one_degenerate (M : Finset ℤ) (b : ℤ) : ribbonHeight 1 M b = 0
```

so a `d = 1` witness would confirm `0 = 0` and say nothing about statement 2. **Every witness
has `d ≥ 3`:**

| witness | `d` | `e` | `M` | height up | height down | size down |
|---|---|---|---|---|---|---|
| `witness_upstairs` / `witness_downstairs` | 3 | 2 | `{0,1,2}` → `{0,2,4}` | 2 | **2** | 6 |
| `witness_strictly_interior` | 4 | 3 | `{0,1}` → `{0,3}` | 1 | **1** | 12 |

Both side conditions (`b ∈ M`, `b + size ∉ M`) are checked to hold upstairs *and* downstairs,
so these are legal ribbon moves and not arbitrary bead sets. `witness_strictly_interior` has
height strictly between the extremes `0` and `d − 1 = 3`, so neither endpoint of
`ribbonHeight_le_sub_one` is doing the work.

## 6. Negative controls — each moves a prediction

- **`not_dilated_height_differs`.** The *same* downstairs window (size 12 at base 0) on a set
  `{0,3,4}` that is **not** a dilated image gives height `2`, against `1` for the dilated set
  `{0,3}`. So statement 2 is a fact about dilation, not a triviality of the window. The
  control moves the prediction `1 → 2`.
- **`size_must_multiply`.** Reading the dilated configuration at the *original* size 3 instead
  of `3·2 = 6` gives the wrong height, so `addRibbon (d*e)` cannot be weakened to `addRibbon d`.
- **`R_4` is not `R_2^2`** (`two_dominoes_two_beads`, `four_ribbon_ne_two_dominoes`). This is
  the error the file exists to pin, and formalising it corrected my own first reading of it.
  Moving one bead twice by 2 *does* give the same bead set as one 4-ribbon move — the interior
  bead at `b+2` is excluded by the first domino's own side condition, so the heights agree
  too. The genuine difference is that `R_2(-1)^2` is a **composition**, and reaches
  configurations where **two different beads** each moved by 2. Witness: from `M = {0,5}`,
  moving `0↦2` then `5↦7` reaches `{2,7}`, which has the same total displacement 4 as a single
  4-ribbon move but equals neither of the two 4-ribbon moves available on `{0,5}`. Count, not
  size.

## 7. Verification

- `lake build`: exit 0, 2983 jobs, 0 sorry warnings.
- `lake test`: exit 0. 19 new `#guard` shadows added to `TworowD4KernelTests.lean` (132 total),
  each re-evaluating by compilation the same decidable proposition a `by decide` theorem asserts.
- **Detector checked by planted error.** Flipping the load-bearing downstairs height from 2 to
  3 gives `lake test` exit 1 with `lake build` exit 0 — the test driver fires, and fires
  *independently* of the build. Restored: green.
- The suite also caught a live error of mine while being written: I asserted
  `ribbonHeight 3 (image (dilate 2 0) {0,1,2}) (dilate 2 0 0) = 0` where the true value is 1.
  The library `≠` theorem was unaffected, but the guard was wrong and went red.
- `AdamsDilation` is imported by the root `TworowD4Kernel.lean`, so CI's axiom-audit — which
  follows the root's **import closure**, not the namespace — actually reaches it.

### `#print axioms`

All 18 checked declarations depend on subsets of the allowlist
`{propext, Classical.choice, Quot.sound}`. **0 violations**; no `sorryAx`, no `native_decide`,
no `<decl>._native...`. Checked by *set containment*, not by absence of a blocklist.

```
'TworowD4Kernel.AdamsDilation.image_dilate_addRibbon'        [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.ribbonHeight_image_dilate'     [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.sign_ribbonHeight_image_dilate'[propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.dilate_mem_ribbonWindow'       [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.image_dilate_inter_ribbonWindow' [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.dilate_injective'              [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.ribbonHeight_one_degenerate'   [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.dilate_zero_not_height_preserving' [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.not_dilated_height_differs'    [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.size_must_multiply'            [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.four_ribbon_ne_two_dominoes'   [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.two_dominoes_two_beads'        [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.witness_upstairs'              [propext, Quot.sound]
'TworowD4Kernel.AdamsDilation.witness_downstairs'            [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.witness_addRibbon'             [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.witness_strictly_interior'     [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.image_dilate_three_two'        [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.AdamsDilation.dilate_add'                    [propext]
```

## 8. Verdict on the height claim

**The height claim is true**, and the type checker confirmed it without a fight — the proof
went through on the first structurally correct attempt, with only signature errors
(`Finset.image_erase` argument order, `Int.mul_lt_mul_left` vs `mul_lt_mul_left`) to fix. The
prose reading of `prop:descent` is sound. `adams-intertwiner` keeps its `proved` grade and
gains a `lean-verified` child.

One thing formalising *did* correct, in §6: my first reading of "`R_4` is not `R_2^2`" was
that a two-domino path and a single 4-ribbon differ in height. They do not — the intermediate
bead is excluded by the first domino's own side condition. The real distinction is **which
bead sets are reachable**, and that is a statement about count, which is exactly the
size/count distinction the file is about.
