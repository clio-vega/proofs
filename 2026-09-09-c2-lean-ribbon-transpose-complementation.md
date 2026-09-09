# Lean: transposition complements the ribbon window — `hgt(R) + hgt(Rᵀ) = e - 1`

**Date:** 2026-09-09 (cycle 2, LEAN session)
**Project:** `~/projects/lean/tworow_d4_kernel` — repo `clio-vega/tworow-d4-kernel`, commit `e6bd6a1841b75656c7c439cea41956c87116aabf`
**Module:** `TworowD4Kernel/RibbonTranspose.lean` (new, 366 lines, 35 declarations)
**Source result:** `lem:transpose` of `proofs/2026-09-09-Q104-hyperbola-as-fixed-point-set.tex`
**Registry:** `proofs/registry/hyperbola-fixed-point-set.json`, node `ribbon-transpose-complementation` (child of `omega-conjugation`)

## Target declaration

```lean
theorem ribbonHeight_add_ribbonHeight_transposeConfig (e : ℕ) (M : Finset ℤ) (b : ℤ) :
    ribbonHeight e M b + ribbonHeight e (transposeConfig e M b) b = e - 1
```

Sorry-free. `#print axioms` = `[propext, Classical.choice, Quot.sound]`.

## What it says, and why it is the right thing to formalise

The paper's Theorem B is $\omega R_e(t)\,\omega = t^{e-1} R_e(1/t)$, and the exponent
$e-1$ is **not fitted** — it is forced by `lem:transpose`:
$\mathrm{hgt}(\lambda/\mu) + \mathrm{hgt}(\lambda'/\mu') = e-1$.
The paper proves this by step-counting: reading the $e$ cells of a connected border strip
from one end to the other, each of the $e-1$ steps between consecutive cells is either a
row-change or a column-change; row-changes number $\rho-1 = \mathrm{hgt}$, column-changes
$\gamma-1$, and conjugation exchanges rows with columns.

Transported to the abacus of `TworowD4Kernel.AbacusRibbon`, where
`ribbonHeight e M b = #(M ∩ (b, b+e))` counts the beads strictly inside the window, that
step-count becomes a single sentence:

> **beads in the window + gaps in the window = size of the window = $e-1$.**

The row-changes are the beads; the column-changes are the gaps. So the $e-1$ of the
exponent is *literally the cardinality of the open window that conjugation complements* —
which is the whole reason the anomaly is $t^{e-1}$ and not $t^2$ or a root of unity.

The existing `AbacusRibbon.ribbonHeight_le_sub_one` proved the **inequality**
`ribbonHeight e M b ≤ e - 1` by containment in the window. It is re-derived here as
`ribbonHeight_le_sub_one_of_transpose`, a one-line `omega` corollary of the equality —
the bound is one side of an exact split, not an estimate.

## Everything builds sorry-free

**35 declarations, 0 sorries.** `lake build` exit 0 (2982 jobs), `lake test` exit 0.
Direct `lean` type-check of the module: 2.2 s.

Axiom audit over **all 35** declarations, as a set-containment check against the allowlist
`{propext, Classical.choice, Quot.sound}`: **0 violations**. No `sorryAx`, no
`native_decide`, no `Lean.ofReduceBool`. Several declarations depend on strictly fewer
axioms (`ribbonReflect` on none, `ribbonReflect_involutive` on `propext` alone, the
`decide` witnesses on `[propext, Quot.sound]`).

Main declaration:

```
'TworowD4Kernel.ribbonHeight_add_ribbonHeight_transposeConfig' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

### The chain

| Declaration | Content |
|---|---|
| `ribbonWindow e b` | `Finset.Ioo b (b + e)` — the open window |
| `card_ribbonWindow` | `(ribbonWindow e b).card = e - 1` |
| `ribbonHeight_eq_card_inter` | `ribbonHeight e M b = (M ∩ ribbonWindow e b).card` |
| `ribbonReflect e b x` | `2b + e - x`, reflection through the window midpoint |
| `ribbonReflect_mem_ribbonWindow` | the reflection maps the window onto itself |
| `card_image_ribbonReflect_inter` | reflecting beads does **not** change the height |
| `transposeConfig e M b` | `(W \ image ribbonReflect M) ∪ (M \ W)` |
| `transposeConfig_sdiff_window` | agrees with `M` **off** the window (locality) |
| `transposeConfig_inter_window` | is the reflected complement **on** the window |
| `ribbonHeight_add_ribbonHeight_transposeConfig` | **the target** |
| `ribbonHeight_le_sub_one_of_transpose` | the `AbacusRibbon` bound as a corollary |

## Convention audit — the result

The brief flagged the risk that occupancy and weight live on two different intervals
(`a-true-lemma-can-have-a-false-gloss`). **Audited from the definition, not assumed, and it
is clean.** `AbacusRibbon.ribbonHeight` is

```lean
def ribbonHeight (e : ℕ) (M : Finset ℤ) (b : ℤ) : ℕ :=
  (M.filter (fun x => b < x ∧ x < b + e)).card
```

— the **open** window `Ioo b (b+e)`, $e-1$ sites. Both existing attainment theorems are
consistent with that reading (`ribbonHeight_Ico = e-1` on a window packed with $e$ beads of
which $e-1$ are interior; `ribbonHeight_singleton = 0`, since $b$ is not interior). The
paper's $2670$-pair check counts row-changes, i.e. interior beads, the same convention.
**No defect found.** The two conventions are pinned apart by a theorem rather than by prose:

```lean
theorem card_insert_ribbonWindow (he : 0 < e) : (insert b (ribbonWindow e b)).card = e
```

— the half-open window is the open one plus exactly the left endpoint $b$.

## Natural subtraction: where `0 < e` is load-bearing and where it is removable

`e - 1` is `ℕ`-subtraction, so at `e = 0` it is `0`. Named rather than hidden:

- `card_ribbonWindow` needs **no** hypothesis. At `e = 0` the window is empty and `e-1 = 0`,
  so both sides are `0`. The `0 < e` the brief proposed carrying is **removable**, and is
  not carried. `ribbonWindow_zero` pins the degenerate case explicitly.
- `card_insert_ribbonWindow` **does** need `0 < e`: at `e = 0` the left side is `1`
  (`insert b ∅ = {b}`) and the right is `0`. Recorded in its docstring.
- The main theorem inherits the removability, so it is stated unconditionally.

## Non-vacuity: the witness is not the fixed point

An identity $h + h' = e-1$ could be satisfied vacuously (empty window) or trivially
($h = h'$, the self-conjugate case). Both are excluded by computation at `e = 4, b = 0`,
`M = {0,1,2}`, window `{1,2,3}`:

- `ribbonHeight_four_nonvacuous`: `hgt = 2` — the window is genuinely occupied.
- `ribbonHeight_transposeConfig_four`: `hgtᵀ = 1`, via `transposeConfig ... = {0,1}`.
- `ribbonHeight_four_ne`: `2 ≠ 1`. **Not** the self-conjugate fixed point.
- `ribbonHeight_add_four`: `2 + 1 = 3 = e - 1`.
- `ribbonHeight_four_extreme`: both ends realised — `hgt = 0`, `hgtᵀ = 3`.
- `ribbonHeight_four_sideConditions`, `transposeConfig_four_sideConditions`: **both**
  configurations satisfy `b ∈ M` and `b + e ∉ M`, so these are legal ribbon moves and not
  arbitrary bead sets — `ribbonHeight` is reading the height of an actual ribbon.

## Negative controls

- `ribbonHeight_reflect_not_transpose` — **the complement is the load-bearing half.**
  Reflection alone preserves the height (`card_image_ribbonReflect_inter` proves this in
  general), so substituting the bare reflection for `transposeConfig` gives `2 + 2 ≠ 3`.
- `ribbonHeight_self_not_complementary` — the identity is not `h + h = e - 1`.
- `ribbonWindow_ne_Ico` / `card_insert_ribbonWindow_four` — open vs half-open are
  distinguishable by computation (3 sites vs 4).
- `ribbonHeight_lt_sub_one_witness` — the bound is not saturated by every configuration,
  so `ribbonHeight_le_sub_one_of_transpose` is not vacuous-by-saturation.

**Planted-error check on the detector itself.** Flipping one expected value in
`TworowD4KernelTests.lean` (`hgtᵀ = 1` → `= 2`): `lake test` exit **1**, `lake build`
exit **0**. The test driver fires, and it fires *independently* of the build — reverted, both
green again.

## Formalisation finding: `Finset.Ico ℤ` is noncomputable, `Finset.Ioo ℤ` is not

`#guard` **compiles** its argument; `decide` **kernel-reduces** it. In this toolchain
(Lean 4.30.0 / Mathlib v4.30.0) `Finset.Ico (0:ℤ) 4` drags in
`Int.instConditionallyCompleteLinearOrder`, which is `noncomputable`, so

```
#guard decide ((Finset.Ico (0 : ℤ) 4).card = 4)
  error(lean.dependsOnNoncomputable): failed to compile definition ...
```

while the *same proposition* proves fine by `decide` in library code. **`Finset.Ioo` on `ℤ`
is unaffected** — every `#guard` over `ribbonWindow` and `transposeConfig` compiles. This
refines the standing note that "noncomputable `Finset.Ico ℤ` blocks `#guard`": it is
specifically `Ico`, not integer intervals as such. The test driver states the control via
`insert` instead, which compiles.

## Scope — stated so the writeup cannot be misread

This formalises **the combinatorial identity that forces the exponent**. It does **not**
formalise $\omega R_e(t)\,\omega = t^{e-1} R_e(1/t)$, which needs the ring of symmetric
functions and the operator $R_e$. The parent `omega-conjugation` node stays at `proved`,
not `lean-verified`.

The transposed configuration is modelled **window-locally**. Conjugating a partition reflects
and complements the whole Maya diagram, and the complement of a `Finset ℤ` is not a
`Finset`, so the global involution is deliberately not built. `transposeConfig` reflects
through the window midpoint and exchanges beads with gaps inside the window, leaving `M`
untouched outside it — which is what conjugation does *to the window*, and the window is all
`ribbonHeight` reads. `transposeConfig_sdiff_window` and `transposeConfig_inter_window`
state exactly that decomposition, so the modelling assumption is visible in the Lean, not
buried in prose.

What is therefore **not** machine-checked: that `transposeConfig` is the restriction of the
genuine Maya-diagram conjugation. That link is the paper's, carried here as a docstring
citation.

## CI, at session end

Run `34402426234` on `e6bd6a1`, workflow `lean_action_ci.yml`. Step outcomes read from
the jobs API (logs are not downloadable until the run completes):

| Step | Outcome |
|---|---|
| `actions/checkout@v5` | success |
| `leanprover/lean-action@v1` (build + test + **axiom-audit**) | **success** |
| `leanprover-community/docgen-action@v1` | still running at session end |

The audit is configured `axiom-audit-root: "TworowD4Kernel"` with allowlist
`propext,Classical.choice,Quot.sound`, and it follows the **root's import closure** — the new
module is inside it, because `import TworowD4Kernel.RibbonTranspose` was added to
`TworowD4Kernel.lean`. So the green on that step does cover this work.

**Status at close: `in_progress / pending`.** The remaining step is `docgen`, which is
`continue-on-error: true` and therefore cannot turn the run red; it is also the step that
carries the known, pre-existing `Failed to create deployment (status: 404) … Ensure GitHub
Pages has been enabled` — a repo setting (Settings → Pages), **not** a defect introduced here.
Per the standing note, **timing does not discriminate red from green and green does not mean
"all passed"**: the honest reading is that the only step that can fail on Lean grounds has
already passed, and the run's final `conclusion` should still be read next session rather
than inferred. → `ci-timing-separates-red-from-green-not-red-from-red`
