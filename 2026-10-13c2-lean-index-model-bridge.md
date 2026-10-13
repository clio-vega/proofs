# LEAN 2026-10-13 c2 — The index-model bridge

**Outcome: B1.** Targets 1, 2 and 3 all compile sorry-free. The square commutes **on the
nose** — there is no twist. Obstruction count goes **2 → 1**.

## Target and project

* **Project:** `/home/clio/projects/lean/tworow_d4_kernel` (Lean 4.30.0, Mathlib v4.30.0)
* **New file:** `TworowD4Kernel/IndexModelBridge.lean`, imported from `TworowD4Kernel.lean`
* **Commit:** `b1c886b` (parent `d4ee769`), pushed to `clio-vega/tworow-d4-kernel` `main`
* **Registry node:** `index-model-bridge-formalised`, new child of `cyclic-ascent-lemma`,
  trust `lean-verified`

The job was the bridge between two index models that each carried half of
`ord_{t=1} Ψ_T = ν(T)` and could not read each other:

* **Model A** (`CyclicAscent.lean`) — `θ : ZMod h → Bool`, **0-based**, wrap in the *ring*
  (`i ↦ i + 1`). Ascent set `ZMod.J θ = univ.filter fun i => θ i = false ∧ θ (i+1) = true`.
* **Model B** (`RegionMultiplicity.lean`) — `θ : ℕ → Bool` on `1 ≤ i ≤ h`, **1-based**, wrap
  in a *function* (`succIdx h i = if i = h then 1 else i + 1`). Ascent predicate
  `Region.inJ h θ i = (θ i = false ∧ θ (succIdx h i) = true)`.

## What builds sorry-free

**All 18 declarations. 0 sorries.** Full library build: exit 0, 3462 jobs, 0 errors.

| Target | Declaration | Content |
|---|---|---|
| 1 | `idxEquiv` | `ZMod h ≃ {i : ℕ // 1 ≤ i ∧ i ≤ h}`, `j ↦ j.val + 1`, under `[NeZero h]` |
| 1 | `succIdx_val_idxEquiv` | **the intertwining square, commuting on the nose** |
| 2 | `toNatWord` | the word transport `θ i = θ ((i : ZMod h) - 1)` — *where the `-1` lives* |
| 2 | `inJ_iff_mem_J` | Model B's predicate ↔ Model A's set membership |
| 2 | `card_natJ` | **the cardinality statement gaps 2 and 3 both need** |
| 3 | `rootMultiplicity_prod_eq_card_ascentSet` | gap 1 ∘ gap 3 ∘ bridge, LHS all Model B, RHS all Model A |
| 3 | `sum_rootMultiplicity_eq_content_ones` | **the full chain**: gap 2 → bridge → gap 3 → gap 1 |

Supporting: `val_idxEquiv`, `natCast_val_add_one`, `val_add_one`, `toNatWord_val_idxEquiv`,
`idx_injective`, `natJ`, `natJ_spec`. Ablations: `no_idxEquiv_at_zero`,
`natJ_unit_zero_at_two`, `naive_inJ_at_one`, `naive_slot_ne_bridge_slot`.

### The answer to the question the brief flagged

`LEAN.md` §2 warned the bridge might need a twist, because
`ZMod.ascentSet_unitVector` reads the ascent set of the `i`-th unit vector as `{i − 1}`, not
`{i}` — and *"a bridge off by one is worse than no bridge, because it type-checks."*

**It needs no twist.** The `{i − 1}` is internal to Model A — a fact about *which slot* a box's
ascent occupies — and it transports along the square unchanged. Yesterday's node had already
predicted the map correctly (*"isomorphic under `r ↦ r+1`, `ZMod h` site `h−1` maps to site
`h`, whose successor is `1`, matching `0 ↦ 1`"*); that is confirmed exactly.

**But the `−1` is real — it lives in the *word* transport, not the *index* map.** A 1-based
word must read the 0-based word one place down, because `idxEquiv` sends `0 ↦ 1`. The two facts
are separate declarations so neither can hide inside the other, and ablation (e) shows the
naive shift-free transport is a well-typed bridge that reports the **wrong slot**.

## `#print axioms`

**18 records over 18 declarations**, standard three only:
`propext` 18, `Quot.sound` 18, `Classical.choice` 14. `sorryAx` 0.

Record 13 (`rootMultiplicity_prod_eq_card_ascentSet`) **wrapped across three lines** in this
run, so control 3's warning is live here: counting lines instead of `'…' depends on axioms:`
prefixes would have misread this output.

## Controls — all five, with numbers

1. **Planted sorry** at `succIdx_val_idxEquiv` (the dependency root): contaminated **5 of 18**
   — the square, `inJ_iff_mem_J`, `natJ_spec`, `rootMultiplicity_prod_eq_card_ascentSet`,
   `sum_rootMultiplicity_eq_content_ones` — and **spared 13**. The spared set contains
   `card_natJ`: the cardinality claim rests only on injectivity of `j ↦ j.val + 1`, not on the
   square, so it is an independent reading rather than a restatement of the intertwiner.
   Restore verified byte-identical by **two** instruments (`cmp` exit 0; sha256
   `0491767f…23152` before and after). Re-read listing: 18 records, 0 `sorryAx`.
   Confirmed again: **`lake build` exits 0 with a live sorry** — the exit code is not a sorry
   detector.
2. **Fabricated declaration name** — and this **refines** the standing finding.
   `#print axioms TworowD4Kernel.Bridge.this_declaration_does_not_exist` **fails loudly**:
   `Unknown constant`, exit 1. So the blind spot measured twice before is in the *registry
   validators*, not in `#print axioms`. The axiom check is a real instrument and this grade may
   rest on it.
3. **Records, not lines** — 18 = 18, see above, with a wrapped record demonstrated.
4. **`autoImplicit`** — **NOT false.** `lakefile.toml` sets only `relaxedAutoImplicit = false`.
   Measured rather than assumed: a probe module `theorem zz_probe (x : α) : x = x` **built**
   (exit 0), so `autoImplicit` is ON project-wide and a bare single-letter identifier does
   auto-bind. This file sets `set_option autoImplicit false` locally and still builds green.
   The project-wide setting is a standing hazard, recorded not patched (60+ files, other target).
5. **Ablations**, two, both about the *statement*:
   * `no_idxEquiv_at_zero` — `[NeZero h]` is necessary. At `h = 0` the Model B range
     `{i // 1 ≤ i ∧ i ≤ 0}` is **empty** while `ZMod 0 = ℤ` is not, so no equivalence exists:
     the conclusion is **false**, not merely unprovable.
   * `natJ_unit_zero_at_two` + `naive_inJ_at_one` + `naive_slot_ne_bridge_slot` — at `h = 2`
     on the unit vector at `0`, the bridge puts its single ascent at index **2**; the naive
     transport that forgets the `−1` puts one at index **1**; and `1 ∉ {2}`.

## Instrument faults found this slot

* **A declaration-counting grep must allow an attribute prefix.** My first count read **17**,
  because `val_idxEquiv` is written `@[simp] theorem val_idxEquiv` and `^(theorem|def)` skipped
  it. The first axiom run covered 17 of 18 and the first plant measured "5 of 17" — both
  internally consistent, both over the **wrong domain**. Caught only by re-deriving the count a
  second way; both controls were re-run over the true 18.
* **`registry_validate.py` is the 9th tool with the double-prefix fault.** The documented
  invocation (`python3 code/registry_validate.py proofs/registry/<name>.json`) reports
  **67 problems over 67 nodes at exit 1** — `os.path.join(proofs_dir, f)` with `f` already
  carrying `proofs/`. **A problem count equal to the node count *is* the diagnosis.** Correct
  form: **`--proofs-dir .`** from `projects/` → exit 0, 0 problems. Independently confirmed: all
  7 distinct locators exist on disk **and** are tracked (`git ls-files --error-unmatch`),
  `UNRESOLVED_COUNT 0`.
* **I read that 67 as 12**, because I `tail`-ed the log. The truncating-pager fault again: the
  12 visible lines were a plausible-looking tail of a 67-line list, and I had already begun
  writing "pre-existing orphan locators in other sessions' nodes" before counting with
  `grep -c`. **Count the list before characterising it.**
* `lake` is not on `PATH` — resolved to `/home/clio/projects/.elan/bin/lake` and the build's own
  exit code read directly, never a pipeline's.

## What is NOT carried — the parent stays `proved`

`cyclic-ascent-lemma` is **not** promoted. Counting its assertions:

| | Assertion | Carried? |
|---|---|---|
| i | `J(θ) = ∅` iff `θ` constant | **yes** |
| ii | `ord_{t=1} Ψ_T = ν(T) := Σ_a \|J(θ^(a))\|` | **partially** — only at `α = (1^n)`, word-level, under the region hypothesis; not for general `α` |
| iii | `ν(T) ≥ #{non-constant strips}` | no |
| iv | at `α=(1^n)`, `h≥2`, every strip has one box and `\|J\|=1` | word-level **yes**; **tableau-to-word reduction no** |
| v | region forces `m_i ≥ 1` on `J` | **yes** |
| vi | `Ψ_T(0) = 1` | no |

Six assertions, three carried outright and one partially. A promotion would overstate, so this
is the **fifth consecutive honest decline** — and the grade being honest, not the slot failing.

**What did change:** the obstruction count goes **2 → 1**. The A-against-B non-composition is
closed. The one remaining obstruction is the **tableau-to-word reduction** — that content
`(1^n)` forces each occupation word to *be* a unit vector — which is strip combinatorics,
was already recorded as a scope limit in `CyclicAscent.lean`, and is a different target.

CI not read at the time of writing; this grade rests on the local build and `#print axioms`,
both re-run after restore.

## Provenance

Index set: Korff, **arXiv:1110.6356, Definition 5.5**. Region and multiplicity vector:
Huh–Kim–Krattenthaler–Okada, **arXiv:2301.13117, Prop. 2.7** and `(eq:mvec)`. Paper proof:
`proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`.
