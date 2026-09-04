# Lean snapshot — `TworowD4Kernel.AbacusRibbon`

**Date:** 2026-09-03 (cycle 2, LEAN session)
**Project:** `clio-vega/tworow-d4-kernel`, local `/home/clio/projects/lean/tworow_d4_kernel`
**Commit:** `6fb6c4a` (pushed)
**Module:** `TworowD4Kernel/AbacusRibbon.lean`, imported from the root `TworowD4Kernel.lean`
**Paper proof:** `proofs/2026-08-31-Q59-commutator-rigidity.tex`, Definition (lines 77–85) and
Lemma `lem:dict`(iii). External: Uglov `arXiv:math/9905196`; Leclerc–Thibon `arXiv:q-alg/9512031`.

## Why this target

Yesterday's DREAM note recorded `N_e^{(h)}` as **removing** an `e`-ribbon and read off
`R_e(-1) = p_e^⊥`. Both primary sources say it **adds** — the Definition in the Q59 paper, and
`fock_ell.py::R_e`, which moves a bead `b ↦ b + e`. The claim was wrong in prose while both
sources were correct, and it survived a full cycle. This session turns the direction into a type.

## Result: sorry-free

| | |
|---|---|
| `lake build` | **2974 jobs**, `Build completed successfully` |
| sorries | **0** (`grep -c sorry` = 0) |
| new declarations | **26 user-facing** (44 including compiler-internal) |
| audited declarations in the root's import closure | **591** total / **162** user-facing |
| axiom offenders | **0** |
| `sorryAx` occurrences | **0** |

### `#print axioms`, compared as a set

The audit was not eyeballed. A `run_cmd` over `env.constants` walked every declaration in the
`TworowD4Kernel` namespace **that is in the root's import closure**, called `Lean.collectAxioms`,
and filtered against the allowlist `{propext, Classical.choice, Quot.sound}`:

```
AUDITED_DECLS_ALL        591
AUDITED_DECLS_USERFACING 162
OFFENDERS                0
SORRY_DECLS              0
```

Zero offenders means every audited declaration's axiom set is a **subset** of the standard three;
no declaration reaches outside it. (`ABACUSRIBBON_ALL 44`, `ABACUSRIBBON_USERFACING 26` for this
module alone.)

**One number in the brief does not reproduce, and I am not going to round it.** The brief expected
the audited count to rise "from 549". My traversal reads 591 total, of which this module
contributes 44, leaving 547 for the prior closure — not 549. The two counts are produced by
different tools (mine is `env.constants` filtered by `Name.isPrefixOf`; the 549 came from CI's
lean-action), so they are not directly comparable and I did not re-run the old measurement. The
delta I can vouch for is **+44**.

## What is proved

Definitions, on `Finset ℤ` — the `WindowBridge` bead convention, reused rather than reinvented:

- `addRibbon e M b = insert (b + e) (M.erase b)` — intended domain `b ∈ M`, `b + e ∉ M`
- `removeRibbon e M c = insert (c - e) (M.erase c)` — intended domain `c ∈ M`, `c - e ∉ M`
- `ribbonHeight e M b = (M.filter (fun x => b < x ∧ x < b + e)).card`

Theorems:

1. **Inverse pair.** `addRibbon_removeRibbon`, `removeRibbon_addRibbon` — mutually inverse on
   their domains.
2. **Legality of the round trip.** `addRibbon_mem_self` (`b + e` occupied after adding),
   `addRibbon_notMem_self` (`b` vacated; this is the one that needs `0 < e`).
3. **Orientation.** `sum_addRibbon`: `∑ (addRibbon e M b) = ∑ M + e`. `sum_removeRibbon`:
   `∑ (removeRibbon e M c) = ∑ M - e`. This is the abacus shadow of `|μ| = |λ| + e`, and it is the
   only statement in the file that the opposite convention **cannot** satisfy — the sign of `e` is
   the direction. It is what pins `R_e(-1) = M_{p_e}` rather than `p_e^⊥`.
4. **Height is read off untouched beads.** `ribbonHeight_addRibbon`.
5. **Degree bound.** `ribbonHeight_le_sub_one` (`≤ e - 1`), `ribbonHeight_lt` (`< e` for `0 < e`) —
   the `deg_t R_e(t) ≤ e - 1` the paper asserts.
6. **Attainment, for every `e`, not just bounded.** `ribbonHeight_Ico` gives exactly `e - 1` on
   `Finset.Ico b (b+e)`; `ribbonHeight_singleton` gives exactly `0` on `{b}`. Side conditions
   discharged separately (`*_sideConditions`) so both witnesses are genuine ribbons.
7. **Non-vacuity** for every hypothesis set: `addRibbon_removeRibbon_nonvacuous`,
   `ribbonHeight_nonvacuous` (an *intermediate* height, `e = 3`, `h = 1` — neither endpoint),
   `ribbonHeight_Ico_nonvacuous`, `sum_addRibbon_nonvacuous`.
8. **Deliberate negative.** `exists_addRibbon_not_removeRibbon_of_mem` and `addRibbon_collision`:
   drop `b + e ∉ M` and the round trip fails — `e = 2`, `M = {0,2}`, `b = 0` collapses to `{0}`.
   A file that can only state the true direction has not tested anything.

## Two things the formalisation found that the paper proof does not say

Both are hypotheses that turn out to be **inert**, which is the kind of thing a type checker
notices and prose does not.

- **The inverse lemmas need no positivity.** `addRibbon_removeRibbon` and `removeRibbon_addRibbon`
  follow from `Finset.erase_insert` and `Finset.insert_erase` alone. `0 < e` enters exactly once,
  in `addRibbon_notMem_self` — and that is the load-bearing place, because without it the "removal"
  of `addRibbon e M b` at `b + e` is not a legal removal, merely an algebraic identity.
- **The height bound needs no ribbon side conditions.** `ribbonHeight_le_sub_one` holds for every
  `M` and `b`: it is `M ∩ (b, b+e) ⊆ Ioo b (b+e)` against `Int.card_Ioo`. So `b ∈ M`, `b + e ∉ M`
  are not used. This yields a **second, independent** route to the `h ≤ e - 1` half of the
  registry node `Q63-exponent-range-and-cross-runner` — by interval cardinality, where
  `Q63-window-bridge` argues by `ZMod N` labels.

## What is NOT formalised — the boundary of the module

- **`lem:dict`(iii) itself.** The bijection between `e`-ribbons addable to the Young diagram of `λ`
  and beads `b ∈ M(λ)` with `b + e ∉ M(λ)`, and the claim that the ribbon's cells have contents
  `b+1, …, b+e`, are statements about partitions and diagrams. Here the bead side is taken as the
  *definition*. Closing this needs Young diagrams given a type.
- **`M(λ)` itself.** A genuine Maya set is cofinite in `ℤ_{<0}`, hence infinite; these are
  `Finset`s. Every statement above is local to the moved bead and insensitive to the tail, so the
  `Finset` model is faithful for exactly these lemmas — and would **not** be for a statement about
  `|λ|`. `sum_addRibbon` is a shadow of `|μ| = |λ| + e`, not a proof of it.

## Registry

`proofs/registry/fock-ribbon-sign-operator.json`, new node
`root/Q75-multiplication-defect/Q75-anchor-MN/abacus-ribbon-move-direction`, `trust:
lean-verified`, `lean` field listing all 24 top-level declarations, `sources: [math/9905196,
q-alg/9512031]`. `registry_validate.py`: clean.

---

# Job 2 — turning the `test` check back on

`lakefile.toml` declared no `test_driver`, so lean-action reported the `test` check as
**skipped** on every run, which in the badge is indistinguishable from a passing check. Only
`axiom-audit` was detecting anything. Commit `5b2b038`.

## What was built

A test driver, `TworowD4KernelTests` (`testDriver = "TworowD4KernelTests"`), holding 15
`#guard`s on the *same decidable propositions* as the `by decide` witness theorems — the three
`WindowBridge` witnesses the brief named, plus the twelve new `AbacusRibbon` ones.

Two design decisions, both load-bearing:

- **A `lean_lib`, not a `lean_exe`.** I tried the executable first. A `lean_exe` driver has to
  link the entire Mathlib import closure: `lake build tests` reached **5931/5934** targets and
  **8101** generated `.c` files before I stopped it (timed out at 2 min, again at 90 s). Worse,
  `lake exe cache get` ships **oleans only**, so CI would compile all of Mathlib to object code
  on *every run*. This is the "needs more scaffolding than the session has" case the brief
  anticipated, and it has a clean way out.
- **Deliberately NOT in `defaultTargets`.** So `lake build` does not build it. This is what makes
  the check *independent*: a broken guard turns `test` red while `build` stays green. A detector
  that can only fire when another detector has already fired carries no information beyond it.

Cost: `lake test` is **3 s** on a warm tree.

## Ranked by planting an error — locally, confirmed

I planted **the DREAM note's own error**: the orientation guard
`∑ (addRibbon 3 {0,1} 0) = 4` changed to `= -2`, i.e. asserting the bead sum *falls* by `e` on
adding a ribbon.

```
lake build   →  Build completed successfully (2974 jobs)          [GREEN]
lake test    →  error: TworowD4KernelTests.lean:66:0: Expression   [RED, exit 1]
```

**This is the step I actually saw fail**, and it fired on `test` alone with `build` green — which
is the property the brief asked to be demonstrated, not merely "red overall". Reverted on `main`.

## Ranked in CI — CONFIRMED, and it took two runs

Runs on this repo turned out to take **~1m40 with a warm Mathlib cache**, not the ~44 min I had
assumed from the 2026-09-03 07:31 run, so the CI half fitted after all.

**First attempt failed for the wrong reason, run `33817046570` (`main`).** `lake test` died with

```
error: no such file or directory (error code: 2)
  file: .../tworow-d4-kernel/TworowD4KernelTests
```

The `TworowD4KernelTests.+` glob requires the directory `TworowD4KernelTests/` to exist. It existed
**locally** — left behind when I switched the driver from a `lean_exe` to a `lean_lib` — and git
never tracked it, being empty. Same family as `lake-glob-plus-excludes-root`: the glob resolves
against the working tree, not the commit. Fixed in `ab4c3c6`; the driver is a single file and the
glob is now `["TworowD4KernelTests"]`.

That run was not a loss: **it proved the wiring.** The step reported
`test … outcome=failure`, **not `skipped`** — lean-action is now genuinely running `lake test`.

**The real ranking, run `33817248757`, branch `ci-test-negative-control`.** Same planted error
(`= 4` → `= -2`). Step outcomes from the log:

```
build        outcome=success   "Build completed successfully (2974 jobs)"
test         outcome=failure   error: TworowD4KernelTests.lean:66:0: Expression
axiom-audit  outcome=skipped
```

**Red on the `test` check specifically, with `build` green in the same run.** That is the property
the brief asked for, and it is now confirmed on CI and not only locally.

### A third thing, which I was not looking for

`axiom-audit … outcome=skipped` in that log is not a configuration error — the steps run in
sequence, and a failing `test` **short-circuits every step after it**. So **a failing test disables
the axiom audit.** The two detectors are in *series*, not in parallel. Nothing is wrong today
(a red run is a red run either way), but it means one cannot read "axiom-audit skipped" as
"axioms are fine" — and if `test` ever becomes flaky, the axiom audit silently stops running
behind it. Worth a `continue-on-error` on `test`, or reordering, but that is a workflow change I
did not make in a session already carrying one.

### Status of `main` at session end — stated exactly

`main` is at `ab4c3c6`. Locally, on that exact tree: `lake build` **green (2974 jobs)** and
`lake test` **green (3 s)**; both run before the push. Its CI run **`33817241928` was still
`in_progress` after 18 minutes** when the session ended, and GitHub does not serve logs for a
running job, so **I did not see `main` go green in CI.**

The likely cause is benign and is itself a correction: the sibling run `33817248757` finished in
**1m38** because a failing `lake test` short-circuits everything after it — including the
**docgen** step. A *passing* run reaches docgen, which builds documentation for the whole Mathlib
closure. **So the "~44 minute run" figure I have been carrying is docgen, not Lean.** The Lean part
of this CI is under two minutes.

Next session: confirm `33817241928` (or a later `main` run) concluded `success` with
`test outcome=success` and `axiom-audit outcome=success`, then delete branch
`ci-test-negative-control`.

## Incidental finding

`Finset.Ico (0 : ℤ) 3` cannot appear inside `#guard`: `#guard` *compiles* its argument, and the
`LocallyFiniteOrder ℤ` instance is noncomputable in Mathlib v4.30.0
(`instConditionallyCompleteLinearOrder`). Kernel reduction inside `decide` is unaffected, so
`ribbonHeight_Ico` and its corollaries are fine — this is a compilation restriction, not a
mathematical one. The guard uses `{0,1,2}`, and the identification
`Finset.Ico (0:ℤ) 3 = {0,1,2}` is itself checked by `decide` rather than asserted.

## Rider, unchanged and still Robin's call

Run `33728508057` concluded `success` while carrying hard `X Error` annotations from
`actions/deploy-pages` (`status: 404`, "Ensure GitHub Pages has been enabled") — the docgen step
errors on every run and the job passes anyway, because of the deliberate `continue-on-error: true`
in `lean_action_ci.yml`. I did **not** touch it. Either enable Pages on
`clio-vega/tworow-d4-kernel` or drop the docgen step; it is a repository setting, not a Lean
question. Also still open: `actions/cache`, `actions/deploy-pages` and `actions/upload-artifact@v4`
are pinned to Node 20 and force-run on Node 24.
