# CI detector repair — making a green badge exclude a `sorry`

**Date:** 2026-09-01
**Session type:** LEAN (detector repair; no new formalisation, by instruction)
**Project:** `clio-vega/tworow-d4-kernel` (`/home/clio/projects/lean/tworow_d4_kernel`)
**Target:** not a declaration — the CI workflow `.github/workflows/lean_action_ci.yml`

## The confirmed defect

Run `33370053897` (branch `ci-negative-control`, 2026-08-31, 21m54s) built

```lean
theorem ci_negative_control : 2 + 2 = 5 := by sorry
```

and finished **`completed success`**. The log carried only

```
warning: TworowD4Kernel/NegativeControl.lean:8:8: declaration uses 'sorry'
end-action id=__leanprover_lean-action.build;outcome=success
```

Predicted outcome 2 from the 08-31 note, arrived exactly: **a green badge on this repo did
not exclude a `sorry`.** The mechanism is not subtle — `sorry` is a *warning* to `lake build`,
never an error, so the build step cannot be the detector, and no amount of watching it go
green will change that.

## Correction to the previous diagnosis

The 08-31 note and the session brief both blamed the log lines `nanoda-allow-sorry: true` and
`axiom-audit;outcome=skipped`. That diagnosis was wrong, and it matters because acting on it
would have produced a fix that changed nothing.

Reading `leanprover/lean-action`'s `action.yml` directly:

| input | default | gate |
|---|---|---|
| `nanoda` | `"false"` | step runs `if: inputs.nanoda == 'true'` |
| `nanoda-allow-sorry` | `"true"` | only consulted *when nanoda runs* |
| `axiom-audit` | `"false"` | step runs `if: inputs.axiom-audit == 'true'` |
| `axiom-audit-allow` | `"propext,Classical.choice,Quot.sound"` | — |

`nanoda` itself defaults to false, so the nanoda step **never ran**; `nanoda-allow-sorry: true`
was an unused input default being echoed into the log, not the cause. Setting it to `false`
— the brief's first-preference fix — would have left the badge exactly as blind as it was.
`axiom-audit` was `skipped` for the same reason: it defaults off.

The workflow had `- uses: leanprover/lean-action@v1` with **no `with:` block at all**. Every
check the action offers was off by default. That is the whole bug.

## The repair (commit `a7225b0`)

```yaml
      - uses: leanprover/lean-action@v1
        with:
          axiom-audit: true
          axiom-audit-allow: "propext,Classical.choice,Quot.sound"
          axiom-audit-root: "TworowD4Kernel"
```

`axiom-audit` (pinned by the action to `leanprover-community/axiom-audit` `v0.1.2`, commit
`46024e00`, verified after clone) fails the run if any declaration under the root namespace
transitively depends on an axiom outside the allowlist.

This is the right instrument for a reason that closes the loop with Job 3: **it is an
allowlist, not a blocklist**, which is the only sound shape here. `sorryAx` and the
per-declaration `<decl>._native.native_decide.ax_N_M` that Lean 4.30.0 generates for
`native_decide` cannot be enumerated in advance, so no blocklist can catch them. The allowlist
configured is the same standard three that every `#print axioms` in this project is asserted
equal to — so CI now enforces in the loop exactly the property the registry claims by hand.

`nanoda` was left off deliberately: it re-typechecks the environment with an external Rust
checker and its cost on a Mathlib-dependent project is unmeasured here. A guard that makes the
badge permanently red is exactly as useless as one that never fires (the docgen lesson from
08-31), so it is not being switched on unverified. `axiom-audit` is scoped to the project root
namespace and is the targeted check.

## A blind spot in the repair, found and closed (commit `5fea80f`)

`axiom-audit` takes **one** `--root`. `TworowD4Kernel/QuantumInteger.lean` declared
`namespace Clio.QuantumInteger` while all sixteen other modules are under `TworowD4Kernel`.
So the repair as first pushed would have silently skipped `C4_coefficient_identity` — **one of
the six `lean-verified` grades** — while reporting success. A detector with a blind spot over
the thing it is meant to check is the bug being fixed, not a cosmetic issue.

Renamed to `TworowD4Kernel.QuantumInteger`; rename only, no proof content touched. The module
is a leaf (nothing imports its declarations), verified by grep before the change. Locally:

```
Build completed successfully (2942 jobs).
'TworowD4Kernel.QuantumInteger.C4_coefficient_identity' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

Registry `lean` field updated to the new declaration name in
`proofs/registry/fock-ribbon-sign-operator.json`.

## Job 3 — `lean_evidence` audit: no change needed

All **six** `lean-verified` nodes (`C4-coefficient-identity` in
`fock-ribbon-sign-operator.json`; `c1`, `c2`, `c3-boundary`, `c4-number-lemma`, `c4-boundary`
in `three-row-even-jstar.json`) already carry a `lean_evidence` field asserting **equality**
with `{propext, Classical.choice, Quot.sound}`, and each already states the allowlist rule
explicitly. The 08-31 c2 session did this correctly. Nothing to fix — recorded because a
verified-and-already-correct check is a result, not a non-event.

## Grades unchanged

No trust grade moved this session. The cold build on 08-31 was clean (2971 jobs, zero `sorry`
warnings, all `#print axioms` exactly the standard three). **What was broken is the detector,
not the claim** — that distinction has now held for seven days and is worth saying out loud.

`registry_validate.py` on `fock-ribbon-sign-operator.json` reports one violation:
`Q63-level-ell-telescope` claims `proved` while child `Q63-rigidity-dichotomy` is `computed`.
This is **pre-existing and untouched** — my edit is a one-line declaration-name string, which
cannot cause a parent/child trust-boundary violation. Resolving it means changing a trust
grade, which is a mathematical judgement and out of scope for a Lean session. Left open and
flagged.

## Verification — red on a planted fault

The 08-31 docgen fix was correct and both `main` runs went green afterwards, and the badge was
*still* not a Lean detector. **Green after a fix is not evidence. Red on a planted fault is.**
So the repair is verified with the same control that exposed the defect.

### The first control was malformed, and that is itself a finding

Run `33550681860` (branch rebuilt on the fixed workflow, `sorry` planted in
`namespace TworowD4Kernel`) was **cancelled before completion, deliberately** — checking the
tool's source while it ran showed it would have gone **green for a reason having nothing to do
with the detector**.

`axiom-audit` builds its environment by **importing the root module**. That is exactly why the
tool carries a `--modules-from <dir>` flag, documented "for a root that does not transitively
import the whole library." `NegativeControl.lean` was compiled by the `TworowD4Kernel.+`
lakefile glob but was **not imported by the root aggregator**, so the planted `sorry` sat
outside the audited environment. Being in the right *namespace* is not sufficient; the module
must be *imported*.

Had I not checked, I would have read that green run as "the repair failed" and withdrawn a
claim that was in fact fine — the mirror image of the original bug. **A negative control that
cannot fail tests nothing.** Corrected in the control branch by importing the module from the
root (`33551104884`).

### Coverage check on `main` — the audit does see all six

The gap above is a general one: `lake build` coverage comes from the lakefile glob, but
`axiom-audit` coverage comes from the root module's **hand-maintained import list**. These are
different sets, and commit `c952afb` exists because that list had already silently omitted five
modules once. So it was checked rather than assumed — import closure of `TworowD4Kernel`
computed locally:

```
closure size: 18 (root + 17 modules)
NOT reachable from root: none
```

All six `lean-verified` declarations live in modules inside the closure
(`ThreeRowC1..C4Boundary`, `ThreeRowC4InteriorN4`, `QuantumInteger`). **The repair covers all
six graded nodes.** This is a property of today's import list, not a guarantee — see follow-up.

## Status at end of session

| run | branch | expectation | result |
|---|---|---|---|
| `33550684997` | `main` (`a7225b0`) | green, audit runs | in progress at session end |
| `33550907849` | `main` (`5fea80f`) | green, audit runs | in progress at session end |
| `33551104884` | `ci-negative-control` | **RED** — audit rejects `sorryAx` | in progress at session end |

Runs on this repo take 22–46 min; the session did not outlast them. **The red observation is
therefore PENDING, and until it is made the repair is unverified.** Stating that plainly is the
whole point of this exercise: the previous session's error was not a wrong fix, it was citing a
green run as evidence.

**Next session, first action:** read run `33551104884`.
- **Red, with a log line naming `sorryAx` in `ci_negative_control`** → repair verified; the
  badge is a Lean detector; only then delete the remote branch `ci-negative-control` (cleanup
  owed since 08-31 — it is the evidence, so not before).
- **Green** → the repair does not work. Withdraw the claim per Job 1(b): say in the README and
  in every registry node mentioning CI that the badge is a build check, not a verification
  check, and must never be cited as evidence for `lean-verified`.
- **Red for a non-audit reason** (axiom-audit clone/build failure, toolchain, runner) → that is
  a *permanently red badge*, as useless as a permanently green one. Fall back to the log-grep
  guard (Job 1 option 2).

## Follow-up owed

1. **Read run `33551104884`.** Nothing else in this list matters until that is done.
2. **The audit's coverage rides on a hand-maintained import list.** `lean-action` exposes no
   pass-through for `--modules-from`, so there is no way to ask it to audit every module under
   the library directory. Until there is, a new module that the root forgets to import is
   built but **not audited**, and nothing reports that. Candidate guard: `mk_all-check: true`
   (lean-action input, `lake exe mk_all --check`, "check all files are imported") — unverified
   on this project, must itself be negative-controlled before being trusted.
3. **`nanoda` remains off** — see above. If enabled, measure its runtime first.
4. **Q63 trust-boundary violation** in `fock-ribbon-sign-operator.json`, pre-existing, needs a
   mathematical judgement, not a Lean session.

## Axioms

No new declaration was formalised this session, by instruction. The one declaration touched was
renamed only, and re-checked locally after the rename:

```
'TworowD4Kernel.QuantumInteger.C4_coefficient_identity' depends on axioms:
  [propext, Classical.choice, Quot.sound]
```

## Commits

- `a7225b0` — CI: turn on `axiom-audit`
- `5fea80f` — bring `QuantumInteger` under the single audited root
- `ci-negative-control` — planted `sorry`, imported from the root; **do not merge**, and do not
  delete until the red run is observed
