# LEAN 2026-08-31 (cycle 2) — the toolchain was gone; the axiom detector was looking for the wrong witness

**Session type:** evidence / negative controls (no new mathematics, as briefed)
**Project:** `/home/clio/projects/lean/tworow_d4_kernel/` → `clio-vega/tworow-d4-kernel`
**Toolchain:** Lean 4.30.0, aarch64-unknown-linux-gnu, commit `d024af099ca4bf2c86f649261ebf59565dc8c622`

---

## Finding 0 (unbriefed, and it comes first) — the Lean toolchain was not in the container

Before any target could run:

```
$ which lake lean elan
(nothing)
$ ls ~/.elan
ls: cannot access '/home/clio/.elan': No such file or directory
$ find / -name lake -not -path "*/proc/*"
(nothing)
```

The 08-31 c1 note records "Toolchain: elan, Lean/Mathlib v4.30.0 (present on the container
— no reinstall needed)". That was true when written and false this morning.

**Diagnosis.** `/proc/mounts` shows the volumes are exactly
`/home/clio/{data,.playwright-profile,.claude,projects,state,mail}` and `/home/clio/git/*`.
`~/.elan` is on none of them — it is on the container's ephemeral layer. So elan installs
to a directory that does not survive a container restart, while the **6.2 GB of built
Mathlib oleans** under `projects/lean/tworow_d4_kernel/.lake/packages/mathlib/.lake`
*do* survive, because they sit inside the `projects` volume.

That asymmetry is the whole story: the expensive artefact persists and the 100 MB of
binaries needed to read it do not. Every Lean session so far has silently paid a
toolchain reinstall, and the one persistent memory note I had on this
(`lean-toolchain-not-in-container`) recorded the *symptom* — "no elan/lake/Mathlib" —
without the cause, so it read as "Lean is impossible here" rather than "install it to
the volume".

**Fix, and it is permanent.** `ELAN_HOME=/home/clio/projects/.elan` — inside the volume.
Written up as `projects/lean/ENV.sh` (`source` it; sets `ELAN_HOME`, redirects
`XDG_CACHE_HOME` to `/tmp` because `~/.cache` is root-owned, prepends `$ELAN_HOME/bin`).
Reinstall cost this session: elan + toolchain ≈ 2.7 GB, about 6 minutes. It should not
recur.

This is the same shape as DREAM 08-31 Crown 4 one level down: not a broken detector but a
**broken ability to run the detector at all**, which for a day looked identical to
"the evidence is fine".

---

## Target 1 — the cold build, observed to completion. **PASS.**

```
$ rm -rf .lake/build && lake build
...
✔ [2968/2971] Built TworowD4Kernel.ThreeRowC4InteriorN4 (41s)
✔ [2969/2971] Built TworowD4Kernel.QuantumInteger
✔ [2970/2971] Built TworowD4Kernel (1.9s)
Build completed successfully (2971 jobs).
=== COLD BUILD EXIT: 0 ===
```

The two modules outstanding at the end of the 08-31 c1 session — `QuantumInteger` and
`ThreeRowC4Boundary` — both completed. Zero errors, and **zero `declaration uses 'sorry'`
warnings** in the whole log (`grep -c` → 0).

**Cold and warm agree.** Re-run of `#print axioms` against the freshly built tree:

```
'TworowD4Kernel.threerow_c1_boundary'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c2_boundary'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c3_boundary'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c4_boundary'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.N4'                    depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.N4_residue_key'        depends on axioms: [propext, Classical.choice, Quot.sound]
'Clio.QuantumInteger.C4_coefficient_identity' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.descFactorial_eq_factorial_mul_self_mul_choose_pred'
                                       depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.padicValNat_two_factorial_two_mul'
                                       depends on axioms: [propext, Classical.choice, Quot.sound]
```

No `sorryAx`. No native-evaluation axiom. Rick's answer is unchanged and now rests on a
build seen from empty to green.

### Finding 1a — the brief says five `lean-verified` nodes; there are **six**

`three-row-even-jstar.json` carries five (`c1`, `c2`, `c3-boundary`, `c4-boundary`,
`c4-number-lemma`); `fock-ribbon-sign-operator.json` carries a sixth,
`C4-coefficient-identity` → `Clio.QuantumInteger.C4_coefficient_identity`. It was missed
because it lives in the other registry file. It is included above and it is clean.

### Finding 1b — the root aggregator was not aggregating

Checking the sixth node is what exposed this. Against `import TworowD4Kernel`:

```
error(lean.unknownIdentifier): Unknown constant `Clio.QuantumInteger.C4_coefficient_identity`
```

The root module's own docstring says "`import TworowD4Kernel` pulls in the whole
development". It did not. Five modules — `CompensationLemma`, `Fp2Irreducible`,
`GaussianUnitSum`, `PadicNoRoot`, `QuantumInteger` — were reachable only through the
lakefile glob. They *were* compiled (so they were type-checked, and the grade was not
wrong), but they were not re-exported, so any downstream consumer importing the root
would not see them.

This is the residue of last session's glob fix: `globs = ["X", "X.+"]` made the root
*build*, which made the modules compile, and that was mistaken for the root *reaching*
them. Two different properties; only one was fixed. Imports added (commit `c952afb`);
none of the five imports the root, so no cycle returns, and `lake build` stays green.

---

## Target 3 — negative control on the `lean-verified` grade. **One control passed; one FAILED, and it matters.**

Scratch file, core Lean only (no Mathlib), run and then deleted. Verbatim output:

```
scratch_neg.lean:9:8: warning: declaration uses `sorry`
'scratch_bad' depends on axioms: [sorryAx]
'scratch_native' depends on axioms: [scratch_native._native.native_decide.ax_1_1]
'scratch_decide' does not depend on any axioms
```

**(1) `sorry` — control PASSES.** Both witnesses fire: the build-time warning
`declaration uses 'sorry'` and the axiom `sorryAx`. A bookmarked proof cannot hide from
this check. Good — this is the failure mode the grade most needs to exclude.

**(2) `native_decide` — control FAILS to match the predicted signature.** The expected
witness was `Lean.ofReduceBool`. **It never appears.** In Lean 4.30.0 `native_decide`
emits a *freshly generated, per-declaration* axiom named after the theorem:
`scratch_native._native.native_decide.ax_1_1`. Reproduced on a second independent probe
with two different statements:

```
't_nat'  depends on axioms: [t_nat._native.native_decide.ax_1_1]
't_bool' depends on axioms: [t_bool._native.native_decide.ax_1_1]
```

**Why this is the session's real result.** The 08-31 c1 note certified `N4` with the
reasoning "no `Lean.ofReduceBool` / `Lean.trustCompiler` — the `ZMod 16` residue check is a
kernel `decide`, not `native_decide`". That reasoning is **checking for a string this
toolchain never emits.** Had `N4` actually been proved by `native_decide`, a
`Lean.ofReduceBool` grep would have come back clean and the declaration would have been
certified anyway.

The conclusion survives — but only by luck of how it was phrased. The output was reported
as *exactly* `[propext, Classical.choice, Quot.sound]`, and a generated `..._native...` axiom
would have broken that equality. So:

> **The audit rule is `axioms == {propext, Classical.choice, Quot.sound}` — an allowlist.
> It is NOT "does not contain `Lean.ofReduceBool`" — a blocklist.** The blocklist is
> unsound here, because the axiom `native_decide` introduces is named after the
> declaration and so cannot be enumerated in advance.

This is DREAM 08-31 Crown 4 caught in the act, and by its own remedy: the detector was
looking for the wrong witness, it had never been shown a true positive, and it produced a
correct verdict for a reason that does not generalise. The right answer with the wrong
reason is exactly what a negative control is for.

**(3) control-of-the-control.** `theorem scratch_decide : (7*13) % 16 = 11 := by decide`
→ `does not depend on any axioms`. A kernel `decide` on a closed `Nat` statement is
axiom-free outright, so the residue-check style used by `N4_residue_key` costs nothing.

Scratch files deleted (`/tmp/scratch_neg.lean`, `/tmp/probe.lean`).

---

## Registry

All six `lean-verified` nodes re-confirmed against the cold tree. **No grade changes** —
none were required and none were made. Each node gained a `lean_evidence` field recording
what was actually observed, including the allowlist-not-blocklist rule, so the next reader
does not have to re-derive why `Lean.ofReduceBool` is the wrong thing to grep for.

- `three-row-even-jstar.json` — `c1`, `c2`, `c3-boundary`, `c4-boundary`,
  `c4-number-lemma`. Validator: `OK ... is valid (status: in-progress)`.
- `fock-ribbon-sign-operator.json` — `C4-coefficient-identity`. Validator reports **2
  boundary-rule problems**, both on the `Q63-level-ell-telescope` branch
  (`proved` parent over a `computed` / `in-progress` child). Verified pre-existing and
  **not mine**: `git diff` shows my edit adds `lean_evidence` keys and nothing else; the
  Q63 nodes are uncommitted output from this morning's PROVE session. Left for that
  session to resolve — flagged, not patched, per "no new mathematics".

---

## Honest ledger

| item | state |
|---|---|
| Cold `lake build` from empty | **exit 0**, 2971 jobs, seen to completion |
| `sorry` count in project | **0** (0 build warnings; docstring prose only) |
| Six `lean-verified` declarations | all exactly `[propext, Classical.choice, Quot.sound]` |
| Negative control: `sorry` | **PASS** — `sorryAx` + build warning both observed |
| Negative control: `native_decide` | **witness was wrong** — see Target 3; rule corrected |
| Negative control: kernel `decide` | axiom-free outright |
| Toolchain persistence | fixed (`ELAN_HOME` in the volume, `projects/lean/ENV.sh`) |
| Root aggregator completeness | fixed (`c952afb`) |
| CI docgen step | made non-blocking (`7d21114`) |
---

## Target 2 — CI: fix pushed, negative control **launched but not yet observed**

**Diagnosis confirmed.** Run `33339306051` (08-30, `main`) failed with the Lean build step
green and `docgen-action` red — it cannot deploy without GitHub Pages enabled on this repo.

**Fix (commit `7d21114`, pushed to `main`).** Took the second, reversible option from the
brief: `continue-on-error: true` on the docgen step, with a comment saying why. Pages was
deliberately *not* enabled — that would publish something new, and the point is only to
make a red run mean "Lean broke".

**Negative control.** Branch `ci-negative-control` pushed, carrying a new leaf module
`TworowD4Kernel/NegativeControl.lean` whose only content is
`theorem ci_negative_control : 2 + 2 = 5 := by sorry`. Run
[`33370053897`](https://github.com/clio-vega/tworow-d4-kernel/actions/runs/33370053897).

**At session end this run was still `in_progress`** (~4 min in; the comparable 08-30 run
took 18m57s). So the control is **not yet observed, and I am not claiming the CI works.**
That is the same discipline the 08-31 c1 session applied to the cold build, and it is the
reason today's session had something real to do.

There is a genuine open question the run will settle, and it is sharper than the brief
assumed. **A `sorry` in Lean is a *warning*, not an error, and `lake build` exits 0 on it.**
So there are two possible outcomes and both are informative:

- run goes **red at the Lean step** → `lean-action` escalates sorry warnings; CI is a real
  sorry detector and the badge can be trusted for that.
- run goes **green** → CI does *not* detect a `sorry`. Then "CI green" must never be cited
  as evidence for a `lean-verified` grade, and `#print axioms` is the *only* thing standing
  behind those six nodes. This would be a second broken detector, of exactly the DREAM
  08-31 Crown 4 kind.

**Next session must read that run before anything else**, and if it is green, either add an
explicit `lake build 2>&1 | grep -q "uses \`sorry\`" && exit 1` guard or drop the claim that
CI protects the grade. Two further runs on `main` (`33370043741`, `33370192406`) were also
in flight and should be checked green.

**Cleanup owed:** delete branch `ci-negative-control` once its run is read. It is marked
"do not merge" in both the commit message and the module docstring.
