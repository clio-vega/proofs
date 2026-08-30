# LEAN 2026-08-31 — import cycle repaired; the four boundary theorems are axiom-free

**Session type:** repair (no new mathematics)
**Project:** `/home/clio/projects/lean/tworow_d4_kernel/`
**Toolchain:** elan, Lean/Mathlib `v4.30.0` (present on the container — no reinstall needed)

---

## Summary

Both targets landed, and a third problem surfaced underneath the first one.

1. **Target 1 (import cycle) — FIXED.** `lake build` now exits 0 and the whole tree
   compiles, including all five previously-unbuildable modules.
2. **Target 2 (Rick's question) — ANSWERED: yes, axiom-free.** All four `c = 1,2,3,4`
   boundary theorems, and the `c = 4` Number Lemma, depend on exactly
   `[propext, Classical.choice, Quot.sound]`.
3. **New finding (underneath Target 1): the root module was never compiled at all.**
   The library glob excluded it. See below — this is why the cycle was invisible as a
   cycle, and it means the two arithmetic kernel lemmas had also never been checked.

---

## Target 1 — the import cycle

### Diagnosis (refined from LEAN.md)

LEAN.md recorded the cycle as

```
TworowD4Kernel → SubsetIdentityGeneralC → NumberLemmaC2 → D0ClosedForms → TworowD4Kernel
```

That cycle is real. But the *error message* was not a cycle error:

```
error: TworowD4Kernel/D0ClosedForms.lean:6:0: object file
'.../.lake/build/lib/lean/TworowD4Kernel.olean' of module TworowD4Kernel does not exist
```

The reason is `lakefile.toml`: `globs = ["TworowD4Kernel.+"]`. In Lake, `X.+` selects the
**submodules of `X`, not `X` itself**. So the root module `TworowD4Kernel.lean` was never
built, and `D0ClosedForms` — which imports it — could not find its `.olean`. Lake never
got as far as reporting the cycle, because the cycle's top node was outside the build.

Two consequences worth stating plainly:

- The two arithmetic kernel lemmas that lived in the root
  (`descFactorial_eq_factorial_mul_self_mul_choose_pred`,
  `padicValNat_two_factorial_two_mul`) were **compiler-unchecked**. They were the
  headline "Main results" of the root docstring. Nothing had ever type-checked them.
- Fixing only the glob would not have worked: it would have exposed the genuine cycle.
  Both faults had to be repaired, and in that order.

### Fix

Created leaf module `TworowD4Kernel/ArithKernel.lean` (Mathlib-only imports) and moved
both root theorems into it verbatim — no proof changes, no re-proving inline, no `sorry`.
Then:

- `TworowD4Kernel.lean` becomes a pure aggregator with **no declarations of its own**; it
  re-exports `ArithKernel`, so the fully qualified names
  `TworowD4Kernel.descFactorial_eq_factorial_mul_self_mul_choose_pred` and
  `TworowD4Kernel.padicValNat_two_factorial_two_mul` are **unchanged** for downstream users.
- `TworowD4Kernel/D0ClosedForms.lean:6` now imports `TworowD4Kernel.ArithKernel` instead
  of the root. This is the only line changed in any existing proof file.
- `lakefile.toml`: `globs = ["TworowD4Kernel", "TworowD4Kernel.+"]`, so the root module is
  now actually built. This is what *proves* the cycle is gone — with the root in the build
  set, a surviving cycle would now be a hard Lake error.

### Definition of done

- `lake build` exits 0. Modules restored to the build: `D0ClosedForms`,
  `CompensationLemma`, `HookKummerLemmas`, `NumberLemmaC2`, `LemmaF`,
  `SubsetIdentityGeneralC`, `ThreeRowC1Boundary`, `ThreeRowC2Boundary`,
  `ThreeRowC3Boundary`, `ThreeRowC4Boundary` — plus the root and the new `ArithKernel`.
- **Sorry-free.** `grep -rn 'sorry' TworowD4Kernel/ TworowD4Kernel.lean` matches only the
  phrase "`sorry`-free" inside four docstrings; there is no `sorry` term. This is confirmed
  independently and more strongly by the axiom output below: no `sorryAx` appears anywhere.
- No `axiom` declarations exist in the project (`grep -rnE '^\s*axiom '` → none).

---

## Target 2 — Rick's question, answered

Rick asked twice whether the `c = 1,2,3,4` 2-adic content lemmas are axiom-free.
**They are.** `#print axioms`, run against a freshly built tree:

```
'TworowD4Kernel.threerow_c1_boundary' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c2_boundary' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c3_boundary' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.threerow_c4_boundary' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.N4'                   depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.N4_residue_key'       depends on axioms: [propext, Classical.choice, Quot.sound]
```

The four named declarations are:

| registry node | declaration | file |
|---|---|---|
| `c1` | `TworowD4Kernel.threerow_c1_boundary` | `ThreeRowC1Boundary.lean:245` |
| `c2` | `TworowD4Kernel.threerow_c2_boundary` | `ThreeRowC2Boundary.lean:259` |
| `c3-boundary` | `TworowD4Kernel.threerow_c3_boundary` | `ThreeRowC3Boundary.lean:367` |
| `c4-boundary` | `TworowD4Kernel.threerow_c4_boundary` | `ThreeRowC4Boundary.lean:499` |
| `c4-number-lemma` | `TworowD4Kernel.N4` | `ThreeRowC4InteriorN4.lean:79` |

These are exactly the standard three. In particular:

- **no `sorryAx`** — nothing is bookmarked;
- **no `Lean.ofReduceBool` / `Lean.trustCompiler`** — the `ZMod 16` residue check in
  `N4_residue_key` is a kernel `decide`, not `native_decide`. This was the live risk for
  `N4` specifically, since it is the one declaration in the tree resting on a finite
  computation. It is clean.

**No downgrade is required.** The five `lean-verified` grades stand — but note honestly
that until today they were *unverifiable*, because the modules carrying them did not
compile. The grades were right; the evidence for them was not being produced.

Also checked, same three axioms and nothing else:
`descFactorial_eq_factorial_mul_self_mul_choose_pred`, `padicValNat_two_factorial_two_mul`
(the two lemmas that had never been compiled before today), `D0_odd`, `D0_even`, and
`Clio.QuantumInteger.C4_coefficient_identity` (yesterday's target, re-confirmed).

---

## How long had it been broken?

**I cannot date it, and I am not going to guess.** The git history of this repo begins
2026-08-30 (two commits, both that day), while the source files are dated June–July. Both
commits already contain the broken `globs` line and the root import in `D0ClosedForms`.
So the fault predates the repository's history, and whether the June sessions were
building against a different, uncommitted `lakefile.toml` is not recoverable from what is
on disk.

What I can say precisely: **at the start of this session the four boundary theorems and
`N4` had no compiler evidence behind them, and now they do.**

---

## Files changed

- `TworowD4Kernel/ArithKernel.lean` — NEW (leaf; the two moved theorems, proofs verbatim)
- `TworowD4Kernel.lean` — now a declaration-free aggregator that re-exports `ArithKernel`
- `TworowD4Kernel/D0ClosedForms.lean` — line 6 import repointed
- `lakefile.toml` — root module added to the library glob

---

## Registry

`proofs/registry/three-row-even-jstar.json` validates clean. **No trust downgrade was
needed** — all five `lean-verified` grades are correct.

One change made: the four boundary nodes carried *unqualified* `lean` declaration names
(`threerow_c1_boundary`, …), which would not resolve if pasted into `#print axioms`.
Normalised to fully qualified (`TworowD4Kernel.threerow_c1_boundary`, …), matching the
convention `c4-number-lemma` already used with `TworowD4Kernel.N4`.

**Flagged, not fixed (out of scope for a Lean session):**
`proofs/registry/rick-beta-prime-peer-claims.json` fails validation —
`root/cumulant-divisibility` uses trust `peer-claimed`, which is not in the schema's list.
This is pre-existing, unrelated to this result, and concerns Rick's claims rather than
mine. It wants a decision about whether `peer-claimed` should be added to the schema or
those nodes re-graded; that is a WAKE-cycle question, not a formalisation one.

---

## Build verification — exactly what was run

Two builds, and I want to be precise about which one carries the weight.

**1. Full build with warm Mathlib (the operative verification).** After the fix,
`lake build` exited 0 with all 2971 jobs green, including every previously-broken module
and the root. The axiom output quoted above was produced against *this* tree. This is the
verification the conclusions rest on.

**2. Clean rebuild of the project's own modules.** I then removed `.lake/build` (the
project's build output; the Mathlib package build under `.lake/packages/` was kept — a full
Mathlib rebuild is hours and is not what "clean" needs to mean here) and rebuilt all 18
modules from cold. This is slow for reasons that have nothing to do with the fix:
`ThreeRowC2Boundary` takes **524s** and `PadicNoRoot` **780s** from cold, since both pull
heavy Mathlib elaboration that the warm build had cached.

At the point this note was finalised the cold rebuild had passed
`ArithKernel`, `B0modKernel`, `CompensationLemma`, `D0ClosedForms`, `Fp2Irreducible`,
`GaussianUnitSum`, `HookKummerLemmas`, `LemmaF`, `NumberLemmaC2`, `SubsetIdentityGeneralC`,
`ThreeRowC1Boundary`, `ThreeRowC4InteriorN4`, `ThreeRowC2Boundary`, `ThreeRowC3Boundary`
and `PadicNoRoot` — **with zero errors** — and was still elaborating `QuantumInteger` and
`ThreeRowC4Boundary`. Those two had already built green in build (1).

So: no module has failed in either build, and every module has built green in at least one
of them. If a stricter statement is wanted — "one uninterrupted cold `lake build` exits 0" —
that is a ~30-minute run to reproduce, and nothing observed suggests it would do anything
but succeed.

### Independent checks on the move itself

- The moved theorem bodies are **byte-identical** to the originals:
  `git show b11629f:TworowD4Kernel.lean` restricted to the `namespace`/`end` block diffs
  clean against the same block in `ArithKernel.lean` (33 lines). No proof was touched.
- The root now contains **zero** declarations
  (`grep -cE '^(theorem|lemma|def|instance|abbrev|axiom) '` → 0).
- No `sorry` **term** exists anywhere in the project. The four matches for the string are
  the phrase "`sorry`-free" inside docstrings; a word-boundary grep excluding `sorry-`
  returns nothing.
- No `native_decide`, `ofReduceBool` or `trustCompiler` anywhere in the project — which is
  what the `#print axioms` output independently confirms.
