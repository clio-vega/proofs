# LEAN 2026-10-22 — the ε-profile is a Möbius transform of `1_M`

**Status: sorry-free.** 17 named declarations, `#print axioms` = the standard three on all
17, checked against a planted-sorry positive control.

- **Project:** `projects/lean/tworow_d4_kernel` (Lean 4.30.0, Mathlib `v4.30.0`)
- **Module:** `TworowD4Kernel/ProfileMoebius.lean` — commit `46831c7` in `clio-vega/tworow-d4-kernel`
- **Paper proof:** `proofs/2026-10-22-boundary-valuation-two-term.tex`, `thm:profile`
  (+ `cor:q546no`, `cor:q546yes`), PROVE 2026-10-22, answering **Q546**
- **Registry:** `proofs/registry/motegi-typeBC-multiplication.json`, new child node
  `eps-profile-moebius-inversion-arithmetic` at `trust: lean-verified`
- **Pre-registration:** `proofs/2026-10-22-lean-PREREG.md`, written before the first build

## Target substitution — declared

The brief named the *corollary* ("`f_M` determines `M`") and said: if the prove slot lands
Q546 first, formalise **its** statement instead and say so. It did, this morning. So the
target here is the prove slot's `thm:profile` — the **explicit inversion**
`1_M(n) = Σ_{k≥1} μ(k) f_M(kn)` — which strictly contains the brief's statement
(`eq_of_profile_eq` is derived from it below, not assumed).

The brief also predicted Möbius would resist in Lean and suggested downward induction on
`max M` as the likely-easier route. **That prediction was wrong**, and interestingly so:
Möbius went through in one sitting, because the paper's reindexing
`Σ_k μ(k) Σ_{kn∣n'} G = Σ_{n∣n'} G · Σ_{k∣n'/n} μ(k)` is *literally* `Finset.sum_comm`
followed by one divisor-sum identity. I never attempted the induction route, so I am
grading it **untested**, not harder.

## What Mathlib does and does not have

Mathlib has Möbius inversion for the sum-over-**divisors** transform only
(`ArithmeticFunction.moebius_mul_coe_zeta`, `coe_mul_zeta_apply`). The paper's transform is
a sum over **multiples**. `moebius_inversion_over_multiples` is the transposed statement,
proved *from* the divisor form via `sum_divisors_moebius`. This was the session's actual
work; everything else is bookkeeping on top of it.

## The declarations

### Engine
```lean
theorem moebius_inversion_over_multiples (B : ℕ) (G : ℕ → ℤ) {n : ℕ}
    (hn1 : 1 ≤ n) (hnB : n ≤ B) :
    ∑ k ∈ Icc 1 B, (μ k) * (∑ a ∈ (Icc 1 B).filter (fun a => k * n ∣ a), G a) = G n
```
Stated for an arbitrary `G : ℕ → ℤ` over the window `[1,B]`, with **no support
hypothesis** — the window *is* the support. That framing is what removed the finiteness
bookkeeping the paper carries in prose ("the sum being finite because `M` is finite").

### Paper displays
```lean
def profile (M : Finset ℕ) (d : ℕ) : ℕ := (M.filter (fun a => d ∣ a)).card

theorem indicator_eq_sum_moebius_profile {B : ℕ} {M : Finset ℕ} (hM : M ⊆ Icc 1 B)
    {n : ℕ} (hn1 : 1 ≤ n) (hnB : n ≤ B) :
    ∑ k ∈ Icc 1 B, (μ k) * (profile M (k * n) : ℤ) = if n ∈ M then (1 : ℤ) else 0

theorem eq_of_profile_eq {M M' : Finset ℕ} (hM : ∀ a ∈ M, 1 ≤ a) (hM' : ∀ a ∈ M', 1 ≤ a)
    (h : ∀ d, 1 ≤ d → profile M d = profile M' d) : M = M'
```
`eq_of_profile_eq` is the brief's named target, now a corollary.

### Application layer
```lean
lemma profile_typeB_pair  : profile {1} 1 = 1 ∧ profile {1} 2 = 0        -- (1,0)
lemma profile_typeC_pair  : profile {2} 1 = 1 ∧ profile {2} 2 = 1        -- (1,1)
lemma profile_typeBC_pair : profile {1,2} 1 = 2 ∧ profile {1,2} 2 = 1    -- (2,1)
lemma profile_additive_typeBC (d : ℕ) : profile {1,2} d = profile {1} d + profile {2} d
theorem two_points_suffice :      -- decided over the whole powerset of {1,2}
    ∀ M ∈ ({1,2} : Finset ℕ).powerset, ∀ M' ∈ ({1,2} : Finset ℕ).powerset,
      profile M 1 = profile M' 1 → profile M 2 = profile M' 2 → M = M'
theorem two_points_insufficient : -- {1} vs {3}: share (1,0), separate first at d = 3
    profile {1} 1 = profile {3} 1 ∧ profile {1} 2 = profile {3} 2 ∧
      profile {1} 3 ≠ profile {3} 3
```
`profile_additive_typeBC` is the `f_BC = f_B + f_C` cell that killed Q544's multiplicative
guess: at `d = 2` it reads `1 = 0 + 1 ≠ 2·0`.

## What is NOT formalised — and why the parent node stays at `proved`

The registry node `eps-profile-determines-the-multiplicity-set` claims strictly more than
this file proves. I added a **child** at `lean-verified` and left the parent at `proved`
rather than promoting it, because three things are outside the Lean file:

1. **The root-of-unity dictionary** `ε_d^a = 1 ↔ d ∣ a`. The Lean `profile` is *defined* as
   `#(M ∩ dℤ)`, so the paper's first display is definitional here. (Standard order theory;
   the paper itself discards it in one line. Not a gap — a deliberate boundary.)
2. **`f_M(ε) = 0` for `ε` of infinite order** — the "and nothing more" clause.
3. **The entire root-system half**: rank-blindness `M_m(B_2) = M_m(B_3) = {1}`, and "every
   crystallographic root system has `M_m ⊆ {1,2}`".

(3) is the one that matters for honest grading. `two_points_suffice` is Lean-verified as an
**arithmetic** statement about subsets of `{1,2}`; that root systems *land* in that class is
exactly premise (3), unformalised. So the **"Q546 false as asked"** half is *not*
Lean-verified at all — its falsifier is (3). Promoting the parent would have made the
file's green build look like it had certified the refutation. It has not.

## Lesion register — pre-registered vs actual

Pre-registered in `proofs/2026-10-22-lean-PREREG.md` before the first build. Every verdict
below is settled by a **proof** — removability by a proof of the weakened statement,
irremovability by an exhibited **witness** — never by the absence of a build error.

| # | Hypothesis | Pre-registered | Actual | Settled by |
|---|---|---|---|---|
| L1 | `M, M'` finite | IRREMOVABLE, and **not a binder** (type-level); false over `Set.ncard` by junk value | **CONFIRMED, both halves** | `finiteness_irremovable` — witness: positive integers vs positive evens, every fibre infinite ⇒ `ncard = 0` |
| L2 | `0 ∉ M` | **REMOVABLE** (against the brief) | **CONFIRMED** | `eq_of_profile_eq_zero_allowed` — full proof, `0`'s membership read off at `d = B+1` |
| L3(a) | range `∀ d ≥ 1` ↓ `[1, max (M ∪ M')]` | REMOVABLE/equivalent | **CONFIRMED** | `eq_of_profile_eq_restricted`, via `profile_eq_zero_of_lt` |
| L3(b) | range ↓ `[1, max M]`, one side | IRREMOVABLE | **CONFIRMED** | `range_max_left_insufficient` — witness `M = ∅`, `M' = {5}` |
| L3′ | L2 and L3(a) **jointly** | each survivable, **jointly fatal** | **CONFIRMED** | `range_union_insufficient_when_zero_allowed` — witness `M = {0}`, `M' = ∅` |

### The one worth keeping

**L3′ was the prediction I most wanted to test, and it held.** `0 ∉ M` deletes. The range
restriction to `[1, max (M∪M')]` deletes. Delete *both* and the theorem is **false**:
`M = {0}`, `M' = ∅` has max `0`, so the range `[1,0]` is empty, the hypothesis is vacuous,
and the conclusion fails. A lesion pass that removes one hypothesis at a time — which is
what a lesion pass normally *is* — grades both REMOVABLE and never sees this.

So: **a lesion register over single deletions cannot see a jointly-load-bearing pair.** The
two hypotheses are not independent; `0 ∉ M` is removable *only because* `d` is unbounded,
and the unboundedness is removable *only because* `0 ∉ M`. Each one's removability proof
silently consumes the other. Pre-registering the *interaction* is what caught it — I would
not have computed the pair otherwise.

Note also that L2's confirmation was **against the brief's explicit expectation** ("expect
this one to bite"). The brief was right that `d ∣ 0` puts `0` in every fibre; that is
precisely what makes `0` the one element an unbounded `d` can isolate.

## Axiom check — with positive control

All 17 named declarations:
```
propext, Classical.choice, Quot.sound          (13 declarations)
propext, Quot.sound                            (4: profile_eq_zero_of_lt,
                                                profile_typeB/C/BC_pair, two_points_insufficient)
```
Positive control, same file, same command:
```
theorem planted : (1:ℕ) = 1 := by sorry
'planted' depends on axioms: [sorryAx]
```
The detector fires. No `native_decide` anywhere (all `decide` calls are kernel-reduced).

## Declaration census — environment-derived, SET-compared both ways

Environment scan over the namespace: **21**. Source-side `grep` of declaration keywords:
**19**.

- environment ∖ source = `{profile.eq_1, setProfile.eq_1}` — auto-generated equation
  lemmas for the two `def`s
- source ∖ environment = **∅**

The totals disagree and the record is still clean; had I compared only the totals I would
have had a 21-vs-19 "discrepancy" to chase. No anonymous `example`s in the file, so the
scan is not blind in the 10-19 direction.

## Instrument faults hit this session

1. **`registry_validate.py` takes `--proofs-dir`, not `--files-dir`.** I reached for
   `--files-dir .` — which is *trustcheck's* flag — and argparse exited 2. The failed run
   printed **zero** problem lines, so `grep -c "^  - "` returned **0**: an argparse failure
   is indistinguishable from a clean pass under line-counting. My standing note says "count
   problem lines, never the exit code"; this is the exact case where that rule fails, and
   the repair is to require a *positive* control line (`N problem(s) in …`) before believing
   a zero. Correct invocation: `--proofs-dir .` (1 real problem, now fixed).
   ⇒ *a remembered flag belongs to a script, not a task* — second occurrence (10-20).
2. **A source `grep` for `sorry` gave a false positive** on the prose word "admitted" in a
   docstring. Harmless here, but it is the same shape as every other
   grep-as-detector fault: the token matched text that was not code. The axiom check is the
   instrument; the grep is not.
3. **`lake` not on `PATH`** — handled by `source lean/ENV.sh` per the carried-forward
   correction, no time lost. The `cd` into a subdirectory silently broke later
   `source ENV.sh` calls three times (relative path); cost ~3 wasted invocations.
4. The `unusedVariables` linter flagged nothing. One `flexible tactic` lint on
   `simp at hj`, fixed to `rw [Nat.mul_zero] at hj`.
5. **`git add` on an ignored path short-circuits an `&&` chain, and
   `git rev-parse --short HEAD` then reports success.** `git add -A <files> && git commit`
   — one of the listed files was a `.bak` ignored by design, so `add` exited non-zero, the
   commit never ran, and my follow-up `echo "COMMIT: $(git rev-parse --short HEAD)"`
   printed the **prove slot's** pre-existing hash. I reported a commit that did not happen.
   Caught only by `git show --stat HEAD`, which named someone else's files.
   ⇒ ***a hash read after a failed write names the previous writer.*** Verify a commit by
   its **file list**, never by the hash, which exists whether or not you wrote it. Same
   shape as 10-19's `echo`-after-failed-build (status carries no information) and 10-20's
   clean `git status` — and note the containing repo was *also* the wrong one
   (`.gitignore:6` again, 2nd occurrence: `tworow_d4_kernel/` is ignored by the outer repo
   **and** is its own repo, so the outer `git status` was clean with a new file on disk).
