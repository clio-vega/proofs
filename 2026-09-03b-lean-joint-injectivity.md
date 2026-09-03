# LEAN 2026-09-03 (cycle 2) — joint injectivity of `uglovLabel`, and the attained bounds

**Project:** `clio-vega/tworow-d4-kernel`, module `TworowD4Kernel/WindowBridge.lean`
**Commit:** `be6f0ae` on `main` (parent `5037b78`)
**Paper proof:** `proofs/2026-09-01-Q63-wedge-fock-normalisation.tex`, `cor:range` / `thm:closed` / §`sec:conv`
**External source:** Uglov, `arXiv:math/9905196` §3, eq. `eq:kcd`
**Predecessor note:** `2026-09-03-lean-window-bridge.md`

## Result

`lake build` → **Build completed successfully (2973 jobs)**, exit 0. **Zero sorries** — the
only occurrences of the string in `TworowD4Kernel/` are the words "sorry-free" inside four
docstrings; no declaration uses `sorry`. Both jobs in the brief landed. Job 3's premise
turned out to be false; see below.

## Job 1 — joint injectivity (the module's own overclaim)

**Target, as set:**

```lean
theorem uglovLabel_injective :
    Function.Injective (fun p : Fin n × Fin l => uglovLabel p.1 p.2)
```

This is now proved. It is the honest content of "the label *is* the residue", and it was
**not** implied by the two axis-wise lemmas the module already had.

**Mathlib was searched first, and the answer is "no, and here is precisely why".**
`finProdFinEquiv : Fin m × Fin n ≃ Fin (m * n)` (`Mathlib/Logic/Equiv/Fin/Basic.lean`) is
the right object but is not directly reusable, for two reasons now written into the module
docstring rather than left to be rediscovered:

1. it sends `(x, y) ↦ y + n * x`, so matching `uglovLabel` needs the factors in the order
   `ℓ * n` and then a transport along `Nat.mul_comm`;
2. its codomain is `Fin (nℓ)`, which is not `ZMod (nℓ)` without a `NeZero` hypothesis that
   `uglovLabel_injective` does not otherwise need.

Its `left_inv` field is proved by `Nat.add_mul_mod_self_left` and `Nat.add_mul_div_left`.
The private helper `add_mul_inj_of_lt` uses **those two core lemmas directly** — the same
mathematics, without the packaging. Five lines. I take this to be the short proof the brief
asked for; wrapping the equivalence would have been longer than the thing it wraps.

**The trap was real and is discharged.** `n` and `ℓ` need not be coprime, so this is not
`ZMod.chineseRemainder`; it is base-`n` positional notation, and `c < n` is used
essentially, in `Nat.mod_eq_of_lt hc` and `Nat.div_eq_of_lt hc`. I did not leave that as an
assertion — see the third witness below.

### Non-vacuity, and the control (all three are named theorems, not `example`s, so the CI axiom audit counts them)

| Declaration | What it measures |
|---|---|
| `uglovLabel_injective_nonvacuous` | `n = 2` **and** `ℓ = 3`; the pairs `(0,1)` and `(1,2)` differ in **both** coordinates and get distinct labels. Neither axis-wise lemma can see this pair. `decide`. |
| `uglovLabel_not_injective_mod_n` | **The control.** The same formula read mod `n` instead of mod `n·ℓ` is **not** injective — the runner is invisible mod `n`. `decide`. A proof of joint injectivity that would also typecheck for this map has an unmeasured kernel. |
| `add_mul_inj_of_lt_needs_lt` | The hypothesis `c < n` is load-bearing: drop it and the arithmetic core is **false** — `0 + 2·1 = 2 + 2·0` with `0 ≠ 2`. `decide`. |

The variation was named before the check, per `degenerate-evidence-has-a-kernel`: the
direction tested is *the modulus*, and the control moves exactly that and nothing else.

## Job 2 — "exactly once", which makes the bounds attained

Job 1 landed, so Job 2 was in scope. Six new declarations:

- `card_filter_window_labels` — the integers of an `N`-window whose residue lies in a
  prescribed `L` number **exactly** `L.card`, not merely at most. Proved by
  `Finset.card_bij` against `card_window_eq_one` (`WindowLemma`). This is the surjective
  half the upper-bound argument never needed.
- `uglovLabel_bijective`, `existsUnique_label` — the label map is a bijection onto
  `ZMod (nℓ)`; every bead has **exactly one** label. This is where joint injectivity is
  load-bearing: axis-wise injectivity gives no uniqueness statement here at all.
- `existsUnique_window_of_label` — §`sec:conv`'s "exactly once" remark, in full.
- `exists_card_eq_sub_one_of_runner_ne`, `exists_card_eq_sub_one_of_residue_ne` — **both**
  halves of `cor:range` upgraded from bounds to **attained** bounds, in *every* window and
  for all `n, ℓ` with `NeZero (n·ℓ)`, not just at the `n = 2, ℓ = 3` witness that pinned
  the constant yesterday.

`NeZero (n·ℓ)` is a real hypothesis, not decoration: at `N = 0` the ambient `ZMod 0 = ℤ` is
infinite and the counting statements are false.

**Honest limit, unchanged, and it must not drift:** this is sharpness of the **Lean**
statement over an arbitrary `Finset ℤ` meeting the hypotheses. It is **not** the paper's
attainment claim, which is about a genuine bead configuration and therefore needs the
abacus.

## `#print axioms`, taken through the root import and compared as sets

No new module was created, so the root import closure is **unchanged at 20 modules** (19
`import` lines in `TworowD4Kernel.lean` plus the root itself). `WindowBridge` entered that
closure on 2026-09-02, so the new declarations are in the audit's view by construction —
and every line below was nonetheless produced from a file whose only import is
`TworowD4Kernel`, the root, not the module.

Compared **as sets** against `{propext, Classical.choice, Quot.sound}` by script, not by
eye — `native_decide` emits `<decl>._native...` and never `Lean.ofReduceBool`, so a
blocklist would pass it silently.

```
<: TworowD4Kernel.uglovLabel_injective:                    [Quot.sound, propext]
<: TworowD4Kernel.uglovLabel_injective_nonvacuous:         [Quot.sound, propext]
<: TworowD4Kernel.uglovLabel_not_injective_mod_n:          [Quot.sound, propext]
<: TworowD4Kernel.add_mul_inj_of_lt_needs_lt:              []
== TworowD4Kernel.card_filter_window_labels:               [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.uglovLabel_bijective:                    [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.existsUnique_label:                      [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.existsUnique_window_of_label:            [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.exists_card_eq_sub_one_of_runner_ne:     [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.exists_card_eq_sub_one_of_residue_ne:    [Classical.choice, Quot.sound, propext]
== TworowD4Kernel.card_le_sub_one_of_runner_ne:            [Classical.choice, Quot.sound, propext]   (regression, unchanged)
== TworowD4Kernel.card_le_sub_one_of_residue_ne:           [Classical.choice, Quot.sound, propext]   (regression, unchanged)

VERDICT: ALL WITHIN STANDARD THREE
```

`==` is set equality, `<:` a strict subset. No declaration emitted a `._native` suffix:
every witness is closed by `decide`, none by `native_decide`.

## Job 3 — the badge. **The brief's premise was false, and the deliverable already existed.**

The brief asked me to set `nanoda-allow-sorry: false` and then plant an error and watch it
fire. Checking the workflow before acting:

1. **`nanoda-allow-sorry` does not appear in `.github/workflows/lean_action_ci.yml`** except
   inside a comment. It has not been in the configuration since 2026-09-01. The string in
   the run log is **lean-action echoing its own default inputs**, not this repo's
   configuration — which is exactly why the diagnosis was seductive, and why the caveat
   survived five days in the registry.
2. `axiom-audit: true` with `axiom-audit-allow: "propext,Classical.choice,Quot.sound"` and
   `axiom-audit-root: "TworowD4Kernel"` has been on since 2026-09-01.
3. **The red run already exists**, and it is red *because of the audit*. Run
   **`33551104884`** (branch `ci-negative-control`, head `f53fb50`, 2026-09-01),
   `completed / failure`:

   ```
   axiom-audit: 1 declaration(s) under 'TworowD4Kernel' use disallowed axioms:
     TworowD4Kernel.ci_negative_control → [sorryAx]
   allowed: [propext, Classical.choice, Quot.sound]
   ##[error]axiom-audit check failed
   ```

   A planted `sorry` **does** turn this repo red. The detector has been shown to fail. That
   is the thing that makes a green badge mean something, and it was already done.

4. **Why `test` / `lint` are still skipped**, read off the same log rather than guessed:

   ```
   lake check-test failed -> will not run lake test
   lake check-lint failed -> will not run lake lint
   ```

   lean-action runs `lake check-test` / `lake check-lint` first and skips when they fail;
   they fail because `lakefile.toml` declares no `test_driver` and no `lint_driver`. That
   is a **lakefile** question, not a workflow flag, and I deliberately did not open it in a
   formalisation session. Recorded, not fixed.

So Job 3 required no change to the workflow. What it required was **reading the file before
editing it** — `recorded-facts-calcify`, fired on a fact I had carried for five days in my
own registry. The caveat is struck in the node and replaced with the run ID above.

### In flight at session end

Run **`33728508057`** (branch `main`, head `be6f0ae`), queued 2026-09-03T07:31Z, ~45 min on
this repo. Close it with:

```
gh run view 33728508057 --repo clio-vega/tworow-d4-kernel
gh run view 33728508057 --repo clio-vega/tworow-d4-kernel --log | grep axiom-audit
```

The audited-declaration count should exceed run `33688549430`'s **549** by the new
declarations. **The `lean-verified` grade does not depend on this run**: it rests on the
local `lake build`, the set-compared `#print axioms` above, and run `33551104884` for the
evidence that the detector can go red.

## Sorry count

**Zero.** Nothing is bookmarked; there is nothing outstanding inside the module.

## Next target — the abacus, and it is named, not started

The `κ`-coordinate ↔ `k`-coordinate change is now the module's **only** stated gap. The
paper's `a` is `#{d' > d : κ ∈ M_{d'}} + #{d' < d : κ + e ∈ M_{d'}}`; its equality with the
`k`-coordinate window count is the last paragraph of the proof of `thm:closed`, solving
`eq:kcd` to get `μ' = μ` for `d' > d` and `μ' = μ − 1` for `d' < d`.

Per the brief I wrote **no code** for this. Design paragraph only:

> The type wants to be a bead set per runner — `M : Fin ℓ → Set ℤ` (or a `Finset ℤ` with a
> boundedness hypothesis), together with the charge datum that fixes the window. The
> solving step is not an induction and should not become one: it is the statement that for
> fixed `(c, d)` and `κ`, the map `d' ↦ (the unique k ≡ uglovLabel c d' in the window
> anchored at κ)` has `μ' = μ` exactly when `d' > d`. That existence-and-uniqueness is
> **already available** — it is `existsUnique_window_of_label`, proved today. So the abacus
> lemma should be statable as a `Function.Bijective` between the paper's `κ`-side index set
> and the `k`-side window fibre, with today's `existsUnique_window_of_label` as its engine,
> and no new arithmetic. That is the shape to try first next session.

This is a multi-session piece of work. Starting it badly today would have been worse than
not starting it.
