# Lean: the window bridge — from "each label occurs once in a window" to `a ≤ ℓ − 1`

**Date:** 2026-09-03 · **Session:** LEAN · **Status:** target met, sorry-free

## Target

| | |
|---|---|
| **Project** | `projects/lean/tworow_d4_kernel`, repo `clio-vega/tworow-d4-kernel` |
| **Module** | `TworowD4Kernel/WindowBridge.lean` (new), imported from the root |
| **Commits** | `d5707cf` (mathematics), `5037b78` (non-vacuity guard) on `main` |
| **Paper** | `proofs/2026-09-01-Q63-wedge-fock-normalisation.tex`, §`sec:closed`, Corollary `cor:range` |
| **External source** | Uglov, `arXiv:math/9905196` — the indexation `k = c + n(d−1) − Nμ`, recalled as eq. `eq:kcd` in §`sec:conv` |

The brief named the bridge from *"each label occurs once in the window"* to the paper's
`0 ≤ a ≤ ℓ − 1`. It is done, in the upper-bound direction, for **both** halves of
`cor:range` — the `h ≤ e − 1` half fell out of the same lemma with the two label axes
exchanged, at the cost of one extra injectivity lemma.

## What builds sorry-free

All of it. **Zero sorries; the module ships none.** `lake build` → *Build completed
successfully (2973 jobs)*, exit 0, no linter warning from this module.

| Declaration | Statement |
|---|---|
| `window_injOn` | `Set.InjOn (fun k : ℤ => (k : ZMod N)) ↑(Finset.Ico p (p + N))` |
| `card_le_card_of_window_labels` | `S ⊆ Ico p (p+N)` and `∀ k ∈ S, (k : ZMod N) ∈ L` ⟹ `S.card ≤ L.card` |
| `uglovLabel` | `Fin n → Fin l → ZMod (n*l)`, `(c,d) ↦ c + n·d` |
| `uglovLabel_rep_lt`, `uglovLabel_rep_inj` | the representative is `< n*l`, and equal labels have equal representatives |
| `uglovLabel_injective_runner` | injective in `d`, with `c` fixed |
| `uglovLabel_injective_residue` | injective in `c`, with `d` fixed |
| `card_le_sub_one_of_runner_ne` | **`a ≤ ℓ − 1`** |
| `card_le_sub_one_of_residue_ne` | **`h ≤ e − 1`** |
| *(anonymous `example`)* | non-vacuity witness, closed by `decide` |

Mathlib did most of it, as the brief predicted: `Finset.card_le_card_of_injOn`,
`ZMod.intCast_zmod_eq_zero_iff_dvd`, `Int.eq_zero_of_abs_lt_dvd`, `ZMod.val_cast_of_lt`.
The bridge in this form is not in `v4.30.0` (searched `Data/Finset/Card.lean`, `Data/ZMod/`).

### The bound is not vacuous, and the constant is sharp

An upper bound whose hypotheses are unsatisfiable proves nothing, and one whose constant is
slack is misleading. Both are excluded in the module by a **witness**, not by argument:
`n = 2`, `l = 3`, `N = 6`, `c = 0`, `d = 0`, window `[0, 6)`. The two labels `(0,1)` and
`(0,2)` are the residues `2` and `4`, realised by the beads `2` and `4`, so `S = {2,4}` and
`S.card = 2 = ℓ − 1` — **attained**.

The same witness is the negative control on the constant: it refutes `S.card ≤ ℓ − 2`, so
the statement can see a failure in the direction it is tested. Closed by `decide`, not
`native_decide`; axiom sets unchanged.

Note this is sharpness *of the Lean statement*, over an arbitrary `Finset ℤ` meeting the
hypotheses. It is **not** the paper's attainment claim, which would need a genuine bead
configuration and therefore the abacus — see gap (2) below.

### Which informal step each Lean hypothesis stands in for

The brief asked for this explicitly.

- **`hS : S ⊆ Finset.Ico p (p + N)`** stands in for *"the beads `k_r` with `p < k_r < p+N`"*
  in the proof of `thm:closed`. It is **weaker** than the paper's hypothesis — the paper's
  open window `(p, p+N)` sits inside the half-open `Ico p (p+N)` — so the Lean theorem is
  the stronger one and a caller holding the paper's hypothesis applies it directly.
- **`hlab : ∀ k ∈ S, ∃ d' ≠ d, (k : ZMod (n*l)) = uglovLabel c d'`** stands in for
  *"the beads in the window with `c`-label `c` are those with label `(c,d')`, `d' ≠ d`"*.
  It encodes the consequence of `eq:kcd` that the label **is** the residue mod `N`.
- **`S`** stands in for the bead set itself. The abacus is not formalised; `S` is an
  arbitrary `Finset ℤ`, and the abacus enters only through `hS` and `hlab`. This was the
  brief's instruction — *formalise the interface, not the abacus* — and the interface turned
  out to be two hypotheses wide.

## The gap that remains — what still separates this from `cor:range`

The brief said this gap would not be zero today and that locating it is the deliverable.
It is not zero. It has three components, and they are recorded in the module docstring and
the registry node as well as here.

1. **The coordinate change (the real gap).** The paper *defines*
   `a = #{d' > d : κ ∈ M_{d'}} + #{d' < d : κ+e ∈ M_{d'}}` in the **κ-coordinate**. What
   Lean bounds is the **k-coordinate** window count `a′`. Their equality is the last
   paragraph of the proof of `thm:closed`: solving `eq:kcd` for the occurrence of label
   `(c,d')` in `(p, p+N)` gives `μ′ = μ` when `d′ > d` and `μ′ = μ − 1` when `d′ < d`,
   whence `κ′ = κ` and `κ′ = κ + n` respectively. **That solving step is still on paper.**
   It is a genuinely different piece of work — an arithmetic identification of two indexings,
   not a counting argument — and formalising it means giving the abacus a type.
2. **Surjectivity — the "exactly once" half.** §`sec:conv` asserts each label other than
   `k`'s own occurs in the window *exactly* once. Only **"at most once"** is proved, because
   only that half bounds. The existence half is what would make the bounds **attained**, and
   `cor:range` as stated does not need it — but any later claim that `a` *achieves* `ℓ − 1`
   will. `WindowLemma.card_window_eq_one` already has the exactly-once fact for a *single*
   residue; what is missing is its transport across the label bijection, which needs (1).
3. **The lower bounds `0 ≤ a`, `0 ≤ h`** are trivial for a cardinality and are not stated.

A fourth, smaller thing worth recording because it is a *deliberate* narrowing:
`uglovLabel` is proved injective only **along each axis separately**, never jointly. Joint
injectivity is the honest content of "the label is the residue", and it needs division by
`n`. Each of the two counts varies one coordinate with the other fixed, so axis-wise
injectivity is exactly what the mathematics uses — but the module therefore does *not*
establish that labels and residues correspond bijectively, only that each axis embeds.

Also deliberate: indices are `0`-based (`Fin n`, `Fin l`) where the paper's are `1`-based.
That is the constant shift `−n` of `ZMod N`, a bijection, and every statement here is about
distinctness and counting, so nothing depends on it.

## `#print axioms`

Taken through the **root** import (`import TworowD4Kernel`), i.e. through the same closure
`axiom-audit` follows:

```
'TworowD4Kernel.window_injOn' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.card_le_card_of_window_labels' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.uglovLabel_rep_lt' depends on axioms: [propext, Quot.sound]
'TworowD4Kernel.uglovLabel_rep_inj' depends on axioms: [propext, Quot.sound]
'TworowD4Kernel.uglovLabel_injective_runner' depends on axioms: [propext, Quot.sound]
'TworowD4Kernel.uglovLabel_injective_residue' depends on axioms: [propext, Quot.sound]
'TworowD4Kernel.card_le_sub_one_of_runner_ne' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.card_le_sub_one_of_residue_ne' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Compared **as sets** against `{propext, Classical.choice, Quot.sound}`, not by substring —
`native_decide` emits `<decl>._native.native_decide.ax_N_M` and never `Lean.ofReduceBool`,
so a blocklist would pass it silently. Four are **equal** to the standard three; four are a
**strict subset** (no `Classical.choice`). Nothing outside the allowlist. 8/8 checked.

## Audit coverage

`TworowD4Kernel/WindowBridge.lean` was imported from the root **before any mathematics was
written** (rule 1), because `axiom-audit` follows the *import closure* of `axiom-audit-root`
and the `TworowD4Kernel.+` glob compiles a module without putting it in the audit's view.
Root closure recomputed: **20 modules, was 19**, `WindowBridge` among them. Both new theorem
groups sit inside `namespace TworowD4Kernel`, so they are inside the audited set.

### CI confirmation is PENDING at the time of writing — not confirmed, not assumed

The brief asked me to confirm **from the run log, not the badge**, that `axiom-audit` ran and
that its declaration count went up from 532. **I could not: the run was still in progress when
the session ended.** Runs on this repo take ~45 min and the push landed ~55 min in.

So this is an open item, not a result. Run IDs `33688249633` (commit `d5707cf`) and
`33688549430` (commit `5037b78`), branch `main`. To close it:

```
cd ~/projects/lean/tworow_d4_kernel && source ~/projects/lean/ENV.sh
gh run view 33688549430 --log | grep -iE 'axiom-audit|AXIOM_AUDIT|audited'
```

What must be seen, and what would falsify: `outcome=success` on the `axiom-audit` step — **not
`skipped`**, which is what it silently was before 09-01 and is what the sibling steps (`test`,
`lint`, `nanoda`, `mk_all`) still are — and an `audited N declaration(s)` line with **`N > 532`**.
If `N = 532`, the new module is compiling but is *not* in the audit's view, and the
`lean-verified` grade on `Q63-window-bridge` is not carrying its claimed evidence.

What is confirmed locally, and is the evidence the grade actually rests on today: `lake build`
green (2973 jobs, 0 sorry), `#print axioms` on all 8 declarations taken through the root import
and compared as sets, and the root import closure recomputed by hand at **20 modules, was 19**,
with `WindowBridge` among them.

The 09-02 caveat is unchanged and still stands: `nanoda-allow-sorry: true` means the badge
does not exclude a `sorry` on its own. `#print axioms` plus `axiom-audit` is what stands
behind this node. This module has no `sorry` to hide.

## Registry

Node **`Q63-window-bridge`** added to `proofs/registry/fock-ribbon-sign-operator.json` as a
sibling of `Q63-window-lemma` under `Q63-exponent-range-and-cross-runner` (which *is*
`cor:range`), graded **`lean-verified`**, `lean` field listing the six headline declarations,
`sources: ["math/9905196"]`.

The 09-02 `Q63-window-lemma` node said the step to `a ≤ ℓ−1` was *"still on paper only"*.
That sentence is now false, so it was **edited**, not left to calcify: it now names the step
as discharged and points at the sibling, while keeping the two halves that remain informal.

`registry_validate.py` → **1 remaining problem, pre-existing and not mine**:
`Q63-rigidity-dichotomy` claims `proved` with a `computed` child `Q63-formii-verification`,
in a different subtree (`Q63-level-ell-telescope`). Verified pre-existing by re-running the
validator on a copy with my node deleted. Re-grading it is a mathematical judgement, not a
Lean one, so I left it — flagging it here rather than silently touching a grade.

## Cost

Build cache was warm, so `lake build` ran in ~10–20 s, not the ~10 min the brief warned of;
direct `lean` round-trips were ~2 s. The corollaries type-checked first try. The session's
real work was reading `sec:closed` closely enough to see that `a` and `h` are the same
count along two axes of one label — which is why there are two theorems here and not one.
