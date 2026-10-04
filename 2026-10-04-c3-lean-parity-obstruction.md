# LEAN 2026-10-04 c3 — the parity obstruction, sorry-free

**Target declaration(s):** `TworowD4Kernel.ParityObstruction.*` (16 declarations)
**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/ParityObstruction.lean` (290 lines)
**Paper proof:** `proofs/2026-10-04-width-vector-M-convexity.tex`, `lem:P` and (H1)
**Registry:** `proofs/registry/cylindric-lorentzian.json`, new node
`A-parity-obstruction-lean` at `lean-verified`, child of
`A-width-vector-M-convexity-REFUTED` (stored `trust` **`proved`** before and after — a Lean
child does not promote its parent)
**Commits:** `tworow-d4-kernel @ 17be689a9de89c35fff3f3ff8846c68c3de0d3ab` (the Lean file),
`proofs @ 45ebe58fa5c096d36917abad10ff270dac2f18ba` (this note + the registry)
(resolved by `git log -1 --format=%H` **inside** `projects/lean/tworow_d4_kernel`, whose
`origin` is `github.com/clio-vega/tworow-d4-kernel`; the hash is not line-scoped and not
typed from memory)
**Lean/Mathlib:** Lean 4.30.0, Mathlib `v4.30.0`, Lake 5.0.0

---

## Result: everything builds, nothing is sorried

`lake build` (unpiped) exits **0**, 3183 jobs. **Sorry count: 0.** No `native_decide`, no
local axiom, no `#exit`.

All 16 declarations are sorry-free. `#print axioms` on all 12 theorems returns exactly

```
[propext, Classical.choice, Quot.sound]
```

— the standard three, no more.

## What the session actually did

The paper states `lem:P` and then *measures* it: 13,522 slices where the width-vector set is
M-convex, 22,451 where it fails, the failures coinciding exactly with the
multiple-coordinate-sum slices. That coincidence is not evidence for a mechanism; it **is**
the mechanism, and every ingredient was already sorry-free in
`TworowD4Kernel.HalfWidthL1`. This session closed that loop in Lean. No counterexample
search appears anywhere in the file.

### The honest accounting of what is new

| declaration | status | hypotheses **consumed** |
|---|---|---|
| `l1Dist_eq_sum_devSeq_add_two_mul_defect` | **new** | **none** |
| `l1Dist_emod_two_eq` | **new** | **none** |
| `exists_lt_of_sum_lt` | **new** | **none** |
| `l1Dist_eq_of_per` | new | `Per n m lam` only — **not** `Per n m nu` |
| `two_dvd_l1Dist_sub` | new | `Per lam`, slice `u m nu = u m nu'` |
| `widthSum_eq` | **a rename, not a result** | `Per lam`, `Per nu` |
| `sum_widthVec`, `sum_widthVec_eq` | new (plumbing to `Fin m → ℤ`) | as above |
| `l1_ex_sub` | **new** | `y' i < y i`, `y l < y' l` |
| `no_sum_lt_of_mConvex`, `sum_eq_of_mConvex` | **new — the real content** | one-sided exchange only |
| `not_mConvex_of_two_sums` | new | — |
| `mConvex_inS` | new (definition check) | — |
| `l1Dist_eq_of_mConvex`, `widthSet_not_mConvex` | new | `Per lam`, `Per nu`, `Per nu'` |
| `two_le_abs_sum_sub` | new | the above + slice condition |

Two things worth stating plainly rather than burying:

**`widthSum_eq` is a rename.** `HalfWidthL1.twoHalfWidth` is *defined* as `(∑ᵢ wᵢ) - m`, so
the paper's `∑ᵢ wᵢ = G + m - ‖y‖₁` with `G = n - m` is literally
`HalfWidthL1.twoHalfWidth_eq_sub_l1` with the definitional `- m` moved across — proof
`omega`. The brief told me to check this *before writing anything*, and it was right to:
the honest deliverable was to declare it anyway (so the paper's own expression occurs
verbatim in the development), say in its docstring that it is a rename, and spend the
session on the obstruction instead.

**The splitting was already written, unnamed.** `‖y‖₁ = ∑ᵢ yᵢ + 2k(ν)` sat as an anonymous
`have hsum` inside the proof of `HalfWidthL1.twoHalfWidth_eq_sub_two_mul_defect`, with
`Per lam` and `Per nu` in scope *and used by neither*. Extracting it as a top-level
declaration with **no hypotheses** is what makes the parity statement available at all. The
parity obstruction had been sitting inside a proof of something else for a day.

### The one genuinely new piece of mathematics

`sum_eq_of_mConvex`: **an M-convex set has constant coordinate sum.** By infinite descent —
if `∑ y' < ∑ y` then some `yᵢ > y'ᵢ`, the exchange axiom supplies `l` with `y l < y' l`, and
`l1_ex_sub` shows the surgery `y - eᵢ + e_l` spends **exactly 2** of the `ℓ¹` distance to
`y'` while preserving `∑ y` (`SublevelMConvex.sum_ex`). Induction on `N` with
`‖y - y'‖₁ ≤ N` closes it; the base case uses `‖y - y'‖₁ ≥ |∑(yᵢ - y'ᵢ)| > 0`.

It consumes only the **one-sided** axiom (B-EXC⁻), not Murota's symmetric one. That was
worth checking rather than assuming: it makes the obstruction apply to strictly more sets
than Theorem M's — including `MConvexExchange.InSupp`, for which only the one-sided axiom
is proved in this development.

## A false claim in my own brief

The brief's §2(4) said: *"`MConvexExchange.lean` and `SublevelMConvex.lean` already have the
exchange machinery — reuse their `MConvex` predicate, do not introduce a second one."*

**There is no `MConvex` predicate.** Grep for `def MConvex` / `structure MConvex` /
`abbrev MConvex` / `class MConvex` / `def IsMConvex` over the whole tree with
`command grep -rn`: **0 hits**. Despite the file name `MConvexExchange.lean`, neither file
ever abstracts the predicate — each proves an exchange axiom for one *concrete* set
(`InSupp` in one, `InS` in the other). The instruction named a declaration that does not
exist, and the token `MConvex` appears only in file names, docstrings and namespace names,
which is exactly what makes the claim read as true.

So `MConvex` is **defined here for the first time in this development**. The brief's own §6
is what caught this — *"if this brief asserts anything about what is already formalised,
check that assertion against the declarations, not against this brief"* — and it is worth
noting that the same brief was right about all five `HalfWidthL1` declarations in §1,
verified individually at their stated line numbers (227, 257, 289, 204, 109/128/132, plus
164, 176, 192, 96). **One brief, nine true claims about existing declarations and one false
one, with nothing in the brief's own prose distinguishing them.**

Because the predicate is new, a misstatement of it (an inequality the wrong way round, say)
would let every theorem above typecheck while applying to nothing. `mConvex_inS` is the
guard: it proves `MConvex {y | InS P Q σ j y}` from
`SublevelMConvex.sublevel_symm_exchange`, i.e. the predicate introduced here is the one
**Theorem M** already verifies.

## What is NOT formalised

Stated so the node is not read as covering more than it does. The parent node's statement
has two halves and a sharpness claim; I formalised one half.

1. **The M♮ (M-natural) half.** The parent asserts that an M♮-convex set has coordinate sums
   forming an *integer interval*. Not formalised; no M♮ predicate is defined.
2. **`A-general-m-tent`**, hence the parent's *"two coordinate sums at distance **exactly**
   2"*. `two_le_abs_sum_sub` proves only `2 ≤ |Δ|`, from the slice condition plus parity.
   The "exactly" needs the tent lemma (occurring half-widths are consecutive), which is
   `proved` on paper and **not formalised**.
3. The **converse** (13,522/13,522 single-half-width slices have `W` M-convex) is measured,
   not proved, and the paper says so. Nothing here touches it.

Per the standing distinction: items 1–3 are **proved** (on paper) and **not formalised**.
Not "not proved".

## Verification, by the methods that actually work

- **Build:** `lake build` run **unpiped**, exit 0, 3183 jobs. Not `lake build | tail`, whose
  exit status is `tail`'s.
- **Sorry scan:** `command grep -rn "sorry" TworowD4Kernel/` over the **whole tree**, not one
  directory deep. 11 hits, **all prose in docstrings**, 0 in code. The scanner was
  **validated against a planted canary** (`theorem planted_canary : False := by sorry`
  appended at line 290) — it found it, so the clean result means something. Canary removed
  and the restore re-verified by rebuild.
- **`command grep`, never bare `grep`,** throughout: the shell `grep` wraps
  `ugrep --ignore-files`, and `lean/.gitignore` line 6 is `tworow_d4_kernel/` — the entire
  active development. Bare `grep -r` returns zero for declarations written minutes earlier.
- **Load-bearingness by ABLATION, not by dependency walk.** Four citations deleted one at a
  time; each must break the build:

  | ablation | result |
  |---|---|
  | `sum_devSeq` in `l1Dist_eq_of_per` | build **failed** ✓ |
  | `twoHalfWidth_eq_sub_l1` in `widthSum_eq` (swapped for the sibling `twoHalfWidth_eq_sub_two_mul_defect`) | build **failed** ✓ |
  | `l1_ex_sub` in `no_sum_lt_of_mConvex` | build **failed** ✓ |
  | `sum_ex` in `no_sum_lt_of_mConvex` | build **failed** ✓ |
  | **control:** pristine file | **builds** ✓ |

  The control matters: without it, "all four failed" is consistent with a build that was
  broken for an unrelated reason.
- **Axioms:** `#print axioms` on all 12 theorems, via `lake env lean` on a scratch file.
  Uniformly `[propext, Classical.choice, Quot.sound]`.

## Hypothesis discipline

Everything is stated with `HalfWidthL1.Per`, never `GreedyChain.IsCylindric`. `Per` is the
weaker hypothesis (`Per.of_isCylindric`) and so the stronger theorem; `IsCylindric` would
formalise something strictly weaker than the paper's claim.

Trying to drop each hypothesis before using it paid off three times, at a cost of about a
minute: the splitting and the parity need **nothing** (not `Per lam`, not `Per nu`, not
`0 < m`), `l1Dist_eq_of_per` needs `Per lam` alone, and `sum_eq_of_mConvex` needs neither
`0 < m` (for `m = 0` both sums are `0`) nor any finiteness, boundedness or nonemptiness
assumption on `S`. None of this was in the brief's plan; all of it is in the docstrings.

## Citations carried into the Lean file

`HalfWidthL1`'s ambient model is WZZ's cylindric skew Schur setting,
**arXiv:2401.14632 §5**, transported to bead coordinates by `prop:dictionary` of
`proofs/2026-09-20-c1-cylindric-M-convexity.tex`, following Lam–Postnikov. The module
docstring of `ParityObstruction.lean` cites `lem:P` of
`proofs/2026-10-04-width-vector-M-convexity.tex` by label, and names
`SublevelMConvex.sublevel_symm_exchange` as the external result `mConvex_inS` reuses.

## Two small Lean notes, for next time

- `rw [widthVec, …]` fails with *"Failed to rewrite using equation theorems"* on a `def`
  returning a function; `simp only [widthVec]` works.
- `Finset.sum_sdiff` + `omega` left the pair sums unlinked (omega reported only
  `e - f ≥ 3` over opaque sum atoms). Mirroring the **already-proven** pattern of
  `SublevelMConvex.kfun_ex` — `Finset.sum_subset` to kill the off-support terms, then
  `Finset.sum_pair` — worked first try. When a sum manipulation fights back, copy the shape
  of a proof in the same file that already does it.
