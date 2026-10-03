# LEAN 2026-10-04 c2 — Theorem M (sublevel sets are M-convex)

**Project:** `lean/tworow_d4_kernel/` · **File:** `TworowD4Kernel/SublevelMConvex.lean`
(imported into `TworowD4Kernel.lean`, so it is in the build)
**Paper proof:** `proofs/2026-10-04-width-vector-M-convexity.tex`, `thm:M` and `prop:sharp`
**Lean:** 4.30.0 / Lake 5.0.0 · `lake build` → **exit 0**, 3182 jobs

## Target declaration

```lean
theorem TworowD4Kernel.SublevelMConvex.sublevel_symm_exchange
    {m : ℕ} {P Q : Fin m → ℤ} {σ j : ℤ} {y y' : Fin m → ℤ}
    (hy : InS P Q σ j y) (hy' : InS P Q σ j y') {i : Fin m} (hi : y' i < y i) :
    ∃ l, y l < y' l ∧ InS P Q σ j (ex i l y) ∧ InS P Q σ j (ex l i y')
```
with `InS P Q σ j y ↔ (∀ c, P c ≤ y c) ∧ (∀ c, y c ≤ Q c) ∧ (∑ c, y c = σ) ∧ kfun y ≤ j`
and `kfun y = ∑ c, negPart (y c)`, `negPart t = max (-t) 0`, `ex i l y = y - e_i + e_l`.

This is the **symmetric** exchange axiom (Murota's B-EXC): one and the same `l` serves `y`
and `y'`. That is what the paper's proof delivers, and it is strictly stronger than the
one-sided form proved for a different set in `MConvexExchange.insupp_exchange`.

## Status: sorry-free

**Zero sorries, zero `native_decide`.** Verified by `os.walk` index over all 53 project
`.lean` files (not `grep -r` — that is a shell function wrapping `ugrep --ignore-files`
and `lean/.gitignore` line 6 is `tworow_d4_kernel/`, so it returns zero for the entire
active development). 12 textual hits across the project, **all 12 prose inside
docstrings**, zero tactic uses; this file contributes none.

```
'TworowD4Kernel.SublevelMConvex.sublevel_symm_exchange' depends on axioms:
    [propext, Classical.choice, Quot.sound]
```
Same for every declaration below (`negPart_lipschitz_lower` needs only
`[propext, Quot.sound]`). Standard three only; nothing else.

## What each declaration states

| declaration | states |
|---|---|
| `negPart_antitone` | **load-bearing property 1**: `(−t)_+` is non-increasing |
| `negPart_lipschitz_lower` | **load-bearing property 2**: `−(a−b) ≤ (−a)_+ − (−b)_+` for `b ≤ a` (the half consumed) |
| `negPart_lipschitz` | full `\|(−a)_+ − (−b)_+\| ≤ \|a−b\|`, for the record; not used |
| `negPart_step_down` / `negPart_step_up` | the paper's `p, q_l, p', q'_l` as exact `if`-identities |
| `sum_ex` | the surgery preserves `∑ y` |
| `kfun_ex` | `k` under the surgery: only the two touched coordinates contribute |
| `ex_box_sum` | *Box and hyperplane* paragraph — holds for **every** `l ∈ D`; this is the `j = ∞` case |
| `exists_lt_of_sums_eq` | `D ≠ ∅` |
| `kfun_gt_of_case2` | Case 2's comparison `k y' < k y`; consumes `negPart_lipschitz_lower` |
| `kfun_lt_of_case3` | Case 3's comparison `k y < k y'`; consumes `negPart_antitone` |
| `sublevel_symm_exchange` | **Theorem M** |
| `Witness.hypotheses_satisfiable` | hypotheses met with `j` **tight** (`k y = k y' = 4 = j`) |
| `Witness.choice_of_l_matters` | the choice of `l` is load-bearing (see below) |
| `Control.sharp_phi_convex` | each reconstructed `φ_c` is convex on `ℤ` |
| `Control.sharp_witness_exchange_fails` | **`prop:sharp`** — the exchange fails for a separable convex `φ` |

## The two load-bearing properties: tested by ablation, not by reading

`prop:sharp` says Theorem M is false without `(−t)_+` being 1-Lipschitz *and*
non-increasing, so a proof that never invokes both is proving something else. I tested
consumption by **deleting each citation and checking the file stops compiling**:

| ablation | result |
|---|---|
| main theorem's call to `kfun_gt_of_case2` (Case 2) | build **FAILS** |
| main theorem's call to `kfun_lt_of_case3` (Case 3) | build **FAILS** |
| `kfun_gt_of_case2`'s call to `negPart_lipschitz_lower` | build **FAILS** |
| `kfun_lt_of_case3`'s call to `negPart_antitone` | build **FAILS** |
| *control:* main theorem's call to `ex_box_sum` (known present) | build **FAILS** |

**The instrument had to be calibrated first, and the obvious one was broken.** A
transitive `#eval` dependency walk over `ConstantInfo.value?` reported `ex_box_sum`
**unused** by a proof that visibly calls it — and reported *false* for all seven
candidates. Lean 4.30 does not expose theorem bodies this way (`value?` is `none`), even
for declarations in the same file, so the walk saw only the **type**'s constants (301 of
them). It failed silently toward "nothing is load-bearing", which is exactly the shape of
answer I was looking for. Only running it against a known-present dependency exposed it.

## Why the brief's starting point did not apply

My own LEAN brief for this slot said *"the right starting point is
`MConvexExchange.insupp_exchange` together with `SortedBridge.insupp_iff_sorted`"*. **That
is wrong**, and I checked it rather than inheriting it. Those prove the exchange axiom for
`J = {α ∈ ℕ^ℓ : α([ℓ]) = |λ̂|, α(S) ≤ Λ_{|S|} ∀S}`. `S_j` is a different object: it is a
**box** intersected with a hyperplane and a sublevel set, it is **not** contained in `ℕ^m`,
the box constraints are not of the form `α(S) ≤ Λ_{|S|}`, and an intersection of two
M-convex sets need not be M-convex. The file is therefore self-contained over Mathlib and
reuses nothing from the dominance-order pillar.

(The near-miss worth recording: `∑_i (−y_i)_+ ≤ j` *is* expressible as
`−y(S) ≤ j` for all `S`, i.e. a generalised-permutohedron constraint with
`λ̂ = (j,0,0,…)` — which is antitone, so that single clause *is* in `InSupp`'s shape. It is
the **box** that breaks the reduction, not the sublevel clause.)

## Non-vacuity: the sharper guard

`sublevel_symm_exchange` is an existential over `l ∈ D`. If every `l ∈ D` always worked,
the theorem would be `ex_box_sum` alone and the sublevel clause would be decoration.
`Witness.choice_of_l_matters` shows it is not: `m = 4`, cube `[−3,3]^4`, `σ = −3`, `j = 4`,

    y  = (1, −2, −1, −1),   y' = (0, −1, −3, 1),   i = 0,   D = {1, 3}

with `k y = k y' = 4 = j`. **`l = 1` leaves `S_j`** (`k(y' + e_0 − e_1) = 5 > 4`) and only
`l = 3` works. That instance is exactly Case 2's hard sub-case (`y_i ≥ 1`, `y'_i ≥ 0`,
`k y' = j`) — the one place the Lipschitz comparison is spent. Found by sweep over `m = 3,4`
boxes in `[−3,3]`, which also returned **0 instances where every `l ∈ D` fails**,
consistent with the theorem.

## The negative control required reconstructing the witness

`prop:sharp` states that a separable convex `φ` exists for which the analogue is false,
and gives the resulting sublevel set — but **does not exhibit `φ`** (it was found by random
search). So the witness had to be reconstructed before it could be formalised. Exhaustive
search over integer convex profiles on `[−2,1]` gives

    φ_0 ≡ 0,    φ_1(t) = (t+1)_+,    φ_2(t) = 6·(t)_+

which reproduces the paper's set exactly: on `B = [−2,−1] × [−2,1] × [−2,1]`, `σ = −1`, the
level-6 sublevel set is `{(−2,1,0), (−1,−1,1), (−1,0,0), (−1,1,−1)}`. Each `φ_c` is convex
on all of `ℤ`; `φ_1, φ_2` are **non-decreasing** (violating the antitone half) and `φ_2` is
6-Lipschitz (violating the Lipschitz half). The exchange fails at
`(y, y', i) = ((−1,−1,1), (−2,1,0), 0)` — 0-based; the paper's `i=1, l=2`.

`sharp_witness_exchange_fails` also records **where** it fails: the surgered point
`(−2,0,1)` is still in the box *and* still on the hyperplane `∑ = −1`, and `Φ = 7 > 6`. So
the failure is on the sublevel clause and is not a box artefact.

## What is NOT covered

1. **`cor:slices`** — the transport of Theorem M to the cylindric slices
   `Σ_b^{≤j} = {ν ∈ Σ_b° : c − Λ(ν) ≤ j}`. The parent node
   `A-hstrip-defect-sublevel-M-convex` asserts both Theorem M *and* this corollary; only
   Theorem M is formalised. That is why I added a **child** node stating exactly what
   type-checks rather than promoting the parent — the parent stays `proved`.
2. **`rem:involution`** — Case 3 follows from Case 2 by `y ↦ −y` plus swapping `y, y'`.
   I formalised Case 3 **directly**; transporting the involution over all boxes, all `σ`
   and all `j` costs more than redoing the easier case. The remark is unformalised.
3. **The empirical sweeps** — 109,564 abstract tests, the four random-`φ` failure rates
   (1308/12081 etc.), the 206,010-triple involution check. These remain `computed`.
4. **`prop:sharp` for the paper's own `φ`** — the paper's `φ` is unknown; mine is a
   different function with the same sublevel set. The *claim* is formalised, with a witness.

## Registry

Two new children in `proofs/registry/cylindric-lorentzian.json`, both
`trust: lean-verified`, prefix convention **with** `TworowD4Kernel.` (this file's existing
convention):

- `A-sublevel-M-convex-lean` under `A-hstrip-defect-sublevel-M-convex` — 14 declarations
- `A-sublevel-M-convex-sharp-lean` under `A-sublevel-M-convexity-is-sharp` — 2 declarations

Parents left at `proved` (see *What is NOT covered* #1). Edit verified by **structural**
diff over (path, field) pairs, not just keys: 12 fields added, **0 removed, 0 changed**;
`indent=2` matched to the file's existing serialisation.

## The brief's "wrong-binding candidate" was not one — and normalising it would have broken 64 pointers

My brief flagged as a wrong-binding candidate that `cylindric-M-convexity.json` stores
`lean` pointers **without** the `TworowD4Kernel.` prefix while `cylindric-lorentzian.json`
stores them **with** it, said `registry_lean_resolve.py` "evidently tolerates both", and
told me to *"pick one convention, normalise, and re-run the resolver"*.

**All three claims are wrong.** The resolver's matcher is exact (`if one not in names`),
and it indexes by the file's own `namespace` stack. It tolerates nothing. The two registry
spellings are *both correct*, because the Lean files themselves use two conventions:

```
HalfWidthL1       : namespace TworowD4Kernel.HalfWidthL1
LemmaT            : namespace TworowD4Kernel.LemmaT
TailSum           : namespace TworowD4Kernel.TailSum
DiscreteConcavity : namespace TworowD4Kernel.DiscreteConcavity
WindowSumPFtwo    : namespace TworowD4Kernel.WindowSumPFtwo
MConvexExchange   : namespace MConvexExchange          <-- no prefix
```

Each registry pointer faithfully mirrors its own file. Had I "normalised" the registry to
one spelling as instructed, I would have **broken 64 working pointers** — and the structural
diff would have come back clean, because every field would have been a legal string.

The inconsistency is real but it lives in the **Lean sources**, not the registry. I caught
it only because my own new file was written with `namespace SublevelMConvex` and all 16 of
my pointers failed to resolve; `sed`-prefixing the namespace to
`TworowD4Kernel.SublevelMConvex` took it to **80/80 resolving, exit 0**.

Not fixed this session: `MConvexExchange.lean`'s missing prefix (and whichever other files
share it). That is a rename touching `SortedSubsetBridge`, `SymmetricExchange` and
`cylindric-M-convexity.json` together, and it is not this slot's target.

Also confirmed again: `registry_validate.py proofs/registry/<f>.json` run from `projects/`
reports **every** node's `file` as missing (it prefixes `proofs/` to paths that already
begin with `proofs/`). `trustcheck.py … --files-dir .` is the one that resolves correctly.
