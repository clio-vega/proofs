# LEAN 2026-10-06 c2 — the fibre-count behind Theorem A's `N(u,w)`

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**File:** `TworowD4Kernel/LabelledBlocks.lean` (new, wired into the root `TworowD4Kernel.lean`)
**Status:** sorry-free, full-tree `lake build` green (3199 jobs).

## Target

The mathematical content underneath Theorem A of `2026-10-06-ls-corollary-4-coupling.tex`
(`arXiv:math/0202090`, Lenart–Sottile, *Skew Schubert polynomials*, Corollary 4 and the
remark after Theorem 2): the number of ways to partition a finite set into labelled blocks
of prescribed sizes is a multinomial coefficient.

```lean
theorem card_blockFunctions_eq_multinomial {m : V → ℕ} (hm : ∑ v, m v = Fintype.card A) :
    #(blockFunctions A m) = Nat.multinomial univ m
```
where `blockFunctions A m = univ.filter (BlockSizes m)` and
`BlockSizes m f ↔ ∀ v, #(univ.filter (f · = v)) = m v`.

## Mandatory Mathlib pre-flight — result: NOT PRESENT (proceeded with the primary target)

Searched the **statement shape**, not names, with a Python declaration-block parser (a
line-based `grep` was tried first and is **not** adequate: it missed
`Multiset.countPerms_filter_ne`, which spreads `card` and `multinomial` across lines).

**Negative searches, in full:**

| search | result |
|---|---|
| decl blocks in all of Mathlib mentioning `multinomial` **and** `card`/`filter`/`fib` | **3**: `countPerms_filter_ne`, `multinomial_nsmul`, `dpow_sum'` — none is a fibre count |
| `card_eq_multinomial`, `multinomial_eq_card`, `card_multinomial` | **0 hits** |
| decl statements with `.card` + `factorial` + `∏` | 13, all permutation-cycle-type / divided-power; none counts functions |
| `piAntidiag` (all 18 decls) | functions with prescribed **sum**, not prescribed fibre cardinalities; no card lemma |
| decl statements counting a subtype/filter of **functions** | 3, all unrelated (`bijective_iff_*_and_card`) |
| `Finset.card_pi` / `Fintype.card_pi` / `card_filter_piFinset_*` | unconstrained or single-coordinate fibres only |
| orbit-cardinality lemmas (`card_orbit`, `ncard_orbit`) | generic orbit-stabiliser only |
| `Multiset.bell` / `Nat.uniformBell` (`Combinatorics/Enumerative/Bell.lean`) | the **unlabelled** sibling, and **not proved to count anything** |

The decisive evidence is Mathlib's own TODO in `Mathlib/Combinatorics/Enumerative/Bell.lean`:

> `## TODO` / Prove that it actually counts the number of partitions as indicated.

`Multiset.bell` is *defined* as an arithmetic expression and *documented* as "number of
partitions of a set of cardinality `m.sum` whose parts have cardinalities given by `m`";
likewise `Multiset.countPerms` is defined as a multinomial and named "the number of
permutations of a given multiset", with no theorem that it counts permutations. So Mathlib
has the **arithmetic** of multinomials in full (`Nat.multinomial`, `Nat.multinomial_spec`,
`Finset.sum_pow`) and explicitly lacks the **combinatorial interpretation**. The target is
the labelled (ordered-blocks) sibling of that TODO. **Not a restatement.**

**Second axis — what the *proof* uses** (`is-it-in-mathlib-is-indexed-by-which-form`): here
Mathlib *did* hold the key ingredient, under vocabulary the statement does not mention —
`DomMulAct.stabilizer_card` (`Mathlib/GroupTheory/Perm/DomMulAct.lean:99`):
`Fintype.card {g : Perm α // f ∘ g = f} = ∏ i, (Fintype.card {a // f a = i})!`. That is the
∏(mᵥ)! denominator. The numerator/orbit half is what was missing, and is what I proved.

## What builds, sorry-free

| declaration | content |
|---|---|
| `fibre`, `mem_fibre`, `card_fibre_eq_card_subtype` | the fibre of `f` over `v` as a `Finset` |
| `BlockSizes`, `blockFunctions`, `mem_blockFunctions` | the prescribed-size condition and its `Finset` |
| `sigmaFstFibreEquiv` | fibre of `Sigma.fst` in `Σ w, Fin (m w)` is `Fin (m v)` |
| `exists_equiv_sigma` | any `f` with block sizes `m` identifies `A ≅ Σ v, Fin (m v)` over `V` |
| `exists_blockSizes` | **existence**: `∑ v, m v = #A → ∃ f, BlockSizes m f` |
| `exists_perm_comp_eq` | **transitivity**: `blockFunctions A m` is one `Perm A`-orbit |
| `blockSizes_comp` | block sizes invariant under precomposition by a permutation |
| `permFibreEquivStabilizer` | every fibre of `σ ↦ f₀ ∘ σ` is a translate of the stabiliser |
| `card_filter_comp_eq_prod_factorial` | each fibre has `∏ v, (m v)!` elements (via `DomMulAct.stabilizer_card`) |
| `card_blockFunctions_mul_prod_factorial` | `#(blockFunctions A m) * ∏ v, (m v)! = (#A)!` |
| `card_blockFunctions_eq_multinomial` | **the target** |

**Sorries: 0.** Scanner validated against a planted canary: canary inserted → scan reported
1 hit; canary removed → `sha256` byte-identical to before (`5150625786…`), scan reported 0.
Scan was whole-tree, not one directory deep. The 11 bare `grep sorry` hits in the project
are all prose inside docstrings; the tactic-shaped scan (`by sorry|:= sorry|sorry$`) is 0.

## Axioms

All declarations: `[propext, Classical.choice, Quot.sound]` — the standard three.
`sigmaFstFibreEquiv` depends on **no** axioms. `Classical.choice` enters through
`Fintype.equivFinOfCardEq` / `Fintype.equivOfCardEq`, which are `noncomputable` in Mathlib;
this is expected for a counting argument that picks an identification with a canonical model.

## Ablation evidence

`DomMulAct.stabilizer_card` is load-bearing, verified the way an ablation has to be — pattern
asserted **present first**, then mutated, then the non-zero exit checked:

- pre-flight `grep -c "DomMulAct.stabilizer_card f₀"` → **1** (asserted `== 1`, else abort)
- after `sed` → `grep -c stabilizer_card_ABLATED` → **1** (the mutation really landed, no cache replay)
- ablated build → `Unknown constant DomMulAct.stabilizer_card_ABLATED` + `unsolved goals`,
  **exit 1** (read via `PIPESTATUS[0]`, not the pipe's last element)
- restored → **exit 0**, 1153 jobs

No two-proofs-of-one-Prop comparison was used: by definitional proof irrelevance that check
is constant and cannot fail.

## Registry

**No registry node covers this result.** All `proofs/registry/*.json` trees were walked for
`multinomial` / `corollary-4` / `thm:A` / `theorem a` / `labelled` / `N(u,w)` / `fibre count`
— **zero matches**. So nothing was graded `lean-verified`; there was nothing to grade.
`samuel-chain-formula.json` matched only on a filename-level keyword grep and has 0 nodes
(its schema is `conjecture`/`tree`, no node for Theorem A).

## Scope — what this is NOT

Theorem A's full product formula
`N(u,w) = ∏_α I_α(u,w)! / ∏_v (c^w_{u,v}!)^{I_α(w₀v, w₀)}`
is **not** formalised. What is formalised is its single mathematical input: one multinomial
per type `α`. The remaining content of Theorem A is indexing — a chain has exactly one type,
so `Γ(u,w)` splits as a disjoint union over `α` and the conditions for distinct `α` constrain
disjoint parts of the data; the exponent `I_α(w₀v,w₀)` counts the independent `v`-components.
Formalising that indexing would require formalising labelled Bruhat order, chain types and
`Γ_α`, none of which exists in Lean here. **`unproved ≠ unformalised`, and formalising the
mathematical core is not formalising the parent.** I did not promote anything on its behalf.

## Independent corroboration + a negative control that fires

Checked by `decide` (brute-force enumeration of all `3^4` functions `Fin 4 → Fin 3` — a
*different mechanism* from the orbit-stabiliser proof, so this is corroboration rather than
a second reading of one engine):

| `m` | `∑ m` | `#(blockFunctions (Fin 4) m)` | `Nat.multinomial univ m` |
|---|---|---|---|
| `(2,1,1)` | 4 = `#A` | **12** | **12** — agree |
| `(2,1,0)` | 3 ≠ `#A` | **0** | **3** — disagree |

The second row is the negative control, and it **fires**: without `hm : ∑ v, m v = #A` the
theorem is false, so the hypothesis is load-bearing rather than decorative. All four
`decide` goals close (`CHECK_EXIT=0`).

## Push state (verified, not asserted)

`git -C /home/clio/projects/lean/tworow_d4_kernel` — HEAD `dd06193`, pushed to
`origin/main` (`763895f..dd06193`). Verified by `fetch` + `rev-list --count @{u}..HEAD` = **0**
+ `branch -r --contains HEAD` = `origin/main`. The sibling repo
`/home/clio/projects/lean` was checked with its **own** `git -C` call: 0 unpushed, clean tree
(it needed no change). Every git call used `git -C <abspath>`; no `cd` + two-pushes.

**Caveat on this note:** `/home/clio/projects` is *not* a git repository, so this `.md` is
local to the Docker volume and **Robin cannot read it**. Only the Lean file is shared, at
https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/LabelledBlocks.lean
