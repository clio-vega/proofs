# LEAN 2026-10-22 c2 — the two-term boundary law, and the colength that cannot see its summands

**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (inner git repo)
**Module:** `TworowD4Kernel/BoundaryTwoTerm.lean` (new, 274 lines)
**Source artifact:** `proofs/2026-10-22-boundary-valuation-two-term.tex` — `thm:main`, `cor:inv`, `rem:recovers`
**Commit:** `7581a06` — verified by `git show --stat HEAD`, which lists
`TworowD4Kernel.lean | 1 +` and `TworowD4Kernel/BoundaryTwoTerm.lean | 274 +`.
Not by the hash: `git add` on an ignored path short-circuits `&&` and leaves `rev-parse`
naming the previous commit.

**Result: sorry-free.** 21 source declarations, 0 sorries, standard three axioms only.

---

## What type-checked

| label | declaration | statement |
|---|---|---|
| L1 | `boundary_count_two_term` | `#{a ∈ M : q^[m∉I]·ε^a = 1} = if m ∈ I then fM M ε else gM M q ε`. **No hypothesis on `q`.** Over any `[Monoid K]`. |
| L2 | `boundary_count_indicator` | the same, as `[m∈I]·f_M(ε) + [m∉I]·g_M(q,ε)` — the form the paper displays and the colength sums |
| — | `gM_one` | `g_M(1,ε) = f_M(ε)`: the 10-21 `q := 1` lesion, read off `eq:main` instead of measured |
| S4 | `ne_inv_pow_iff` | `q ≠ (ε⁻¹)^h ↔ q·ε^h ≠ 1` over a `GroupWithZero` — the paper's form vs the multiplicative one |
| S4 | `gM_singleton_eq_zero_iff` | for `M = {h}`: the 10-21 binder `q ≠ ε^{-h}` **is** `g_M(q,ε) = 0` |
| S4 | `thmBD_of_two_term` | **old law from new:** under that binder, `boundaryCount = [m∈I]·f_M(ε)`, i.e. (BD) |
| L5 | `indicator_add_compl` | `[m∈I] + [m∉I] = 1` pointwise — the two summands *partition* |
| — | `filter_nmem_powerset` | `{I ⊆ s : m ∉ I} = powerset (s.erase m)` — the engine for both colengths |
| L3 | `colength_nmem_card`, `colength_nmem` | `Σ_{I ⊆ range N} [m ∉ I] = 2^(N−1)` |
| L3 | `colength_mem_card`, `colength_mem` | `Σ_{I ⊆ range N} [m ∈ I] = 2^(N−1)` |
| **L4** | **`total_cannot_witness`** | **equal totals ∧ distinct summand functions** — the pre-registered statement, verbatim |
| L4+ | `total_blind_supports_split` | the sharper form: equal totals ∧ disjoint, both-nonempty supports |
| S2 | `summands_differ` | conjunct 2 of L4, **hypothesis-free** |
| S1 | `colength_mem_of_nmem` | negative control: drop `m ∈ s` and the total is `0` |
| S1 | `colength_mem_lesion_witness` | numeral control: `0 ≠ 4` at `N = 3`, `m = 5` |
| S3 | `total_cannot_witness_joint_lesion` | drop `hm` *and* `hN`: conjunct 1 is **false**, by `decide` |

**Sorries: none.** Not read off a `grep` for the prose word "sorry" (which false-positives on
"admitted") but from an environment scan: `SORRIED 0` over every constant whose value
`hasSorry`, plus `#print axioms` on all 18 theorems individually.

## Declaration set, both directions

* Source: **21** named `theorem`/`def`.
* Environment (`env.constants`, namespace-prefixed, non-internal): **23**.
* Difference, both ways: env ∖ source = `{fM.congr_simp, gM.congr_simp}` — auto-generated.
  source ∖ env = ∅.

A count alone would not have said which way the discrepancy ran. It was +2 and benign;
on 10-19 it was −8 and benign, on 10-18 it was an over-count. Only the set says.

## Axioms

All 18 theorems: `[propext, Classical.choice, Quot.sound]`. Four of them need only
`[propext, Quot.sound]` (`boundary_count_two_term`, `boundary_count_indicator`,
`indicator_add_compl`, `gM_one`) — the two-term law and complementarity are
choice-free. Nothing outside the standard three; nothing I wanted and could not have.

```
'TworowD4Kernel.BoundaryTwoTerm.total_cannot_witness' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.BoundaryTwoTerm.boundary_count_two_term' depends on axioms: [propext, Quot.sound]
```

## Build

`lake build` on the whole library: `✔ [3477/3478] Built TworowD4Kernel (51s)`,
`Build completed successfully (3478 jobs)`. Read from the log, not from `$?` — a wrapped
command's exit status carries no information in either direction (10-19 reported "exit 0"
on a failed build; 10-06 was the reverse). The only warnings in the run are pre-existing,
in `ParityObstruction` and `IotaAdjunction`; `BoundaryTwoTerm` builds with none.

## Olean artifacts

`.lake/build/lib/lean/TworowD4Kernel/BoundaryTwoTerm.olean` (290424 bytes, 22:02 UTC) and
`.lake/build/lib/lean/TworowD4Kernel.olean` (11304 bytes, 22:04 UTC — rebuilt because the
root aggregator now imports the new module).

---

## The lesion table, graded against what Lean actually said

| # | hypothesis / pair | predicted | **measured** | how Lean said it |
|---|---|---|---|---|
| S1 | `m ∈ s` in L3 | IRREMOVABLE | **IRREMOVABLE** ✓ | `colength_mem_of_nmem`: the total is `0`; numeral control `0 ≠ 4` |
| S2 | `2 ≤ N` in L4 | IRREMOVABLE for conjunct 2 | **REMOVABLE** ✗ | the `unusedVariables` linter fired on `hN` |
| S3 | the pair (S1, S2) | jointly irremovable; S2 survives independently | **jointly irremovable, but only because S1 is** — partly ✗ | `total_cannot_witness_joint_lesion`, by `decide` |
| S4 | any hypothesis on `q` in L1 | ABSENT, not removable | **REMOVED — discharged**, not RELOCATED | `thmBD_of_two_term` derives the binder-carrying law from the binder-free one |
| S5 | `ε` a root of unity | ABSENT | **ABSENT** ✓ | no binder exists to lesion; the `Finset.filter` predicate carries it |

### S2 is worse for the pre-registration than "removable"

I predicted `2 ≤ N` was irremovable *for conjunct 2* and wrote, in the brief, "at `N = 1`
the powerset is `{∅, {0}}`; check by hand whether the witness survives". It survives — and
so would a witness at `N = 0`. But the measurement went further than the prediction's axis.
Conjunct 2 does not use `hm` either. For *any* `I`, exactly one of `[m∈I]`, `[m∉I]` is `1`,
so they can never be equal, and `I = ∅` is in every powerset. **Conjunct 2 is
hypothesis-free** — `summands_differ (N m : ℕ)`, no binders at all.

So the shape of L4 is not what I registered. I registered it as a conjunction of two facts
and expected the hypotheses to be shared. In truth:

- conjunct 1 (**equal totals**) carries all the content and needs `hm`;
- conjunct 2 (**distinct summands**) is free.

Which inverts the rhetoric of the brief. I wrote that conjunct 2 "is the content" and "is
what no census of the total can ever see". It is indeed what the census cannot see, but it
is *cheap* — true for trivial reasons, including when `m` is not in the ambient set at all.
The fact worth having is the one I thought was the easy half: that two functions this
obviously different have the *same total*. `total_blind_supports_split` is the statement I
should have pre-registered — equal totals, plus supports that are disjoint and both
nonempty. That version does need `hm` (for the `I = {m}` witness) and still not `hN`.

### S3: this cell did **not** repeat this morning's shape

This morning (10-22 c1, Q546) `0 ∉ M` and the range bound were each Lean-verified
**removable** and **jointly false**, because each removability proof consumed the other
hypothesis — a register of marginals was not a register of the joint. I armed S3 expecting
the same trap and wrote "this is the cell where I expect to be wrong."

I was wrong, but in the other direction. Here S1 is irremovable *on its own*, so the joint
verdict is forced by the marginal and adds nothing. **A register of marginals was sufficient
in this cell.** The lesson that generalises is therefore narrower than "always lesion the
pairs": pairs can hide a joint failure only when *every* marginal survives. When one marginal
is already fatal, the pair is uninformative. The cheap pre-test is to lesion the marginals
first and only escalate to pairs if they all survive.

### S4: the trap, discharged rather than deferred

`boundary_count_two_term` has no binder on `q`, and I pre-committed to grading S4
`RELOCATED` if I could not locate the paper lemma that *states* the hypothesis — because a
missing binder and a removable binder give the same green build, so the clean compile would
have demonstrated nothing. (That is the 10-19 finding: `pairwise distinct` pre-registered
IRREMOVABLE, found ABSENT, its binder living one lemma upstream.)

I did reach it. The binder is 10-21's `thm:BD` hypothesis `q ≠ ε^{-h}`;
`gM_singleton_eq_zero_iff` shows it *is* the locus `g_M = 0`, and `thmBD_of_two_term`
derives (BD) from the two-term law under it. The removal is therefore type-checked, not
merely unmentioned. Grade: **REMOVED**.

### The free lesion pass I keep not reading — and it fired again

Lean's linters produced **two** verdicts I had not pre-registered:

1. `unusedVariables` on `hN` — which *is* the S2 verdict. I kept the pre-registered
   statement verbatim and silenced the linter for that one declaration with a comment
   saying the warning is the verdict, rather than quietly renaming `hN` to `_hN`. (I had
   in fact written `_hN` on the first pass, which suppressed the finding; restoring the
   name is what put it on the record.)
2. `unusedDecidableInType` on `[DecidableEq K]` in `ne_inv_pow_iff` — graded **ABSENT**,
   and removed. I had assumed every declaration in the file needed decidable equality
   because the `Finset.filter`s do; the inverse-form lemma does not mention a filter.

Two of five cells in the table moved, and the linter supplied both. On 10-20 the only two
lesions that moved were also the ones the linter deleted. That is now twice; it is not an
anecdote.

---

## Scope: what is **not** formalised

The module is named after a valuation and does not formalise one. Explicitly **not**
Lean-verified, and not to be inferred from the declaration names:

- that `ord_{S^ε_m} Ã^X_{I,I}` equals `boundaryCount` **at all** — the factor-by-factor
  valuation argument over `eq:inv`/`eq:T`/`eq:diag`, resting on `lem:unit`;
- `eq:inv`/`eq:T`/`eq:diag` themselves, which for `X ∈ {B,C}` are **Motegi's theorem**
  (`arXiv:2610.15344`, L1305–1335 type `B`, L998–1046 type `C`) and for `X = BC` are my
  own **definition** with no weight function in the literature;
- the `(B_N, ε=−1, q₀=−1)` specialisation that makes `cor:inv` a statement about a
  geometric object rather than about `Finset ℕ`.

This gap is now a **registry node**, `clio-Q544-root-system-layer-NOT-formalised`, not a
remark in a docstring — a dependency is only graded once it is a node (10-20). It is what
caps the parent at `proved`.

## Registry

`proofs/registry/rick-q488-converse.json`, backup `*.bak-1022c2-lean`. Added one parent
with four children, every node carrying its `role`:

| node | trust | role |
|---|---|---|
| `clio-Q544-thm-main-boundary-two-term` | **`proved`** (not promoted) | `premise` |
| ├ `clio-Q544-two-term-arithmetic-layer` | `lean-verified` | `premise` |
| ├ `clio-Q544-thmBD-binder-is-the-locus-gM-zero` | `lean-verified` | `premise` |
| ├ `clio-Q544-colength-cannot-witness-its-summands` | `lean-verified` | `premise` |
| └ `clio-Q544-root-system-layer-NOT-formalised` | `proved` | `premise` |

`python3 code/registry_validate.py proofs/registry/rick-q488-converse.json --proofs-dir .`
→ `OK: ... is valid (status: in-progress)`, **0 problem lines**. All **20** live registries
re-validated: 20 OK, 0 not-OK, so the global clean state at `c301ba1` is preserved.

`trustcheck.py` was again the wrong instrument: in this container it takes subcommands
(`validate|sources|report|ops`) and a required `--deployment`, so a registry path makes
argparse exit 2 with a *usage* message. My problem-line grep counted `1` — which is the
usage error, not a registry fault. A remembered flag belongs to a script, not to a task;
`registry_validate.py --proofs-dir .` is the registry checker here.

## Caveats travelling with the node

`G2` (non-cancellation at fixed `χ`) untouched. The `BC` arms rest on verification routes
A+B only; route C is inapplicable. The shared RTV stable-envelope root under `eq:diag` is
unreduced by this work. And the standing downstream obligation: the **form** `2^{N−1}` must
not be inherited under any `q`-specialisation into §§5–7 of
`2026-10-19-c2-sa-divisors-full-q.tex` — `colength_mem`/`colength_nmem` are theorems about
`Finset ℕ` and carry no `q` at all, which is exactly why they cannot license that inheritance.
