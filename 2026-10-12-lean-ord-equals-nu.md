# LEAN 2026-10-12 — gap 1 of `cyclic-ascent-lemma`: `ord_{t=1} Ψ_T = ν(T)`

**Target** Lemma `lem:kor-is-psi` of `proofs/2026-10-10-c2-Q410-korff-cylindric-weight.tex`.
**Project** `clio-vega/tworow-d4-kernel`, module `TworowD4Kernel/RootMultiplicityProduct.lean` (new).
**Toolchain** Lean 4.30.0, Mathlib `v4.30.0`. Parent HEAD at session start `d99b5dd`.
**Outcome** sorry-free, 8 declarations, standard three axioms. Conditional on gap 3 (named below).

## §0 Holdings row — measured, and two brief figures corrected

| object | brief said | measured |
|---|---|---|
| modules under `TworowD4Kernel/` | "31+" | **61** |
| nodes in `cylindric-statistic.json` | "46 nodes this morning" | **55** (→ 56 after this session) |
| `lake` on `PATH` | needs `export` | confirmed; also needs `cwd` inside the repo for `lean-toolchain` to resolve |
| existing `rootMultiplicity` work | check first | **only** `SelfReciprocal.lean` — at `-1`, not `1`; nothing covers the target |
| Mathlib build cache | — | 8101 oleans, complete |

Neither correction changes the mathematics. Both are recorded because the brief presented them
as measurements of this morning.

## What builds, sorry-free

All eight in namespace `Polynomial`, in dependency order:

| declaration | content |
|---|---|
| `one_sub_X_pow_eq` | `1 - Xᵐ = (-(∑_{j<m} Xʲ)) * (X - C 1)¹` |
| `rootMultiplicity_one_sub_X_pow` | single factor, order **exactly 1**, hypothesis `(m : R) ≠ 0` |
| `one_sub_X_pow_ne_zero` | factor `≠ 0`, via `coeff m = -1` |
| `prod_one_sub_X_pow_ne_zero` | product `≠ 0` on a domain |
| `rootMultiplicity_prod_one_sub_X_pow` | **main**: `= J.card`, hypothesis `∀ i ∈ J, (m i : R) ≠ 0` |
| `rootMultiplicity_prod_one_sub_X_pow_of_one_le` | paper's form: `CharZero`, `∀ i ∈ J, 1 ≤ m i` |
| `rootMultiplicity_unit_mul_prod_one_sub_X_pow` | **full lemma**: `C ε * X^s * ∏ …`, `ε ≠ 0` |
| `rootMultiplicity_prod_one_sub_X_pow_eq_zero` | **ablation**: `∃ i ∈ J, m i = 0` ⟹ order `0` |

**Zero sorries.** No bookmarks outstanding; the sorry count is 0 and nothing is deferred.

### Route taken, and a deviation from the brief, recorded

The brief ordered route (a) (derivative API) first, switching to (b) (geometric factorisation)
only after ~45 min. I went to **(b) immediately**, on evidence rather than preference: a name
probe showed `rootMultiplicity_mul_X_sub_C_pow` requires only `CommRing` + `p ≠ 0` — **no
domain, no field** — and `SelfReciprocal.lean` already exercises exactly that lemma plus
`rootMultiplicity_eq_zero` in this toolchain. Route (b)'s machinery was therefore *measured
present* while route (a)'s (`rootMultiplicity_eq_one`-flavoured) was still a recollection.

This paid the generality dividend the brief's point 3 predicted: the leaf lemma holds over **any
commutative ring** with `(m : R) ≠ 0`. No `CharZero` is needed until the paper-facing corollary,
and the bracketing argument (`0 < k` and `¬(1 < k)` ⟹ `k = 1`) was never needed — the
factorisation gives order exactly 1 outright. *The stronger statement was again the cheaper one.*

### Mathlib names: all six probed before use, all six real

`rootMultiplicity_mul`, `rootMultiplicity_mul_X_sub_C_pow`, `rootMultiplicity_eq_zero`,
`rootMultiplicity_C`, `rootMultiplicity_zero`, `Finset.cons_induction`, `geom_sum_mul_neg`,
`Nat.cast_ne_zero`. One brief guess was **wrong and caught**: the module path
`Mathlib.Algebra.GeomSum` does not exist (it is `Mathlib.Algebra.Ring.GeomSum`), and the lemma
is `geom_sum_mul_neg`, not the `Polynomial.geom_sum` the brief hypothesised.

## The two carried hypotheses — this is the scope audit

**1. Conditional on gap 3.** The parent node's final sentence, *"the region condition forces
`m_i ≥ 1` for `i ∈ J`"*, is **not proved** — here or anywhere. It is carried as an explicit
named hypothesis `(hm : ∀ i ∈ J, 1 ≤ m i)`. So: **gap 1 is `lean-verified` conditional on gap
3.** This does *not* establish `ord_{t=1} Ψ_T = ν(T)` outright, and the registry node says so.

The hypothesis is witnessed **load-bearing**, not assumed: `m i = 0` gives `1 - X⁰ = 0`, so the
product is `0`, and `rootMultiplicity _ 0 = 0` by Mathlib's convention — the count collapses
from `J.card` to `0`. Executable at `J.card = 2`, `m = (0,3)`.

**2. `s : ℕ`, where the paper has `s(P) ∈ ℤ`.** A negative exponent needs `LaurentPolynomial`,
where `rootMultiplicity` is not the available API. The mathematical content lost is nil (`t^s` is
a unit at `t = 1` either way) but **the Lean statement is strictly weaker than the printed
one**. Recorded, not patched — per the session rule against new mathematics.

**Also not carried:** the reduction `ν(T) = |J|` from tableaux to this index-set cardinality.
That is the sibling node (`ZMod.nu_eq_card_of_content_ones`). This module is pure algebra —
tableau-free, cylinder-free, `h`-free.

## Controls — all four, as integers

**`#print axioms`, all 8 declarations.** Standard three only
(`propext`, `Classical.choice`, `Quot.sound`). Three readings agreed: 8 lines / 8
`depends on axioms` / 8 decl-name matches.

```
'Polynomial.one_sub_X_pow_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.rootMultiplicity_one_sub_X_pow' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.one_sub_X_pow_ne_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.prod_one_sub_X_pow_ne_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.rootMultiplicity_prod_one_sub_X_pow' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.rootMultiplicity_prod_one_sub_X_pow_of_one_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.rootMultiplicity_unit_mul_prod_one_sub_X_pow' depends on axioms: [propext, Classical.choice, Quot.sound]
'Polynomial.rootMultiplicity_prod_one_sub_X_pow_eq_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
```

**Negative control — `sorry` planted at the dependency root** (`one_sub_X_pow_eq`):
**5 of 8 contaminated** with `sorryAx`, **3 spared**.

- contaminated: `one_sub_X_pow_eq`, `rootMultiplicity_one_sub_X_pow`,
  `rootMultiplicity_prod_one_sub_X_pow`, `…_of_one_le`, `…_unit_mul_…`
- **spared:** `one_sub_X_pow_ne_zero`, `prod_one_sub_X_pow_ne_zero`,
  `rootMultiplicity_prod_one_sub_X_pow_eq_zero`

The spared set *is* the audit. The ablation theorem being spared proves it does not quote the
lemma whose hypothesis it witnesses as necessary — the two `≠ 0` lemmas being spared confirms
the `coeff` argument is genuinely independent of the factorisation. Restored
**byte-identical** (`md5 c3f6efb6a5f274ffd35299de61e4ba1c` before and after) and the axiom
reading re-read `diff`-identical to pre-plant, `sorryAx` count **0**.

In the planted arm the pretty-printer **wrapped two records across 4 lines each**: raw
`wc -l` read **14** where there are **8** records. The 10-11 warning fired for real;
`grep -c "depends on axioms"` was the reading that survived it.

**Positive control — an independent mechanism, not a second reading.** The explicit
factorisation `(1-t²)(1-t³) = (t-1)²(t³+2t²+2t+1)`, checked by `ring`, with the cofactor
evaluating to **6 ≠ 0** at `t = 1`. That pins `ord_{t=1} = 2` *without using the theorem*; the
theorem independently returns `J.card = 2`. Two mechanisms, one answer. Cardinalities
**0, 1, 2, 3** with `m ∈ {1,2,3}` all check, plus the `m = (0,3)` ablation returning **0**.
`lake test` exit **0**; `lake build` exit **0**, zero warnings in the new module.

### A positive control that passed for the wrong reason — caught

The first version of the `ring` check **compiled and was worthless**. `TworowD4KernelTests.lean`
does not `open Polynomial`, and `autoImplicit` silently promoted the bare `X` to a
**universally quantified ring variable** — under which `(1-X²)(1-X³) = (X-1)²(X³+2X²+2X+1)` is
a true identity in *any* commutative ring and says nothing whatever about `Polynomial.X`. It was
visible only because the *other* five checks in the same block failed loudly (14 errors) on
`Matrix.cons` index synthesis; had I written the `ring` line alone, it would have gone green and
I would have reported an independent verification that verified nothing.

**A control that typechecks has not thereby been read.** The generic form is the dangerous one:
it is *more* likely to compile, because it asks for less. Fixed by `open Polynomial`, with the
trap written into the test file's comment so the next reader cannot re-lay it.

## Instrument fault: the double-prefix, 8th instance

`python3 code/registry_validate.py proofs/registry/cylindric-statistic.json` — the form printed
in **this session's own prompt** — reported **56 problems**, one per node, every one a phantom
`file … not found`. Diagnosed by reading the source rather than guessing:
`registry_validate.py:456` defaults `proofs_dir` to `dirname(dirname(abspath(registry)))` =
`projects/proofs`, then line 299 joins the node's own `proofs/…` prefix → `proofs/proofs/…`.

The tell was an **identity, not a discrepancy**: *problem count 56 == node count 56*. A fault
that flags every node is reporting on itself.

Correct form: `--proofs-dir .` from `projects/` → **0 problems, exit 0**.
`trustcheck.py … --files-dir .` → **valid, exit 0**.

## Registry

Child `ord-at-one-equals-card-formalised` added under `cyclic-ascent-lemma`, `trust:
lean-verified`, all 8 declarations in `lean`, `sources: ["1110.6356"]`. Nodes **55 → 56**.

**Parent not promoted — the third consecutive decline, and the arithmetic says so.**
`cyclic-ascent-lemma` asserts four things; this module carries **one** of them, and carries it
*conditionally on the fourth*. It stays at `proved`.

First validation run surfaced **1 real violation** (`missing key 'children'`) — mine, fixed.

## Citations carried into the Lean source

Module docstring cites `arXiv:1110.6356` Def. 5.5 (Korff's cylindric weight `Ψ`, the external
object), the `.tex` file and lemma label `lem:kor-is-psi`, and states explicitly that `Ψ` is
**refuted** as an admissible weight by Q410/Q417/Q424 — this file formalises the *engine of the
refutation*, not a claim that `Ψ` is admissible.
