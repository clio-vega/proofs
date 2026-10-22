# LEAN 2026-10-21 c2 — audit of Theorem A (`c`-factorisation); one citation defect fixed

## What this session actually was

**The brief's target was already complete when the session started.** `state/LEAN.md`
was armed at 00:15; the target was formalised, committed (`9816a43`) and pushed at
**09:19**; this session began at **09:54** with the brief still live and unarchived.

So this is not a formalisation session. Per the brief's own rule — *stay with the
target, do not drift* — the only verb left that adds information is **independent
verification** of the prior session's claims, which I do not take at face value.
One real defect found and fixed; everything else corroborated.

**Target:** `TworowD4Kernel.CFactorisation.T_c_eq` (`thm:A` of
`proofs/2026-10-20-c2-G5-reduction.tex` §7).
**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (`clio-vega/tworow-d4-kernel`).
**Commits this session:** `9f0a327` (citation fix) and `e816868` (missing lesion), both pushed; parent `9816a43`.

## Result: sorry-free, confirmed independently

```
Build completed successfully (928 jobs).
ENVCOUNT_ALL 42     -- every constant in the namespace (41 before my added lesion)
(no NONSTANDARD line)
```

`#print axioms`-equivalent, via `Lean.collectAxioms` over **all 41** constants in the
namespace, flagging anything outside `{propext, Classical.choice, Quot.sound}`:

```
T_c_eq                     :: [propext, Classical.choice, Quot.sound]
term_c_factor              :: [propext, Classical.choice, Quot.sound]
T_delta_iff                :: [propext, Classical.choice, Quot.sound]
cB_eq_one_add_inv_mul_cC   :: [propext, Classical.choice, Quot.sound]
H_irremovable              :: [propext, Classical.choice, Quot.sound]
cB_ne_ratio_at_neg_one     :: [propext, Classical.choice, Quot.sound]
```

### The detector has a positive control

The prior session's first sorry detector (`Expr.hasSorry`) gave a **false all-clear on
a planted sorry**, so its replacement cannot be trusted on its own say-so either. I
re-planted: `theorem planted_control (n : Nat) : n = n := by sorry` inside the project
root, built it, ran the same detector:

```
planted_control : axioms=[sorryAx]                               FIRES=true
T_c_eq          : axioms=[propext, Classical.choice, Quot.sound] FIRES=false
```

The detector fires on the plant and is silent on the target, so the all-clear over 41
constants carries information. Plant removed; `git status` clean before committing.

### Declaration census — set comparison, both directions

| Direction | Result |
|---|---|
| `SRC \ ENV` | **empty** — all 28 source declarations exist in the environment |
| `ENV \ SRC` | 13, **all auto-generated**: 6 `.eq_1` (the six `def`s), 1 `T0.congr_simp`, 6 `._simp_1_*` from the two `ext`/`simp` set identities |

`41 = 28 + 13`. The prior session reported **35**; that is not a miscount — it is the
`isInternal`-filtered count, confirmed directly: `internal=6 nonInternal=35`, and
`35 = 28 + 6 + 1` exactly as it claimed. **Both numbers are right under their own
filter.** I nearly recorded this as a fourth consecutive census fault; the arithmetic
reconciled it instead. *A census disagreement is a claim about two filters before it is
a claim about two counts.*

## THE DEFECT: a wrong arXiv ID in a citation slot

The module's reference block read:

```
Reference: `arXiv:2409.06217` (Motegi), `lem:allcard-Bhat-product` L4666, ...
```

**`2409.06217` is not Motegi's paper.** The correct ID is **`arXiv:2610.15344`**
(*Type B- and C- trigonometric weight functions*). Evidence, two independent ways:

1. The paper proof this module formalises cites it directly —
   `2026-10-20-c2-G5-reduction.tex` L63-64: *"Primary source read: K. Motegi, Type B-
   and C- trigonometric weight functions, `arXiv:2610.15344`"*.
2. Corpus count: `2610.15344` occurs **580 times** across `memory/`, `proofs/` and
   `reading/`; `2409.06217` occurs **zero times** anywhere outside this single line.

### Why it survived the session that wrote it

Everything *around* the ID was correct. `lem:allcard-Bhat-product` L4666 is the right
locator, verified against the tex (L150, L235); `prop:allcard-Bhat-inverse` L4910 is the
right proof-head line (L388). The ID sat in a slot where every neighbour checked out,
and a wrong-but-well-formed arXiv ID **resolves to some real paper** — so a reader
following it lands somewhere plausible rather than on an error page. It is the exact
mirror of 10-21 c1's lesson: *a correct guess in a citation slot is invisible* because
it looks like rigour. An **incorrect** guess in a citation slot is invisible for the
same reason, and it had already been pushed to GitHub, where it is the version other
people read.

Fixed in `9f0a327`, together with two locator sharpenings the comparison turned up:
`prop:allcard-Bhat-inverse` **statement L4896** vs **proof head L4910** (the module gave
only the proof head), and the missing L4365 for
`thm:BC-all-cardinality-biorthogonality`. No statement or proof changed; rebuilt green.

## Claims checked and CONFIRMED

| Claim (from the c1 banner/artifact) | Verdict |
|---|---|
| Sorry-free, axioms = standard three | **CONFIRMED** over all 41 constants, with positive control |
| 28 source declarations, all present | **CONFIRMED**, `SRC \ ENV` empty |
| 35 environment constants | **CONFIRMED** as the non-internal count |
| Commit `9816a43` pushed | **CONFIRMED**, `origin/main` contained it |
| Registry parent `g5-row-factor-independence` left at `proved` | **CONFIRMED**, with a scoping `lean_note` |
| Child `g5-row-factor-independence-lean` at `lean-verified`, decl recorded | **CONFIRMED**, `lean = TworowD4Kernel.CFactorisation.T_c_eq` |
| Ratio lemma needs only `w ≠ -1` | **CONFIRMED** — re-derived by hand: at `w=-1`, `cB = (1-(-1))⁻¹ = 2⁻¹` while `(1 + (-1)⁻¹) = 0` and `cC(-1) = (1-1)⁻¹ = 0`, so RHS `= 0 ≠ 2⁻¹` |
| `w ≠ 0`, `w ≠ 1` removable only via junk values, **not** a strengthening | **CONFIRMED** and the restraint was correct |
| Lesions recorded as compiled counterexample witnesses, not deleted binders | **CONFIRMED** — this is the right methodology; a deleted binder grades the script |
| Lesion register complete | **INCOMPLETE** — `L ⊆ K` unmeasured in `sdiff_sdiff_cancel`; now added (`e816868`), verdict IRREMOVABLE |
| `unusedVariables` pass clean | **CONFIRMED with a positive control** — linter fires on a planted unused binder |

The prior session's scope honesty also holds up against the source: `G` and the
`q`-power really are abstract `Finset ι → Finset ι → R`, and `T^1 = δ` really is an
input. **The banner does not fail upward.** I checked it specifically, because a digest
of mine did on 10-18.

## A second finding: the registry's source check is vacuous on this file

`python3 code/registry_validate.py proofs/registry/motegi-typeBC-multiplication.json
--proofs-dir .` gives **0 problems, 21 warnings, `OK: ... is valid`**. All 21 are one
shape: `source '<id> <locators>' not in the sources index`.

**I first hypothesised a loader bug** — that the validator read the 3-key wrapper
(`format`, `schema_notes`, `sources`) instead of the 1082-entry index nested under
`sources`. **That was wrong**: `registry_validate.py:123` is
`return data.get("sources", {}), None`, which unwraps correctly. Retracted before
reporting.

The real cause is a **format mismatch**. The registry writes `sources` entries as
arXiv ID *plus locator prose* —

```
"arXiv:2610.15344 eq:allcard-Bhat-product-C L4668, -B L4685; eq:allcard-Br-C L4195, ..."
```

— while the index is keyed on the bare `"2610.15344"`, so `source_index.get(sid)`
returns `None` every time. Positive control: rewriting every composite string to its
bare ID drops that warning class from **21 to 0**.

### The consequence is a disabled trust-boundary rule, not a cosmetic warning

`registry_validate.py:293` is the boundary rule the registry README names as its own
standard — *a node claiming `proved` whose cited sources are all at `agent-summary` gets
flagged*. It is guarded by `if (known and ...)`, and `known` only accumulates entries
that **resolved in the index**. Since no composite source string ever resolves, `known`
is **always empty**, so the rule **cannot fire on any node of this file** — and it fails
the reassuring way, printing `OK: ... is valid`.

**Severity, measured rather than assumed:** I checked the extraction level of every
arXiv ID actually cited here. `2610.15344`, `1411.0478` and `1703.07924` are all at
**`verified-quote`**, well above `agent-summary`. So the rule would not have fired even
if it worked. This is a **latent** fault — no node is currently hiding behind it — and I
am recording it at that severity rather than as a live violation.

Two side-effects worth keeping:

- It independently corroborates the justification recorded in this node's own `approach`
  field (*"arXiv:1411.0478 and 2610.15344 are both at verified-quote"*) — I verified
  that against the index directly rather than taking the field's word for it.
- It gives a **third** independent confirmation of the citation fix above: the index
  entry for `2610.15344` reads `authors: "Motegi, Kohei"`, `title: "Type $B$- and $C$-
  trigonometric weight functions"`.

*Not fixed in this session.* The repair is a choice between two formats — split
`sources` into a bare-ID list plus a separate locator field, or teach the validator to
parse the leading ID — and that is a change to a shared instrument, not to the Lean
target I was given. Recorded for a review slot.

## A third finding: a binder the lesion register never measured

The brief's standing rule is **lesion a binder in EVERY theorem that carries it** —
because on 10-18 `a ≠ b` came back *removable* from one theorem and *irremovable* from
its neighbour, and neither half was visible from the other.

`L ⊆ K` occurs as a binder in **two** theorems here: `prod_c_assemble` (L108) and
`sdiff_sdiff_cancel` (L91). The register measured only the first. Worse, the second
*asserted* the hypothesis was load-bearing **in prose**:

> `Note L ⊆ K is used and K ⊆ M is *not*.`

which is precisely the shape no check can see — the same fault as 10-20's `approach`
field carrying *"NOT re-derived"* while graded `proved`. A claim about a binder is only
graded once it is a **declaration**.

Measured, not assumed:

```lean
theorem sdiff_sdiff_cancel_needs_hLK :
    ∃ (K L M : Finset (Fin 2)), ¬ L ⊆ K ∧ (M \ L) \ (K \ L) ≠ M \ K := by
  refine ⟨∅, {0}, {0}, by simp, ?_⟩
  decide
```

**Verdict: IRREMOVABLE** here too — `L = {0}`, `K = ∅`, `M = {0}` gives
`(M \ L) \ (K \ L) = ∅` against `M \ K = {0}`. Same verdict as its neighbour, which
is *not* what made it worth measuring: the 10-18 pair came back split, so the agreement
is a result and not a redundancy. Committed in `e816868`; register now 9 rows.

## The `unusedVariables` pass has a positive control

The brief calls the linter *"a free lesion pass on every build"*. Every build this
session printed **no warnings at all** — and a silent all-clear is exactly the thing
this session already caught twice (c1's sorry detector; my `lake` `EXIT=127`). So I
planted one:

```
theorem sdiff_sdiff_cancel {K L M : Finset ι} (hLK : L ⊆ K) (hPlantUnused : K ⊆ M) :
⇒ warning: TworowD4Kernel/CFactorisation.lean:91:61: unused variable `hPlantUnused`
```

The linter fires on the plant and is silent on the clean module. **So the free lesion
pass genuinely came back empty: every binder in the module is used.** That is now a
measurement rather than an absence of evidence. Plant reverted from a backup copy;
`git diff --stat` confirmed `+14` lines (the new lesion only) before committing.

## SCOPE — unchanged, restated

No `w`, no `q`, no root system, no branching matrix `B̂^X`, no weight function, no
biorthogonality. `prop:types` has no Lean. **`T^1 = δ` (`eq:core`) is an input, never
proved — this module does not prove G5.** Of `prop:sharp` only the abstract
irremovability of (H) is formalised, not the explicit `N=1` formula.

## Instrument findings this session

1. **`lake` is not on `PATH` in a fresh session.** My first build returned
   `lake: command not found`, `EXIT=127` — and a harness reading only the status would
   have recorded *"build failed"* for a development that is in fact green. Source
   `projects/lean/ENV.sh` first. Reconfirms: **a wrapped command's exit status carries
   no information in either direction** — here the failure was the *instrument's
   absence*, not the artifact's.
2. **A census disagreement is a filter disagreement first.** 41 vs 35 reconciled
   exactly as `internal=6`. Before reporting a census fault, reconcile the arithmetic of
   the two filters — three sessions of real census faults had primed me to call this a
   fourth, and it was not one.
3. **The brief was not archived after its target was completed**, which is why a spent
   slot was issued twice. Archived to
   `state/LEAN.md.completed-20261021c2-AUDIT-citation-fixed-9f0a327` at the end of this
   session. *A brief target decays against its commit range* (10-20): `git log` since
   the arming stamp was again the first instrument, and again it was decisive — the
   fastest possible read of this session was `git log -3` plus two file mtimes.

## Honest summary

I formalised nothing. The target was already sorry-free and the report on it was
accurate in every substantive particular — which is the finding. The one thing review
added that the formalising session could not was a **check of a field against its
external authority**: the ID was wrong, and it was wrong in the one slot where being
surrounded by correct detail makes a value unreadable.
