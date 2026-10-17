# LEAN 2026-10-17 c2 — the collapse multiplier's level set, and `thm:dict`

**Target** (from `/home/clio/state/LEAN.md`, armed by WAKE 2026-10-17 c2): the algebraic
core of Proposition `prop:rhoXminus1` and Theorem `thm:dict` of
`proofs/2026-10-17-typeB-multiplier-and-boundary-entry.tex` (commit `e347ed3`, L410–L455).

**Project** `/home/clio/projects/lean/tworow_d4_kernel` (`clio-vega/tworow-d4-kernel`),
new file `TworowD4Kernel/CollapseMultiplier.lean`, imported from the root module.
Commit `1f2b4de`, pushed; verified as a set difference (`git log --oneline origin/main..HEAD`
is empty, not by reading a "pushed" line).

## Verdict

**Sorry-free, all four parts of the brief landed, plus `thm:dict` proper.**
The brief asked for (1) the identity, (2) the level set, (3) the `rootsOfUnity` form,
(4) the `h = 1` / `h = 2` table, and said that if it finished early I should extend. It
finished, so the file also carries the boundary-entry half of `thm:dict`
(`Kentry`, `rho_sub_one_eq_Kentry_sq`, `numerator_typeB`, `numerator_typeC`) and the
disagreement-set theorem.

## Counts, separately as required

| quantity | value |
|---|---|
| declarations | **22** |
| `theorem` | 20 |
| `noncomputable def` | 2 (`rho`, `Kentry`) |
| `lemma` | 0 |
| `example` | **0** — nothing here is an `example`, so nothing is undeliverable |
| `sorry` | **0** |
| `#print axioms` records | **23** = 22 declarations + 1 planted `sorryAx` |
| target declarations with a nonstandard axiom | **0** |

`#print axioms` for every one of the 22 reads
`[propext, Classical.choice, Quot.sound]` — the standard three, nothing else.

### How the audit was made trustworthy

Three guards, each installed because a previous session was fooled without it.

1. **The audit list was regenerated from the source file by regex**, not hand-written, so a
   declaration cannot be silently omitted from its own audit.
2. **Records counted wrap-aware.** `#print axioms` line-wraps long names: in the first run
   19 records occupied 28 raw lines. A line count would have been wrong. I count
   occurrences of `depends on axioms: [...]` parsed across newlines.
3. **Positive control.** `theorem PLANT_sorryAx : (0 : Nat) = 1 := by sorry` was planted in
   the audit file. It **fired** (`[sorryAx]`), so the scan is not blind. Without this,
   "0 nonstandard records" and "the scan did not run" grep identically.

### The build is an elaboration, not a replay

`lake build` exit 0 can be a trace-hash match. I **deleted the olean** and rebuilt; its
mtime moved, `1792273899 → 1792273919`. Full project build: 3469 jobs, green.

## The declarations

Over any field `K`, with `rho q t = (1-q)^2 / ((1 - q*t) * (1 - q*t⁻¹))` and `t = w_a^h`
(`h = 1` short roots / type `B`, `h = 2` long roots / type `C`):

- `rho_sub_one` — `prop:rhoXminus1`: `ρ - 1 = q t (1 - t⁻¹)^2 / D`.
- `rho_eq_one_iff` — the level set: `ρ = 1 ↔ q = 0 ∨ t = 1`.
- `Kentry h u = 1 - (u^h)⁻¹` — Motegi's charge-violating entry `[K^X(u)]_0^1`
  (arXiv:2610.15344, `eq:Ktri-components` L332–345, `eq:KC/KB-components` L352–375):
  `1 - u⁻²` in type `C`, `1 - u⁻¹` in type `B`.
- `rho_sub_one_eq_Kentry_sq` — `thm:dict`: `ρ^X - 1 = q w^h ([K^X(w)]_0^1)^2 / D`.
- `numerator_typeC` — `q w^2 (Kentry 2 w)^2 = q (w - w⁻¹)^2`.
- `numerator_typeB` — `q w (Kentry 1 w)^2 = q (w + w⁻¹ - 2)`.
- `rho_eq_one_iff_pow`, `rho_eq_one_iff_mem_rootsOfUnity` — the locus in `w` is
  `{q = 0} ∪ μ_h`, the second via Mathlib's `mem_rootsOfUnity'`.
- `rho_eq_one_iff_typeB`, `rho_eq_one_iff_typeC` — the paper's table.
- `typeC_collapses_typeB_does_not_at_neg_one` — the disagreement set.
- `rho_zero`, `rho_one` — the two values on the level set.
- seven lesion declarations, below.

### An honest word about what `thm:dict` costs in Lean

`rho_sub_one_eq_Kentry_sq` is `rho_sub_one` at `t = w^h` **with no further work** — once
`Kentry` is unfolded the two statements are literally the same expression. That is not a
defect of the formalisation; it is what the dictionary *is*. The content is the
**coincidence**, and the coincidence carries information only because `Kentry` was
transcribed from Motegi's source — with 2610.15344's triangularity, not 1703.07924's
transposed one — rather than reverse-engineered from the numerator. I have said this in
the docstring so a later reader cannot mistake a definitional identity for a proved
mechanism.

And the parent registry node's **mechanism gap** is untouched: `ρ`'s numerator `(1-q)^2`
comes from the bulk factor and contains no `K`-matrix; the boundary entry appears only
after subtracting `1` and restricting to `H`. Lean certifies the identity, not the
causality. **That is why I did not set the parent node to `lean-verified`** — doing so
would convert a named gap into a certificate. Two new child nodes carry the Lean field
instead (see *Registry*).

## Lesions — the honest part

The brief required deleting each hypothesis once. **Deleting a binder makes the stored
proof fail with `Unknown identifier`, which grades the script and not the statement** —
that is not a measurement. So every lesion below is settled instead by a *compiled
counterexample* (irremovable) or a *re-proof without the hypothesis* (removable).

| hypothesis | verdict | witness |
|---|---|---|
| `t ≠ 0` in `rho_sub_one` | **irremovable** | `rho_sub_one_needs_ht`: `q = 3`, `t = 0` |
| `1 - q t ≠ 0` in `rho_sub_one` | **irremovable** | `rho_sub_one_needs_hd1`: `q = 1/2`, `t = 2` |
| `1 - q t⁻¹ ≠ 0` in `rho_sub_one` | **irremovable** | `rho_sub_one_needs_hd2`: `q = t = 2` |
| `t ≠ 0` in `rho_eq_one_iff` | **irremovable** | `rho_eq_one_iff_needs_ht`: `q = 2`, `t = 0` |
| `1 - q t ≠ 0` in `rho_eq_one_iff` | **irremovable** | `rho_eq_one_iff_needs_hd1`: `q = t = 1` |
| `1 - q t⁻¹ ≠ 0` in `rho_eq_one_iff` | **REMOVABLE** | `rho_eq_one_iff_of_hd1` — a stronger theorem for free |
| `h ≠ 0` | **never needed** | not a hypothesis of `rho_eq_one_iff_pow` at all |
| `(2 : K) ≠ 0` in the separation theorem | **irremovable, real content** | `typeB_also_collapses_at_neg_one_of_char_two` |
| `q ≠ 0` in the separation theorem | **irremovable** | `typeB_also_collapses_at_neg_one_of_q_zero` |

Four of these are worth saying out loud.

**`t ≠ 0` is needed *because Lean's inverse is total*.** At `t = 0` the expression
`rho q 0` evaluates to `(1-q)^2` — a perfectly good element — and it equals `1` at
`q = 2`. So the lesioned statement typechecks and is *false*. On paper `t = 0` is outside
the domain and nobody would write it down; in Lean the domain restriction has to be a
hypothesis or the theorem is wrong. The counterexample is in `ℚ`, by `norm_num`.

**`1 - q t⁻¹ ≠ 0` is removable from the level set, and I can say *why* there is no
isolated counterexample for it.** If that factor vanishes then `ρ = 0 ≠ 1`, so the left
side is false; and the right side is false too — `q = 0` would make the factor `1`, and
`t = 1` would make it equal to `1 - q t`, contradicting `hd1`. So the equivalence holds
vacuously, and the hypothesis was never spent. `rho_eq_one_iff_of_hd1` proves it.
This is also exactly why the `hd1` counterexample at `q = t = 1` has *both* denominators
vanishing: the two loci meet there and cannot be separated.

**`h ≠ 0` was in the brief's list and is not in the theorem.** At `h = 0` we get
`w^0 = 1`, hence `t = 1`, hence `ρ = 1`, and the right-hand side `q = 0 ∨ w^0 = 1` is
true as well: both sides hold unconditionally and the equivalence is fine. Free
strengthening, recorded in the docstring.

**The characteristic hypothesis is spent on the *separation*, not on either locus.**
This is the finding I did not expect. `rho_eq_one_iff_typeC` — `w^2 = 1 ↔ w = 1 ∨ w = -1`
— needs **no** characteristic hypothesis: in characteristic `2` the two disjuncts merely
coincide and the equivalence stays true. The brief predicted `(2 : K) ≠ 0` would be real
content at this statement, and at *that* statement it is not. Where it is real is the
*discrimination*: `typeB_also_collapses_at_neg_one_of_char_two` proves that if
`(2 : K) = 0` then `rho q (-1) = 1` as well, so the type-`B` multiplier collapses at
`w_a = -1` too, the two loci coincide, and the paper's whole discriminating measurement
carries no information. The brief was right that the hypothesis matters and wrong about
which sentence carries it — ***a characteristic hypothesis belongs to the theorem that
separates two loci, not to the theorem that describes one of them.***

I proved this in generality rather than by exhibiting a characteristic-`2` field, and the
reason is itself a small lesson: **in `𝔽₂` the hypotheses `q ≠ 0` and `1 - q ≠ 0` are
jointly unsatisfiable**, so the smallest concrete counterexample would have needed `𝔽₄`.
A counterexample over a field too small to satisfy the hypotheses is not a counterexample;
it is an empty domain wearing one's clothes.

## Mathlib interaction

The brief predicted the session's time would go into the `rootsOfUnity` connection. It
did not — `mem_rootsOfUnity' (k : ℕ) (ζ : Mˣ) : ζ ∈ rootsOfUnity k M ↔ (ζ : M) ^ k = 1`
does the whole job in one `rw`, and `Units.ne_zero` supplies `w ≠ 0` with no coercion
fight. In this Mathlib (`v4.30.0`) `rootsOfUnity` takes `k : ℕ`, not `ℕ+`, which is why
no positivity hypothesis is needed anywhere.

Where the time actually went: **`field_simp` could not finish `rho_sub_one`** because
`1 - q t⁻¹` is `(t - q)/t` and it had no way to see `t - q ≠ 0`. Supplying that as an
explicit `have` (derived from `hd2` and `ht`) closed it. One other small thing:
`Mathlib.Algebra.Order.Ring.Lemmas` no longer exists in this Mathlib and `linear_combination`
is not transitively imported by `Mathlib.RingTheory.RootsOfUnity.Basic`; both were replaced
with `rw [show … by ring]`.

## Registry

`proofs/registry/motegi-typeBC-multiplication.json` (backup
`.bak-1017c2-lean-pre`). Two **new child nodes** under `longroot-twisted-collapse`:

- `multiplier-minus-one-lean` — `trust: lean-verified`,
  `lean: TworowD4Kernel.CollapseMultiplier.rho_sub_one`
- `boundary-entry-dictionary-lean` — `trust: lean-verified`,
  `lean: TworowD4Kernel.CollapseMultiplier.rho_sub_one_eq_Kentry_sq`

**The parents stay where they were, deliberately.** `multiplier-minus-one` stays `proved`
because its statement bundles three consequences — the `ρ(ρ-1)` Laurent discrepancy,
`ord_q = 1` for every long-root residue, and the `q = 0` trivialisation that derives
`q0-gate-vacuous` — which go through the branching formula and are **not** formalised.
`boundary-entry-dictionary` stays `computed` because of its named mechanism gap. Each new
node's `statement` field carries an explicit **SCOPE** clause saying what Lean does and
does not cover, so the grade cannot be read as broader than it is.

`python3 code/registry_validate.py proofs/registry/motegi-typeBC-multiplication.json --proofs-dir .`
→ `OK ... is valid (status: proved)`, 7 advisory warnings. All 7 are one class: a `sources`
string carrying line locators (`arXiv:2610.15344 L332-345, …`) where `sources.json` indexes
the bare id. **Six of the seven predate this session**; the seventh is mine and follows the
convention already established in the file by its own parent node. Not repaired, because
changing one node's format and not the other six would make the file inconsistent in a way
the validator cannot see.

## Where the note went

Per the brief, the companion note goes to **`memory/for-robin/`** (the canonical path), not
`projects/for-robin/`. The LEAN protocol template writes to the second path; there are
notes in both and an inventory of either loses the rest. That is a **generator fault** and
it is still unfixed — today I only complied by hand.
