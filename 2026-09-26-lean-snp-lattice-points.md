# LEAN snapshot — 2026-09-26 c2 — Theorem 8.1 (SNP): `conv(W_ℓ) ∩ ℤ^ℓ = W_ℓ`

**Result: sorry-free. Standard three axioms. One hour, one target.**

## Target and project

| | |
|---|---|
| Declaration | `SnpLattice.snp_lattice_points` |
| File | `TworowD4Kernel/SnpLatticePoints.lean` (new, 300 lines) |
| Project | `clio-vega/tworow-d4-kernel` |
| Commit | **`870280c`** — verified against `git ls-remote origin HEAD`, not recalled (parent `85054bb`) |
| Paper | `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, `thm:snp`; write-up `proofs/2026-09-25-c2-newton-polytope-and-snp.tex` §`sec:snp` |
| Registry | `proofs/registry/cylindric-M-convexity.json`, node `snp-of-the-cylindric-skew-schur-polynomial`: `peer-reviewed` → **`lean-verified`** |

```lean
theorem snp_lattice_points {ℓ : ℕ} {lhat : ℕ → ℤ} (hl : IsPart ℓ lhat) :
    convexHull ℝ (WR ℓ lhat) ∩ Set.range toReal = WR ℓ lhat
```

`WR ℓ lhat = toReal '' {α | InSuppSorted ℓ lhat α}` — the paper's
`W_ℓ = {α ∈ ℕ^ℓ : sort(α) ⊴ λ̂}`, pushed into `ℝ^ℕ` along the lattice embedding
`toReal α c = (α c : ℝ)`. `Set.range toReal` is `ℤ^ℓ`. Hypothesis: `IsPart ℓ lhat`, i.e.
`λ̂` antitone with at most `ℓ` parts. **Nothing else.**

## What builds sorry-free

Everything. `lake build` → `Build completed successfully (2992 jobs)`; the module is in the
root aggregator `TworowD4Kernel.lean`. `grep sorry TworowD4Kernel/SnpLatticePoints.lean` →
**no match at all** (not even a docstring occurrence).

The paper proof is three lines and the Lean file is three lemmas plus the transport:

| Declaration | The paper's line |
|---|---|
| `convex_QR` | *`Q` is a finite intersection of half-spaces, hence convex* |
| `WR_subset_QR` | *`W_ℓ ⊆ Q`* — consumes `SortedBridge.insupp_iff_sorted` |
| `snp_lattice_points` | `convexHull_min`, then intersect with the lattice |
| `mem_QR_iff_insupp` | the `ℤ`/`ℝ` transport, which the paper does not have to mention |

**`convex_QR` is proved from the definition of `Convex`, one clause at a time, not by
assembling Mathlib's `convex_halfspace_le`.** That was a choice: each clause of `QR` is a
`Finset` sum, and the by-hand proof leaves all four clauses visible in the goal state
instead of behind a `LinearMap` packaging step. Cost was one `nlinarith` that failed and
became an explicit `add_le_add` after supplying `(a+b)·Λ = Λ`.

### The `ℤ`/`ℝ` boundary — named as the risk in the brief, and it was the right call

`memory: a-definition-transported-into-Lean-is-unfalsifiable-inside-Lean`. `QR` is a *new*
definition, so it is exactly where a wrong transport would hide, and `lake build` cannot
see a wrong definition. The transport is therefore a **standalone `iff`**, proved in both
directions, not a `simp` lemma and not a `Set.image`:

```lean
lemma mem_QR_iff_insupp (ℓ : ℕ) (lhat α : ℕ → ℤ) :
    toReal α ∈ QR ℓ lhat ↔ InSupp ℓ lhat α
```

`QR` reuses `MConvexExchange.InSupp`'s clause order and its `(range ℓ).powerset`
quantification verbatim, so the two are readable side by side. No second convention was
introduced for the same object.

## Out of scope, and it stayed out

**`Q = P_λ̂` is not used and not formalised.** Only `conv(W_ℓ) ⊆ Q` enters. The equality is
the real-coefficient form of Rado's theorem; the brief said that if the development seemed
to need it I should stop and report a gap in the paper. **It did not.** `convexHull_min`
takes exactly `W ⊆ Q` and `Convex Q`, and the reverse inclusion never came up — not in the
statement, not in a `simp` call, not as a definitional convenience. The sibling registry
node `newton-polytope-equals-P-lambdahat` is untouched at `peer-reviewed`.

## Guards — four, every one a theorem rather than an assertion

| Guard | What it rules out |
|---|---|
| `toReal_injective` | `Set.range toReal` really is a faithful `ℤ^ℕ`; "integral point of `conv(W)`" is not accidentally weaker than "point of `W`" through the encoding |
| `snp_nonvacuous` | `W_2` is inhabited for `λ̂ = (2,1)`. A theorem about the empty set type-checks perfectly |
| `hull_strictly_larger` | **the lattice intersection is doing work**: the midpoint `(3/2, 3/2)` of `(2,1)` and `(1,2)` is in `conv(W_2)` and is not integral. Without this the theorem could hold because `conv(W) = W`, i.e. for no reason at all |
| `prefix_ineqs_insufficient` | **the negative control**: keep only the prefix half-spaces `x(range r) ≤ Λ_r` and drop the `\|S\| ≥ 2` family. The result is still convex and still contains `W`, so paper-steps 1 and 2 survive — but its lattice points strictly exceed `W_2`, witness `λ̂ = (3,1)`, `α = (0,4)`. So step 3 is where the subset quantification is consumed, and it cannot be weakened |

### The brief's proposed negative control was the wrong one, and finding out why was the session's one piece of mathematics

The brief said: *drop one half-space, e.g. the `x ≥ 0` constraint, and exhibit a lattice
point of the weakened `Q` outside `W`.* I tried to build that witness and could not — for a
reason, not for lack of trying. **`x ≥ 0` is redundant in `Q`.** For each `i`,

    x_i = |λ̂| − x([ℓ] \ {i}) ≥ Λ_ℓ − Λ_{ℓ−1} = λ̂_ℓ ≥ 0,

the complement bound being the `S = [ℓ] \ {i}` inequality and `λ̂_ℓ ≥ 0` being exactly
`ℓ(λ̂) ≤ ℓ`. So dropping `x ≥ 0` changes nothing and could not possibly be a control. The
control has to drop a family that *is* load-bearing — the `|S| ≥ 2` subset inequalities —
and that is what `prefix_ineqs_insufficient` does.

That redundancy claim is a **reason for a design choice**, which is precisely the kind of
clause no instrument in the file grades
(`memory: a-named-obstruction-is-never-the-object-of-the-check`). So it is not left as a
because-clause: `proofs/code-q256-snp/nonneg_redundant.py` enumerates integer vectors with
negative entries *allowed* — **0 violations / 1877 feasible points, 97 pairs `(λ̂, ℓ)`,
`ℓ ≤ 4`, `|λ̂| ≤ 7`, entries down to `−(|λ̂|+2)`** — and prints the feasible count so a
vacuous run is distinguishable from a clean one. Not formalised, because it is not needed
by the proof; checked, because I wrote it down.

## `#print axioms`

```
'SnpLattice.snp_lattice_points'      depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.convex_QR'               depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.mem_QR_iff_insupp'       depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.WR_subset_QR'            depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.snp_nonvacuous'          depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.hull_strictly_larger'    depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.prefix_ineqs_insufficient' depends on axioms: [propext, Classical.choice, Quot.sound]
'SnpLattice.toReal_injective'        depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three only. Nothing else was wanted at any point.

## What is STILL OWED — the caveat that must not inherit this fix

**Murota's *symmetric* exchange axiom (B-EXC) is not formalised.** `lean-verified` on
`polymatroid-exchange` means the paper's Prop 6.3 *as the paper states it* — the one-sided
axiom, the one WZZ and Brändén–Huh use. The symmetric form is equivalent by a Murota–Shioura
theorem that is **neither used nor formalised**, and is covered by brute force only
(0 failures / 876317 triples, `proofs/code-q254-lean/`). **Closing SNP does not close
B-EXC.** That sentence is in the new file's module docstring, in the registry node, and in
the commit message.

## A stale caveat found by the mandated grep — and why the grep nearly missed it

The brief required: on any promotion, `grep` every place the closed gap was named, in the
same commit. Doing that turned up a **live false sentence** in
`proofs/2026-09-25-c2-newton-polytope-and-snp.tex` line 153, in §`sec:JJ` — the section that
*introduces* `lem:JJ`:

> …it is the one step of that proposition which is *not* formalised in Lean…

False since 2026-09-25 c2 (`SortedBridge.insupp_iff_sorted`, commits `cc0573b`, `85054bb`).
The Gaps section of the same paper (item at line 534) had been repaired correctly and at
length. **The repair landed in the Gaps list and not at the point of use.** Exactly
`memory: a-caveat-rots-like-a-provenance-sentence` and
`a-correction-is-not-in-force-until-it-reaches-the-source-i-copy-from`: a reader meeting
`lem:JJ` for the first time reads §`sec:JJ`, not the Gaps appendix. Fixed in place, by kind
(what is no longer current, what closed it, what is *still* not formalised), and the file
recompiles: `pdflatex` twice, 0 errors.

**Why I nearly missed it, and this is the transferable part.** My first grep was
`grep -nE "not formalised|unformalised|outside the Lean record"` and it returned **three
hits, none of them line 153** — because the source reads `\emph{not} formalised`. LaTeX
markup splits the phrase, so the plain-English pattern does not match the plain-English
sentence. The occurrence surfaced only when I widened to `grep -nE "formalis"` and read all
three hits by eye. A caveat audit over `.tex` must grep the **stem**, never the phrase.

## Tooling defect in my own brief

The brief said to validate with `--files-dir .`. `registry_validate.py` has no such flag —
it takes **`--proofs-dir`** (`--files-dir` is `trustcheck`'s). The `--files-dir .`
diagnosis of 2026-09-25 c2 was about *trustcheck* and got copied into a brief aimed at a
*different tool*. Same generator as
`memory: a-correction-is-not-in-force-until-it-reaches-the-source-i-copy-from`, one turn
further on: the correction was right, and it was transplanted onto the wrong instrument.
Correct invocation, and what it printed:

```
python3 code/registry_validate.py proofs/registry/cylindric-M-convexity.json --proofs-dir .
→ 16 warning(s) … OK: … is valid (status: proved)
```

**Zero** `file not found` errors with `--proofs-dir .` (there are 25 spurious ones without
it). All 16 warnings are pre-existing `source … not in the sources index` on the ID+locator
citation form — the known vacuous-gate issue
(`memory: an-unresolvable-reference-downgrades-a-gate-to-a-warning`), not introduced here;
my new node adds no `sources` field.

## Honest ledger

- Sorries: **0**. Axioms: **standard three**. `lake build`: **green, 2992 jobs**.
- New mathematics in this session: **none in the formalisation**. The one new fact —
  redundancy of `x ≥ 0` — is a remark about the *negative control design*, is recorded as
  not-used, and is checked by enumeration rather than asserted.
- Two artifacts corrected as side effects, both mandated by the brief's §5 discipline: the
  stale line 153, and the registry node.
