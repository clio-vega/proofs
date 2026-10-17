# For Robin — LEAN 2026-10-17 c2: the collapse multiplier is certified, the mechanism is not

Sorry-free. 22 declarations, 0 `sorry`, standard three axioms only.
`clio-vega/tworow-d4-kernel` commit `1f2b4de`, pushed:
https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/CollapseMultiplier.lean

Full snapshot note:
https://github.com/clio-vega/proofs/blob/main/2026-10-17-lean-collapse-multiplier.md

Registry node (two new lean-verified children):
https://github.com/clio-vega/proofs/blob/main/registry/motegi-typeBC-multiplication.json

**And a correction to my own protocol, found while writing this.** My brief told me to put
this note in `memory/for-robin/` rather than `projects/for-robin/`, to keep the inventory
in one place. I did — and then discovered that **`memory/` is not a git repository**, so a
note that lives only there can never reach you. The directory that is *canonical for
inventory* and the directory that is *reachable by its addressee* are different
directories, and the brief only reasoned about the first. A copy therefore also sits in the
`proofs` repo as `2026-10-17-for-robin-lean-collapse-multiplier.md`, which is where the
previous eleven for-robin notes actually live and the only version you can read. ***A note
filed in the canonical place and nowhere else is invisible to the person it is addressed
to.***

## What is now machine-checked

Over any field, `ρ = (1-q)^2/((1-q t)(1-q t⁻¹))` with `t = w_a^h`:

- `ρ - 1 = q t (1 - t⁻¹)^2 / D` — `prop:rhoXminus1`;
- `ρ = 1 ⟺ q = 0 ∨ t = 1`, hence in `w` the locus is `{q=0} ∪ μ_h`;
- `ρ^X - 1 = q w^h ([K^X(w)]_0^1)^2 / D` — `thm:dict`, with `[K^X(u)]_0^1 = 1 - u^{-h}`
  transcribed from arXiv:2610.15344 L332–345/L352–375;
- the explicit numerators `q(w-w⁻¹)^2` (type C) and `q(w+w⁻¹-2)` (type B);
- the disagreement set: at `w_a = -1`, `ρ^C = 1` and `ρ^B ≠ 1`.

So the type dependence of the collapse locus is **entirely** the `h`-dependence of `μ_h`.
That is the cleanest statement of yesterday's retraction, and it is now certified.

## Three things I want you to see, because two of them are against me

**1. `thm:dict` is definitionally free in Lean, and I said so in the docstring.**
`rho_sub_one_eq_Kentry_sq` is `rho_sub_one` at `t = w^h` with no further work — once
`Kentry` unfolds, the two statements are the same expression. The content of the
dictionary is the *coincidence*, and the coincidence carries information only because I
transcribed `Kentry` from Motegi's source (2610.15344's triangularity, not 1703.07924's
transposed one) instead of reverse-engineering it from the numerator. A reader who sees
`lean-verified` on that node should not conclude more than that.

**2. I did NOT grade the parent nodes `lean-verified`, on purpose.**
`boundary-entry-dictionary` carries a named mechanism gap — `ρ`'s numerator `(1-q)^2`
comes from the bulk factor and contains no `K`-matrix at all; the boundary entry appears
only after subtracting `1` and restricting to `H`. Lean certifies the identity, not the
causality. Promoting that node would have converted a gap I named yesterday into a
certificate. Same for `multiplier-minus-one`, whose statement bundles three branching-formula
consequences that are not formalised. Two new *child* nodes carry the Lean fields, each with
an explicit SCOPE clause.

**3. The brief put the characteristic hypothesis on the wrong sentence, and I'd have
believed it.** My wake brief said `(2 : K) ≠ 0` was real content in the type-`C` locus
statement. It is not: `w^2 = 1 ↔ w = 1 ∨ w = -1` needs no characteristic hypothesis at
all, because in char 2 the two disjuncts simply coincide and the equivalence stays true.
Where the hypothesis is real is the *separation* — if `2 = 0` then `-1 = 1`, so the
type-`B` multiplier collapses at `w_a = -1` too, both loci coincide, and the whole
discriminating measurement carries zero information. **A characteristic hypothesis belongs
to the theorem that separates two loci, not to the theorem that describes one of them.**

A smaller twin inside that: I wanted to refute the char-2 lesion with a concrete field and
found that **in `𝔽₂` the hypotheses `q ≠ 0` and `1 - q ≠ 0` are jointly unsatisfiable** —
the smallest witness would have needed `𝔽₄`. So I proved it in generality instead. A
counterexample over a field too small to satisfy the hypotheses is not a counterexample,
it is an empty domain wearing one's clothes.

## One free theorem

The brief listed `1 - q t⁻¹ ≠ 0` and `h ≠ 0` as hypotheses to lesion. Both are
unnecessary. `h ≠ 0` never entered the statement (`h = 0` forces `t = 1`, both sides hold).
And `1 - q t⁻¹ ≠ 0` is *removable from the level set* given the other two — if it vanishes,
`ρ = 0 ≠ 1` and the right-hand side is false too, so the equivalence holds vacuously.
`rho_eq_one_iff_of_hd1`. That also explains why the `1 - q t ≠ 0` counterexample has *both*
denominators vanishing: the two bad loci meet at `q = t = 1` and cannot be separated.

The other four hypotheses are irremovable, each with a compiled `ℚ` counterexample rather
than a build failure — a deleted binder fails with `Unknown identifier`, which grades the
script and not the statement.

## Method note, and one standing fault

Guards that fired this session: the audit list was regenerated from the file by regex (so a
declaration cannot omit itself); records were counted **wrap-aware**, which mattered — 19
records occupied 28 raw lines, so a line count would have been wrong; a planted `sorryAx`
confirmed the scan was not blind; and the olean was deleted and its mtime confirmed to move
(`…899 → …919`), so `lake` exit 0 was an elaboration and not a trace replay.

Still unfixed, and it is a **generator fault**, not an instance: the LEAN protocol template
writes session notes to `projects/for-robin/` while the canonical directory is
`memory/for-robin/`. Notes now exist in both and an inventory of either silently loses the
rest. I complied by hand today. Fixing the template is the actual repair.
