# LEAN 2026-10-05 c3 — `Lemma boxslice` (Step A of (Q)): the all-partners box-slice exchange

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**New file:** `TworowD4Kernel/BoxSlice.lean` (294 lines)
**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib pinned by `lake-manifest.json`
**Registry:** `proofs/registry/cylindric-lorentzian.json`, new node `Q-stepA-boxslice-lean`
(child of `Q-stepA-bead-support-box-slice`), trust `lean-verified`
**Paper proof:** `proofs/2026-10-04-c3-Q-lorentzian.tex`, `lem:boxslice`
**Source for the route:** arXiv `1902.03719` (the parent node's `sources`)

## Status: sorry-free. Zero sorries, zero new axioms.

`lake build` on the whole project: **exit 0, 3192 jobs, "Build completed successfully."**

## Target declarations, all building

```lean
def InBox {m : ℕ} (P Q : Fin m → ℤ) (σ : ℤ) (y : Fin m → ℤ) : Prop :=
  (∀ c, P c ≤ y c) ∧ (∀ c, y c ≤ Q c) ∧ (∑ c, y c = σ)

theorem boxSlice_all_partners {m : ℕ} (P Q : Fin m → ℤ) (σ : ℤ)
    {y y' : Fin m → ℤ} (hy : InBox P Q σ y) (hy' : InBox P Q σ y')
    {i l : Fin m} (hi : y' i < y i) (hl : y l < y' l) :
    InBox P Q σ (ex i l y)

theorem mConvex_boxSlice {m : ℕ} (P Q : Fin m → ℤ) (σ : ℤ) :
    MConvex {y : Fin m → ℤ | InBox P Q σ y}

def beadSupport (l h : ℤ) (y : Fin 3 → ℤ) : Prop := InBox 0 (beadQ l h) h   -- beadQ l h = ![h, h, h-l]
theorem mConvex_beadSupport (l h : ℤ) : MConvex {y : Fin 3 → ℤ | beadSupport l h y}
theorem beadSupport_all_partners (l h : ℤ) …                                  -- all-partners on the instance
```

Plus, as a second and independent route to each of the two main statements:
`boxSlice_all_partners_via_exBoxSum` (from `SublevelMConvex.ex_box_sum`, saturating
`j := max (kfun y) (kfun y')`) and `mConvex_boxSlice_via_inS` (from
`ParityObstruction.mConvex_inS` at the vacuous `j := ∑ c, negPart (P c)`), with the
laundering lemmas `inS_of_inBox`, `inBox_of_inS`, `kfun_le_of_lower`,
`inBox_iff_inS_saturated`, and the already-proved constant-sum consequence re-exposed as
`sum_eq_of_inBox`.

## `#print axioms` — all 13 theorems

Every one of
`boxSlice_all_partners`, `boxSlice_all_partners_via_exBoxSum`, `mConvex_boxSlice`,
`mConvex_boxSlice_via_inS`, `inBox_iff_inS_saturated`, `kfun_le_of_lower`,
`sum_eq_of_inBox`, `mConvex_beadSupport`, `beadSupport_all_partners`,
`beadSupport_all_partners_witness`, `bead_hl_load_bearing`, `bead_sum_load_bearing`,
`bead_hl_not_necessary_at_l_zero`:

```
depends on axioms: [propext, Classical.choice, Quot.sound]
```

The standard three. No `native_decide`, no local axiom, no `sorry`.

## Sorries

**None.** The scan, and why I believe it:

- `command grep -rn --include=*.lean -E '\bsorry\b|\badmit\b'` over the **whole** source
  tree, 52 `.lean` files, `.lake` excluded. `command grep` and not `grep` because
  `lean/.gitignore` line 6 ignores `tworow_d4_kernel/` wholesale and the wrapped `grep`
  returns false absences (`memory: grep-here-honours-gitignore-and-my-lean-project-is-ignored`).
- **Validated against a planted canary before being believed.** I wrote
  `theorem canary_sorry : True := by sorry` into the tree, ran the scan, confirmed it was
  found, then deleted it. A clean scan from an unvalidated instrument is not a result.
- Result: `BoxSlice.lean` zero hits; all remaining tree-wide hits are docstring prose
  ("sorry-free", "No `sorry`"), none a tactic.

## Findings

### 1. The all-partners content already existed. (The main finding.)

`SublevelMConvex.ex_box_sum` (`SublevelMConvex.lean:170`, proved 2026-10-04 as the easy
half of `thm:M`) **is** the all-partners statement. Its own docstring says so: "For
**every** `l` with `y l < y' l` the surgered vector stays in the box and on the
hyperplane … This is already the statement that a box meets a hyperplane M-convexly."

My brief for this session checked for overlap against `mConvex_inS`, found none, and
concluded the strengthening "is genuinely not `mConvex_inS`" — which is true and beside
the point. `ex_box_sum` sits in the same file as the `ex` / `sum_ex` / `ex_apply_*`
vocabulary the brief instructed me to reuse, four lemmas above them. **A guard that
names one declaration is discharged by reading that declaration** — this is
`a-guard-that-names-a-field-is-dodged-by-the-next-field` at the granularity of Lean
names rather than JSON keys, and the brief's own pre-flight walk of the project missed it
because it was looking for the word "box slice", not for the statement.

What is left, and is real but small:

1. **Hypothesis weakening.** `ex_box_sum` assumes `InS P Q σ j` — box *and* `kfun y ≤ j` —
   and discards the bound (`⟨hyP, hyQ, hyS, -⟩`). So the box-slice statement was only
   available wrapped in the sublevel machinery (`kfun`, `negPart`) it does not use.
   `InBox` mentions neither. This is what lets Step A of (Q) cite the Lean side without
   dragging in `thm:M`.
2. **A named predicate**, so `MConvex {y | InBox P Q σ y}` is statable at all.
   `mConvex_boxSlice` existed in no form.
3. **The intended instance** `mConvex_beadSupport`, the point where the Lean side meets
   the `.tex` side.

So the honest description of this session is *restate with properly weakened hypotheses,
then instantiate*, not *prove the strengthening*. One slot, correctly sized, but for a
different reason than the brief gave.

### 2. `lem:constsum` was already formalised *and already pointed at*.

The brief anticipated half of this — it told me not to re-prove `lem:constsum` because it
is `ParityObstruction.sum_eq_of_mConvex`. Correct. It also told me to "give
`sum_eq_of_mConvex` its pointer wherever `lem:constsum` is recorded." **That pointer
already exists**, in `A-parity-obstruction-lean.lean`, written 2026-10-04. Task 2's second
half was done before the session began. Printing the stored field cost one line and saved
a duplicate pointer.

### 3. The node's "each spent once" accounting is correct in the direction that matters, and mildly over-stated.

The parent node says each of the four needed inequalities costs "one strict inequality
between `α, β` plus one box constraint, each spent once." Proving it directly rather than
via `ex_box_sum` is what let me check this. The actual ledger:

| goal              | strict ineq. | box constraint |
|-------------------|--------------|----------------|
| `P i ≤ y i − 1`   | `hi`         | `P i ≤ y' i`   |
| `y i − 1 ≤ Q i`   | —            | `y i ≤ Q i`    |
| `P l ≤ y l + 1`   | —            | `P l ≤ y l`    |
| `y l + 1 ≤ Q l`   | `hl`         | `y' l ≤ Q l`   |

Two of the four need a strict inequality; the other two need only `y`'s own box
constraint. Nothing beyond the claimed budget is spent and no hypothesis instance is used
twice, so the claim is safe; two goals are simply cheaper than advertised. No goal needed
more than its row plus `omega`, which is what the brief said to treat as the pass
condition. `hi` and `hl` are each spent once more, jointly, to derive `i ≠ l`.

### 4. I wrote a vacuous control and it closed by `rfl`. I removed it.

The brief asked me to do both routes to `mConvex_boxSlice` and "check they agree". I
stated `boxSlice_all_partners = boxSlice_all_partners_via_exBoxSum` and it closed by
`rfl` on the first build — because `InBox` is a `Prop` and Lean has **definitional proof
irrelevance**, so *any* two proofs of *any* proposition are equal by `rfl`. The theorem
cannot fail; it checks nothing. Banking it would have been a vacuous control filed as
corroboration (`memory: a-pass-count-does-not-report-the-rank-of-the-test`,
`a-silent-control-needs-its-silence-explained`).

Removed, with a note in the file saying why. What the two routes genuinely give is two
independent *derivations* of one statement; the sharing of the statement is visible in the
two signatures and is not a theorem. **"Check that two Lean proofs agree" is never a
check** — this generalises past this session and belongs in the instrument list.

The replacement real check is **ablation**, which does fire (below).

### 5. The first ablation run was a silent no-op that printed exit 0.

My `sed` for the ablation had the wrong indentation, matched nothing, and `lake build`
replayed a cached module: **exit 0**. Read naively that says *the `ex_box_sum` citation is
not load-bearing* — the exact direction the brief warned the `ConstantInfo.value?` walk
fails in, reproduced by a different mechanism. An ablation that does not ablate is
indistinguishable from a dependency that isn't there.

Fixed by asserting the pattern is present in the file *before* mutating, and grepping the
mutated line back out before building. Both ablations then fire:

| citation | ablated to | build |
|---|---|---|
| `SublevelMConvex.ex_box_sum` | `ex_box_sum_ABLATED` | **exit 1**, `Unknown identifier` |
| `ParityObstruction.mConvex_inS` | `mConvex_inS_ABLATED` | **exit 1**, `Unknown identifier` |

Restored, full rebuild exit 0.

Also: `lake` is not on `PATH` in this container (it lives at
`/home/clio/projects/.elan/bin/lake`). The first `lake build` returned
**`command not found`, exit 127** — caught only because I read `PIPESTATUS[0]`, exactly
the 10-03 failure the brief told me to guard.

## Deviation from the brief, and why

The brief said to set `Q-stepA-bead-support-box-slice` itself to `trust: lean-verified`.
**I did not.** That node's statement continues past the box slice: "…whence `log c = 0` on
an M-convex domain is M-concave and `normalizedcoefficients` makes `N(P̃)` Lorentzian."
Neither of those two clauses is formalised. Flipping the parent would have asserted Lean
coverage of the Lorentzian conclusion.

Instead I added a child, `Q-stepA-boxslice-lean`, trust `lean-verified`, scoped in its
`approach` to exactly what builds, and left the parent at `proved`. This matches the
repo's own convention — `A-parity-obstruction-lean` is a `lean-verified` child of a
`proved` parent for the same reason.

## Not attempted (unchanged from the brief)

`sigPos` ↔ count of positive eigenvalues (Mathlib has `sigPos`, no Cauchy interlacing);
anything from today's PROVE slot; `samuel-chain-formula` (no Schubert polynomials in
Mathlib).

### 6. `registry_validate --proofs-dir` takes `.`, not `proofs`, and the wrong value fails *everything*.

The brief's startup section says to run `python3 code/registry_validate.py <registry>`
and notes that "`trustcheck` takes `--files-dir`; `registry_validate` takes
`--proofs-dir`. Both are correct for their own tool — do not normalise them." True about
the *flag names*. But the *value* is the same for both: `.`, because node `file` fields
already begin with `proofs/`. Passing `--proofs-dir proofs` double-prefixes and reports
**111 nodes** as `file not found` — including my new one, which does exist.

So the first run looked like I had broken the registry. The discriminator is the count:
a real violation from a one-node edit cannot be 111. With `--proofs-dir .` both tools are
green, and both target files exist on disk (`2026-10-05-c3-lean-Q-stepA-boxslice.md`,
`2026-10-04-c3-Q-lorentzian.tex`). The flag-name warning in my brief is the right warning
pointed one level too shallow: **the hazard is the argument, not the flag.**

### Control on the validator itself

`one-planted-control-per-claim-the-tool-makes`: before believing `trustcheck`'s "OK", I
planted two mutations in the node I had just added.

| control | mutation | result |
|---|---|---|
| A | `file` → `proofs/THIS-FILE-DOES-NOT-EXIST.md` | **exit 1**, names the node and the path |
| B | `trust` → `totally-bogus-trust-level` | **exit 1**, invalid-trust *and* the boundary rule "premise children must be at least 'proved'" |

Both fire, so the green is licensed for file existence and for the trust vocabulary.
Control B also confirmed that `lean-verified` satisfies the parent's premise-boundary rule,
which is why attaching a `lean-verified` child under a `proved` parent validates.
