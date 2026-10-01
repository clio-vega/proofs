# LEAN snapshot — 2026-10-01 — Lemma T's two sign lemmas

**Project:** `projects/lean/tworow_d4_kernel/`, new module
`TworowD4Kernel/LemmaT.lean` (imported from the root `TworowD4Kernel.lean`, so it is inside
the axiom-audit closure).

**Paper proof:** `projects/proofs/2026-09-30-c2-lemma-T.tex` — `lem:signs` (l.166),
`lem:offJ` (l.238), `lem:decay` (l.157).

**Build:** `lake build` exit 0, **3177 jobs**. **Sorry count: 0.** (`grep -rn sorry` over the
project's non-`.lake` `.lean` files returns 7 lines; all 7 are the word "sorry" inside
docstring prose — *"single `sorry`-free declaration"*, *"No `sorry` in this file"* — and zero
are tactics. `sorryAx` does not appear in any `#print axioms` output.)

## Targets, both closed

| declaration | paper | axioms |
|---|---|---|
| `LemmaT.Setup.lem_signs` | `lem:signs`, l.166 | propext, Classical.choice, Quot.sound |
| `LemmaT.Setup.lem_offJ_lower` | `lem:offJ`, l.238 (the `G(s-1)` cell) | propext, Classical.choice, Quot.sound |
| `LemmaT.Setup.lem_offJ_upper` | `lem:offJ`, l.251 (the `G(s+1)` cell) | propext, Classical.choice, Quot.sound |

24 declarations in all; full `#print axioms` sweep in the file's `Audit` section. Every one is
on the standard three, except `v_step_left`, `v_step_right` (propext only) and `v_concave`,
`abs_r_le`, `abs_l_le` (propext, Quot.sound).

## The modelling decision that made this worth doing

`lem:signs`'s proof turns on *"the cell `(i₀-1, x₁+1)` lies in the box, so `v_{i₀-1}` is
**defined**"*. With `Φ, Ψ : ℤ → ℤ` total, definedness is free and **that step is a no-op** — the
lemma would typecheck with `(INT)` never used, and the formalisation would certify nothing about
the step that was actually in question.

So every hypothesis in `structure Setup` is **restricted to the box**: `Φ_concave`, `Φ_lip` hold
only on `[A,B]`, `Ψ_convex`, `Ψ_lip` only on `[C,D]`, and the cell-support description `supp`
only at cells of the box. `(INT)` is then load-bearing **twice**, visibly:

1. it puts `i₀-1` and `i₁+1` inside the box, so `supp` applies there and yields
   `v (i₀-1) ≤ 0`, `v (i₁+1) ≤ 0` (`Setup.v_lo`, `Setup.v_hi`);
2. it puts `i₀` and `i₁+1` inside `[A+1, B]`, the range where `Φ_lip` is available.

## Findings

**1. `lem:signs` uses neither concavity hypothesis.** It needs only the two Lipschitz bounds and
the support description — `Φ_concave` and `Ψ_convex` are never invoked in `lem_signs` or in the
`key_r`/`key_l` it rests on. The paper's remark says Lemma `lem:signs` "is where `Φ`'s Lipschitz
bound is spent, and it is spent *only* here", which is true and is the stronger half of the
observation; it does not say that *nothing else* is spent there. Worth a sentence in the paper.

**2. `lem:offJ` cannot cite `lem:signs`; it must cite the intermediate inequalities.** The
residual branch (`i' = i₀-1`) needs `r ≥ 1 - φ'(i₀)`, not `r ≥ 0`. In the paper these live inside
`lem:signs`'s proof with no label, and the prose at l.253 reaches back into that proof to reuse
them. They are now separate declarations `key_r`, `key_l`, which is what makes the citation
resolvable at all. This is the same class as the defect the lemma was written to repair: a
pointer from prose to a proposition, where the proposition is not where the pointer says.

**3. The `lem:offJ` bound is tight.** `Control.wit_offJ_tight` exhibits an interior instance in
which both off-`J` cells have value *exactly* `0`, so `≤ 0` cannot be strengthened to `≤ -1`.

**4. The repaired `lem:offJ` argument survives formalisation in both branches.** The original
*"by the decay lemma"* had no object off `J`; the repair routes through the index `i'` one step
inside `J` (where `v` is defined) and closes the residual against `|r| ≤ 1`. Both branches of
both halves are formalised. The asymmetry in the hypotheses is real and is now explicit: the
lower half needs only `A ≤ s-1-D`, the upper only `s+1-C ≤ B`.

**5. `lem:decay` needs no `(INT)`** and is stated for a general sequence concave on an interval
(`decay_left` / `decay_right`), since that is all its proof uses. `decay_right` is `decay_left`
reflected through `c ∘ neg`. `decay_left` is the only induction in the file.

## Negative controls

A hypothesis class nobody instantiates proves nothing. Three live controls, all in
`namespace Control`:

- **`wit` / `wit_interior` — non-vacuity.** `A=0, B=5, C=0, D=2, s=4, i₀=i₁=3, Φ ≡ 0,
  Ψ = |x-1|`, so `J = [2,4]` and `v = (0,1,0)`. Chosen so `(INT)` holds **and both** off-`J`
  cells (`i = 1` and `i = 5`) lie in the box: all four conclusions are exercised, none vacuous.
  `wit_slopes` gives `l = -1`, `r = 1`, so `lem_signs` is not asserting `0 ≤ 0`.
- **`noLo` — the `(INT)` clause `A ≤ i₀-1` is load-bearing.** It violates *exactly that one*
  clause (`noLo_only_lo_fails` records that the other three hold, and `x₁+1 ≤ D`, so `r` is
  still a genuine in-domain increment of `Ψ`), and `noLo_r_neg : noLo.r = -1`. So `lem:signs`
  fails. `noLo_l_nonpos` records that the *other* half still holds — the control says which
  clause feeds which inequality.
- **`Φ_lipschitz_needed` — `Φ`'s Lipschitz bound is load-bearing.** Drops only `Φ_lip`, keeps
  `(INT)` in full, and again gets `r = -1` (`Φ = (-10,0,0,…)` concave on `[0,6]`, `Ψ = -x`,
  `s = 4`, `I = [1,4]`). Written as an inline existential, because `Φ_lip` is a **field** of
  `Setup`, so a counterexample to it cannot be a `Setup`. First attempt at this control *was*
  written as a `Setup` instance with the `Φ_lip` field "discharged" by `absurd rfl rfl` — an
  incoherent control that would have compiled only by being unprovable; caught before building.

## Not formalised

`thm:main` (interior concavity) is **not** promoted. Still open in Lean:

- `lem:notrunc` (no truncation inside `I`) — needs `|δᵢ| ≤ 1`, `|εᵢ| ≤ 1` from `Ψ_lip` at each
  `i ∈ I`, plus the four box inclusions. Expected easy; it is bookkeeping of the same kind.
- `lem:bdry` (`E₋ ≤ r`, `E₊ ≤ -l`) — needs the four-cell case analysis **and** a finite-sum
  formalisation of `E±` over `J \ I`, which the current file has no `Finset` machinery for.
  This is the real cost, and it is indexing, not inequalities.
- `lem:tele` (telescoping `Σ ε - Σ δ = r - l`) — needs `Finset.sum_range_succ_sub` style
  telescoping over `[x₀, x₁]`.
- `thm:main` itself — one `omega` once the four are in place.

A green `thm:main` resting on `sorry`-ed lemmas is worth nothing, so the parent node
`lemma-T-interior-concavity` keeps `trust: proved`, with the four verified lemmas as
`lean-verified` children. Same standard as yesterday's restraint on `A-m2-four-regions`.

## Editorial repair, carried out

`2026-09-30-c1-cylindric-kostka-logconcavity.tex`, `lem:trunc`. Its proof said *"at the two ends
of the support the argument of Lemma `lem:conc` applies verbatim"*; `lem:conc` is stated for a
**bounded** support interval `J = [j₀,j₁]`, and `max(c,0)` need not have bounded support
(`c ≡ 1`; `posPart_support_unbounded` is the Lean witness). The citation's object is undefined on
that case.

The proof is now the quasiconcavity argument that `PFtwo_posPart` has already machine-checked,
transcribed into prose, plus `rem:truncLean` recording the declaration names and axioms. I also
**dropped** the sentence permitting `c` to take the value `-∞` off an interval: it is a separate
statement, it is not formalised, and it is not needed — at the only use site (`w_I`, `w_III` in
`prop:regI`) `c` is a genuine `ℤ`-valued concave sequence, being minus a maximum of finitely many
affine functions. Both papers recompile (`pdflatex` exit 0, 9 pages each).

I had deferred this on 09-30 under "no new mathematics in a LEAN session". That reasoning was
wrong: transcribing a proof Lean has already checked into the prose that claims it is not new
mathematics.

## `by Lemma X` sweep across both papers

Six sites. Five sound, one was the defect above.

| site | citation | object defined on the case? |
|---|---|---|
| c2 l.252 | `by Lemma lem:decay` at `i' ≤ i₀-2` | **yes** — `i' ∈ J` is established two lines earlier, and this is now `decay_left` |
| c2 l.290 | `by Lemma lem:offJ` | yes — `thm:main` assumes interior, which is `lem:offJ`'s hypothesis |
| c1 l.339 | `w_I ∈ PF₂ by Lemma lem:trunc` | yes — `w_I` is the positive part of a `ℤ`-valued concave sequence |
| c1 l.413 | `by Propositions prop:regI and prop:regII` | yes — `regI` covers regions I and III, `regII` covers II and IV; all four summands |
| c2 l.393 | `Lemma lem:bdry` | not a citation-step (a question in the discussion) |
| c1 l.125 | `the argument of Lemma lem:conc applies verbatim` | **NO** — repaired above |

## A false alarm I nearly recorded, and the instrument that caused it

`registry_validate.py` **never reads the `lean` field.** Grep it: the only occurrence outside
the argument parser is inside report formatting. Measured: promoting
`lemma-T-interior-concavity` to `lean-verified` with
`lean: TworowD4Kernel.LemmaT.Setup.thm_main_does_not_exist` prints
`OK ... is valid`, **exit 0**. The validator *is* alive — an invalid `trust` value gives
`1 problem(s)`, exit 1, node named; a missing `children` key likewise, which is how the four new
nodes were caught — but the `lean` pointer is entirely ungraded. A validator's name is a claim
and nothing grades the name.

So I wrote `code/registry_lean_resolve.py`: for every node's `lean` field, does the
fully-qualified name occur as a declaration in the Lean sources? It refuses the planted
violation (exit 1, both nodes named: the dangling pointer *and* a `lean-verified` node with no
`lean` field) on a file that `registry_validate.py` passes.

**And its first run reported 35 dangling pointers out of 160 — all 35 phantom.** The bug was in
my own `end` handling: a *named* `end` closing a **section** popped a namespace frame, so every
declaration after `end Audit` in a file lost its `TworowD4Kernel.` prefix. What makes this worth
writing down is the confirmation step. I hand-spot-checked five of the names with
`grep -rln <name> --include=*.lean lean/` and **all five returned zero hits**, which read as
independent confirmation. They were not independent and they were not confirmation: that
recursive grep returns nothing from `projects/` for *any* of these names — `addRibbon` greps 0
hits from there while sitting at `AbacusRibbon.lean:74` — while the same grep rooted one
directory deeper works. I have not diagnosed why; I only know the instrument is unreliable at
that root and the only check I trust here is the `os.walk`-based index.

The shape: **a cross-check that is broken toward the same answer as the thing it checks is not a
cross-check.** Both instruments were blind in a way that produced *absence*, and absence read as
a finding — [[an-absence-claim-is-ungraded]] and [[an-inventory-fails-toward-alarm]] arriving
together, with the alarming number doing the work of making me believe the spot-check. What
saved it was noticing `addRibbon` in a file listing I had printed for an unrelated reason.

After the fix: **160 `lean` pointers across 12 registry files, every one resolves, exit 0.** That
is the baseline, and it is a measurement, not an assumption.

## Registry

`proofs/registry/cylindric-lorentzian.json`:

- `lemma-T-offJ-cells` → **`lean-verified`**, `lean: lem_offJ_lower, lem_offJ_upper`
- new `lean-verified` children of `lemma-T-interior-concavity`: `lemma-T-signs-lean`,
  `lemma-T-key-slopes-lean`, `lemma-T-decay-lean`, `lemma-T-signs-controls-lean`
- `lemma-T-interior-concavity` itself **not** promoted; its note now says which four lemmas are
  missing.

`python3 code/registry_validate.py --proofs-dir . proofs/registry/cylindric-lorentzian.json`
→ `OK`, exit 0. `python3 code/registry_lean_resolve.py proofs/registry/cylindric-lorentzian.json`
→ 24 pointers, all resolve, exit 0.

## A second thing nobody can see

`code/` is **not a git repository** — `git rev-parse` from inside it fails up to the mount point.
So `registry_validate.py`, and every other registry tool, exists only in this container: the
artifacts whose whole job is to be trusted without re-checking are precisely the ones Robin
cannot read. `registry_lean_resolve.py` is therefore committed into the `proofs` repo at
`registry/registry_lean_resolve.py`, next to the registries it grades, with `code/` keeping the
working copy. Those two will drift; the pushed one is the canonical text.
