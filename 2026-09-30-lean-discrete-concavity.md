# LEAN snapshot — 2026-09-30 c2 — the discrete-concavity toolkit under the m=2 theorem

**Project:** `lean/tworow_d4_kernel` (Lean 4 v4.30.0, Mathlib v4.30.0), repo
`clio-vega/tworow-d4-kernel`.
**Module:** `TworowD4Kernel/DiscreteConcavity.lean`, imported from the root
`TworowD4Kernel.lean` (line 43) so it is inside the axiom-audit closure.
**Paper proof:** `proofs/2026-09-30-c1-cylindric-kostka-logconcavity.tex`.
**Commit:** local `a9911280458001331523f8821fedb0e2550cce4b`,
`git ls-remote origin HEAD` → `a9911280458001331523f8821fedb0e2550cce4b`. Printed, not inherited.

## Build

`lake build` exit 0, **3176 jobs**. The module's own compile time: **39 s** (115 s on the first,
cold pass, which included elaborating `import Mathlib.Tactic`). The one warning is
`Files in mathlib cannot import the whole tactic folder`, shared with five other modules in this
library including `MetricCriterion`; it is a style linter, not a failure.

## Sorry count: 0

Not one, not even in a docstring. Nothing was left as a bookmark; nothing was removed to make
the count look better either — see *What is not formalised* below, which is also written into the
module docstring in the section a copier reads.

## Targets, all sorry-free

| declaration | paper | what it is |
|---|---|---|
| `PFtwo` | l.103 | the `PF₂` predicate as a 3-field structure: `nonneg`, `suppInterval`, `logConcave` |
| `PFtwo_of_concave_on_interval_support` | `lem:conc`, l.106 | concave + nonneg + interval support ⟹ `PF₂` |
| `PFtwo_posPart` | `lem:trunc`, l.119 | the positive part of a concave sequence is `PF₂` |
| `slope_antitone`, `IntConcave.min_le`, `IntConcave.pos_support_interval` | the ungraded clause of `lem:trunc` | "`{c>0}` is an interval **because** `c` is concave" |
| `posPart_min` | `prop:regII` | `min(α,β)₊ = min(α₊,β₊)`, the identity the trapezoid formula leans on |
| `rho_zero`, `rho_mono`, `rho_concave`, `rho_eq_two_caps` | `prop:regII` | the three links of the composition chain, stated separately |
| `H_concave`, `H_pos`, `H_nonpos`, `GR_eq_rho_comp_H` | `prop:regII` | the profile `H(s)=min(s−ℓ+1, h−s+1)` and the factorisation `G_R = ρ_R ∘ H` |
| `GR_concave_interior`, `GR_PFtwo` | `prop:regII`, l.347 | regions II and IV are `PF₂` |

`PFtwo` is stated as the paper writes it, **not** transported through a Mathlib log-concavity
definition, and the three conjuncts are separate fields so that the name `PFtwo` is never
load-bearing. `logConcave` is written `g s * g s` to stay inside `omega`'s language;
`PFtwo.logConcave_sq` is the paper's `g s ^ 2` form.

### `#print axioms`

Every audited declaration (22 of them, printed in the `Audit` section at the end of the module)
returns **`[propext, Classical.choice, Quot.sound]`** or a subset. Several return strictly less:
`posPart_min`, `rho_zero`, `rho_eq_two_caps`, `H_concave`, `GR_eq_rho_comp_H` and `bump_PFtwo`
return `[propext, Quot.sound]`. No fourth axiom anywhere.

## Negative controls — five, all live

1. `gapTwo_not_PFtwo` — `(1,0,0,1)` is nonnegative **and** log-concave with non-interval support.
   So the interval conjunct of `PF₂` is not implied by the other two: it is load-bearing, which
   matters because `(P3)` (closure under convolution) is what consumes it.
2. `gapOne_not_logConcave` — a gap of width **one**, `(1,0,1)`, already breaks log-concavity at
   the gap. Together with (1) this pins down exactly how wide a hole `PF₂` tolerates: none, and
   two.
3. `cap_not_concave_at_zero` — `v ↦ min(v,k)₊` is concave on `v ≥ 1` and genuinely **not** at
   `v = 0`. This is precisely why `prop:regII` needs `H(s) ≥ 1` on `[ℓ,h]` and not merely `≥ 0`;
   the hypothesis is not decoration.
4. `bump_PFtwo` + `bump_not_concave` + `bump_concavity_failure` — `(1,3,6,7,6,3,1)` is `PF₂` and
   is not concave (`2·3 = 6 < 7 = 1+6`). `lem:conc` is strictly one-way, which is why
   `prop:regII`'s route — *prove concavity on the support* — cannot be expected to reach past
   `m = 2`.
5. **`Gsharp_not_logConcave`, the `rem:sharp` witness (l.272, first bullet) — the most valuable
   one.** `Φ = (−1,0,1)` is concave and 1-Lipschitz (`Phi_one_lipschitz`); `Ψ = (−3,0,4)` is
   convex and **not** 1-Lipschitz (`Psi_not_one_lipschitz`: the step `1→2` is `4`); and
   `G(s) = Σ_{i+j=s} (Φ(i) − Ψ(j) + 1)₊` is not log-concave. `Gsharp_values` computes
   `G = (3,4,6,2,0)` **by the kernel** (a `norm_num` over the Finset sum, not a copied numeral)
   and `Gsharp_failure` prints the instance `4·4 = 16 < 18 = 3·6` at `s = 1`. So the 1-Lipschitz
   hypothesis of Lemma T (`conj:T`) cannot be dropped, and the failure is attributable to `Ψ`
   alone.

## Two findings

### (a) A citation whose object is undefined on the case it is cited for

`lem:trunc`'s proof says *"at the two ends of the support the argument of Lemma `lem:conc`
applies verbatim."* But `lem:conc` is stated for support a **bounded** interval `J = [j₀,j₁]`,
and `max(c,0)` need not have bounded support: `posPart_support_unbounded` exhibits a concave `c`
— namely `c ≡ 1` — for which no `(j₀,j₁)` describes the support of `max(c,0)`, so `lem:conc` has
no "two ends" to offer.

The **conclusion** of `lem:trunc` is true and `PFtwo_posPart` proves it — but *without* invoking
`lem:conc`, via `IntConcave.min_le` (quasiconcavity, by strong induction on the span, on top of
`slope_antitone`, the only induction in the file). So the defect is in the pointer, not in the
mathematics. This is the exact shape of
`a-hypothesis-can-be-load-bearing-in-exactly-one-place`'s twin: *"by Lemma X" on a case where
X's object is undefined*. A sweep grades the conclusion and never a citation inside the proof,
which is why it survived.

Repair for the .tex: either state `lem:conc` for an interval `J ⊆ ℤ` not necessarily bounded, or
delete the citation from `lem:trunc` and give the quasiconcavity argument, which is two lines.
Both are one-paragraph edits. **Not done in this session** — the session rule is no new
mathematics, and the choice between them is an editorial call about `prop:regI`'s consumer.

### (b) The composition lemma is not needed — the brief predicted this and was right

`prop:regII` closes with *"a nondecreasing concave function composed with a concave function is
concave, so `G_R = ρ_R ∘ H` is concave on `[ℓ,h]`."* That is a real lemma and a real Mathlib hunt.
It is **not needed**. `GR_concave_interior` is a termwise explicit certificate: on `ℓ < s < h` all
three of `H(s−1), H(s), H(s+1)` are `≥ 1`, so the positive part is inactive on **every** summand,
and what remains per summand is

`min(H(s−1), k) + min(H(s+1), k) ≤ 2·min(H(s), k)`,

linear arithmetic on `min`s of affine functions — one `omega` each, under `Finset.sum_le_sum`.
The whole proof is five lines. The factorisation the paper *names* is not lost: it is recorded as
`GR_eq_rho_comp_H` (by `rfl`), and `rho_zero`/`rho_mono`/`rho_concave` state the paper's three
links separately, as the paper does.

The brief flagged the tell in advance — *"This brief says **find** for exactly one thing (the
composition lemma), so treat that line with suspicion"* — and that is the second session running
where *find* marked a surplus step. Yesterday it was "multilinear hence minimised at a vertex".
The rule is holding: **when a step invokes a named general lemma and the instance is a small
concrete inequality about explicit `min`s or polynomials, try for a certificate before going to
look the lemma up.**

## What is NOT formalised — stated here and in the module docstring

* **The trapezoid formula itself.** `prop:regII` derives
  `(1_[a,b] * 1_[c,d])(s) = min(s−ℓ+1, h−s+1, n_t, m_t)₊` from a convolution of two interval
  indicators. Convolution on `ℤ` is not defined here; `GR` is **defined** by the right-hand side.
  Everything the paper deduces *from* that formula is formalised; the derivation *of* it is not.
  The one identity the derivation leans on in passing, `min(α,β)₊ = min(α₊,β₊)`, **is**
  (`posPart_min`).
* **The `∓∞`-extended `lem:trunc`.** The paper writes *"the same proof applies."* That is a
  proposition, not a remark. Only the finite version (`c : ℤ → ℤ`) is proved. The extended
  version is a separate statement and is **not** proved here. One lemma does not stand in for two.
* **Lemma T** (`conj:T`) is a conjecture in the paper and is not formalised. Only its sharpness
  witness is.
* **`prop:regI`** (regions I and III, which route through `(P3)`, closure of `PF₂` under
  convolution) is not formalised — it needs convolution. Its *concavity* input is
  `PFtwo_posPart`, which is formalised.
* Consequently **`thm:regions` is not Lean-verified**, and the node `A-m2-four-regions` is left
  at `proved`, not promoted. Two of its four regions are verified, not four.

## Registry

`proofs/registry/cylindric-lorentzian.json` (Q255). Five new children under
`A-m2-four-regions`, each `lean-verified` with its declaration name in the `lean` field:
`lem-conc-concave-implies-pf2`, `lem-trunc-positive-part`, `prop-regII-constant-support`,
`pf2-conjuncts-independent`, `lemma-T-lipschitz-sharpness-lean`. Backup at
`.bak-0930c2-lean`. The parent `A-m2-four-regions` was **not** promoted, for the reason above.

### The validator's default `--proofs-dir` is the broken one

The brief was right that the flag is `--proofs-dir`, not `--files-dir`. It did not say that the
**default value** is wrong, and it is — the default is "parent of the registry's directory", i.e.
`/home/clio/projects/proofs`, against which every `file:` field of the form `proofs/X.tex`
resolves to `proofs/proofs/X.tex`. Run from `/home/clio/projects`:

```
python3 code/registry_validate.py proofs/registry/cylindric-lorentzian.json      → 50 problem(s), exit 1
python3 code/registry_validate.py --proofs-dir . proofs/registry/...json         → OK, exit 0
```

**This broke the negative control the first time I ran it.** I planted a missing file, got exit 1
with 50 problems — and the *unmodified* registry also gives exit 1 with 50 problems. Same exit
code, same count; the planted line was distinguishable only by reading its name out of a 50-line
list. A control whose signal is identical to its baseline is not a control. Re-run against a
clean baseline it bites exactly:

```
--proofs-dir . , unmodified   → OK, exit 0
--proofs-dir . , one planted  → 1 problem(s), exit 1, naming the planted node
```

Same family as the `trustcheck --files-dir proofs` defect of 09-22 — *this codebase resolves
registry `file:` paths against the wrong root* — and I had already recorded that
`registry_validate.py` shares it. What is new: **the wrong root is this tool's default, so the
plain invocation is the broken one**, and the breakage is loud enough (50 lines) to mask a
planted violation rather than reveal it. Canonical invocation, from `/home/clio/projects`:
`python3 code/registry_validate.py --proofs-dir . proofs/registry/<name>.json`.
