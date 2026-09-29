# LEAN 2026-09-29 c1 — B-EXC: Murota's symmetric exchange axiom

**Outcome: closed, sorry-free.** The last claim in the cylindric M-convexity registry that
was carried by brute force alone is now machine-checked.

## Target and project

Project: `clio-vega/tworow-d4-kernel`, at `/home/clio/projects/lean/tworow_d4_kernel`.
Baseline HEAD `870280c` (verified against `git ls-remote origin HEAD` →
`870280c7654c43638f995d7c62c7fcf609354c8e`, not quoted from recall).

New file: `TworowD4Kernel/SymmetricExchange.lean` (registered in `TworowD4Kernel.lean`).

```lean
theorem SymmetricExchange.insupp_symm_exchange {ℓ : ℕ} {lhat α β : ℕ → ℤ}
    (hl : IsPart ℓ lhat) (hα : InSupp ℓ lhat α) (hβ : InSupp ℓ lhat β)
    {i : ℕ} (hi : β i < α i) :
    ∃ j, α j < β j ∧ InSupp ℓ lhat (ex i j α) ∧ InSupp ℓ lhat (ex j i β)

theorem SymmetricExchange.sorted_symm_exchange {ℓ : ℕ} {lhat α β : ℕ → ℤ}
    (hl : IsPart ℓ lhat) (hα : SortedBridge.InSuppSorted ℓ lhat α)
    (hβ : SortedBridge.InSuppSorted ℓ lhat β) {i : ℕ} (hi : β i < α i) :
    ∃ j, α j < β j ∧ SortedBridge.InSuppSorted ℓ lhat (ex i j α) ∧
      SortedBridge.InSuppSorted ℓ lhat (ex j i β)
```

`ex i j α` is `α - eᵢ + e_j`, so `ex j i β` is `β + eᵢ - e_j`. The whole content is that a
**single** `j` serves both. The second theorem is stated for the paper's own sorted-form
`J = {α : sort(α) ⊴ λ̂}`, not only the subset-form surrogate.

## Where the proof comes from — say it plainly

**The paper does not contain this proof.** `2026-09-20-c1-cylindric-M-convexity.tex` proves
the one-sided axiom and stops; its `prop:perm-mconvex` does not deliver the symmetric
conjunct. What is formalised here is the classical argument for bases of an integral
submodular system, in the shape armed in `state/LEAN.md` §4. So this session is *not* a
transcription — it is new formalisation of a known polymatroid fact. Recording that
distinction is part of the result.

The Murota–Shioura theorem (the two axioms cut out the same class) is **not formalised and
not used**. It is not needed: B-EXC is proved outright.

## The argument

With `ρ(S) = Λ_{|S|}` and "tight" meaning `α(S) = ρ(S)`:

* `A = tightAvoid ℓ lhat α i` — the **union** of all `α`-tight sets avoiding `i`. Tight by
  `MConvexExchange.tight_union` (the paper's `lem:tight`), via `Finset.sup_induction` over
  the filtered powerset. And `α - eᵢ + e_j ∈ J` ⟺ `j ∉ A`.
* `B = minTight ℓ lhat β i` — the **intersection** of all `β`-tight sets containing `i`.
  Tight by the new `tight_inter`. And `β + eᵢ - e_j ∈ J` ⟺ `j ∈ B`.
* `asum_sdiff_le` — the squeeze:

      α(B \ A) = α(A ∪ B) − α(A) ≤ ρ(|A∪B|) − ρ(|A|) ≤ ρ(|B|) − ρ(|A∩B|) ≤ β(B \ A)

  the middle step being `MConvexExchange.rank_submodular`, where submodularity of `ρ` is
  **derived** from `λ̂` antitone, not postulated.
* `exists_gain_in_sdiff` — `i ∈ B \ A` contributes `β i − α i < 0` to a sum that is `≥ 0`,
  so some other `j ∈ B \ A` has `α j < β j`. That `j` is the witness.

Both memberships then come from the *existing* `MConvexExchange.insupp_ex_of_no_tight`,
applied twice with `(i,j)` swapped — no new surgery lemma was needed, because
`β + eᵢ - e_j` is literally `ex j i β`.

### One Lean-specific obstacle worth recording

`Finset ℕ` has **no `⊤`**, so `Finset.inf` is unavailable and the intersection of a family
cannot be written directly. `minTight` is therefore defined as the complement, inside
`[ℓ]`, of the **union of the complements** — and the closure property that union needs is
exactly `tight_inter`. `Finset.sdiff_union_distrib` and `Finset.sdiff_sdiff_eq_self` carry
the bookkeeping.

## Sorry / axiom status

**Zero sorries**, not even in a docstring. `grep -n sorry TworowD4Kernel/SymmetricExchange.lean` → no match.

`lake build`: **Build completed successfully (2993 jobs)**, exit 0.

`#print axioms` (pasted, not asserted):

```
'SymmetricExchange.insupp_symm_exchange' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.sorted_symm_exchange' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.symm_conjunct_not_automatic' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.symm_exchange_fires_on_control' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.sorted_symm_exchange_witness' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.tight_inter' depends on axioms: [propext, Classical.choice, Quot.sound]
'SymmetricExchange.asum_sdiff_le' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Exactly the standard three. Nothing else was wanted.

## The negative control — the part that makes the theorem have content

`state/LEAN.md` §5 was right to demand this first. B-EXC is stated as Murota states it (a
literal second membership), but that alone does not show the strengthening is real: if the
symmetric conjunct followed from the one-sided one, `insupp_symm_exchange` would be a
restatement of `insupp_exchange` and `tight_inter` / `minTight` / the squeeze would all be
dead weight.

It does not follow. `SymmetricExchange.symm_conjunct_not_automatic`, machine-checked:

| | |
|---|---|
| `ℓ = 4`, `λ̂ = (2,1,1,0)` | `α = (1,0,1,2)`, `β = (0,1,2,1)`, `i = 3` (`β₃ = 1 < 2 = α₃`) |
| `j = 1` | `α − e₃ + e₁ = (1,1,1,1) ∈ J` ✔ one-sided … `β + e₃ − e₁ = (0,0,2,2) ∉ J` ✘ (two largest sum to `4 > Λ₂ = 3`) |
| `j = 2` | `α − e₃ + e₂ = (1,0,2,1) ∈ J` ✔ … `β + e₃ − e₂ = (0,1,1,2) ∈ J` ✔ |

So `insupp_exchange` may legitimately return `j = 1` and `insupp_symm_exchange` may not.
The theorem constrains the witness and the constraint bites.

This is the **minimal** configuration: exhaustive search finds no negative control at all
for `ℓ ≤ 3`; at `ℓ = 4` the smallest `|λ̂|` admitting one is `4`, and `λ̂ = (2,1,1)` is the
unique partition of `4` that does. (Contrast 2026-09-26, where the brief's negative control
turned out to be *impossible* because the dropped half-space was redundant. Here it exists.)

## The §3 pre-check, and why it was the right order

`proofs/code-q254-lean/bexc_symmetric_check.py` — written from the definitions with `J`
built in the **sorted** form, independent of the Lean file. Range `|λ̂| ≤ 9, ℓ ≤ 4`
(`check-bexc-d9-ell4.out`; a wider `|λ̂| ≤ 11, ℓ ≤ 5` run was still going at session end):

```
(EX2) B-EXC, symmetric exchange : 0 failures / 876317 triples (alpha,beta,i)
(AUT) one-sided j that also works for beta : 1372412 / 1378676 quadruples
      => symmetric conjunct AUTOMATIC? NO  (6264 one-sided j fail for beta)
(NB)  predicted witness set D cap (B\A) == actual : 876317 / 876317 triples, 0 mismatches
```

Three things, and the third is the one I would not have thought to ask for:

1. **(EX2)** B-EXC holds — 876317 triples, the same count as the historical baseline.
2. **(AUT)** answered the question §3 actually posed: does the existing one-sided
   construction's `j` already work both ways? **No** — 6264 one-sided witnesses fail
   symmetrically in this range alone. So the construction *had* to change, and I knew that
   before writing a tactic rather than after.
3. **(NB)** — the check I added: the witness set my proof *predicts*,
   `{j : α_j < β_j} ∩ (B \ A)`, equals the true set of good `j` **exactly**, not merely
   nonemptily. A construction validated only by nonemptiness would pass with the wrong
   witness set. 0 mismatches over all 876317 triples said the extremal-sets picture was
   right before I tried to make Lean agree, and the formalisation then went through in one
   pass with three small errors, all of them Mathlib argument-name slips.

## Caveat sweep — `grep` the stem, not the phrase

Four locations were listed in `state/LEAN.md` §1 as naming this gap. Grepping the stems
(`B-EXC`, `Murota`, `876317`) rather than any English phrase found a **fifth** the brief
did not list: the registry node `snp-of-the-cylindric-skew-schur-polynomial`, whose
`approach` ends *"B-EXC (Murota's symmetric axiom) is STILL OWED and is not touched by this
closure."* True on 09-26, false now, and sitting at a node whose trust is `lean-verified` —
i.e. exactly where the next reader would copy from. Fixed in the same commit as the proof.

Updated (all in this commit):

* `TworowD4Kernel/MConvexExchange.lean` docstring — the "not proved here … the paper's
  argument does not deliver it either" paragraph now points at `SymmetricExchange`.
* `TworowD4Kernel/SortedSubsetBridge.lean` docstring, §"What Lean still does not see".
* `proofs/registry/cylindric-M-convexity.json` — `polymatroid-exchange`,
  `polymatroid-exchange-subset-form`, **and** `snp-of-the-cylindric-skew-schur-polynomial`;
  plus a new child `polymatroid-exchange-symmetric-bexc`, trust `lean-verified`,
  `lean: SymmetricExchange.sorted_symm_exchange`.

The two commit messages (`cc0573b`, `85054bb`) that also name the gap are immutable; the
new commit message says so and points forward, which is the only repair available.

## Validators — watched to refuse before recording that they hold

`registry_validate.py` refused first: *"missing key 'children'"* on my new node. Then
`trustcheck.py` refused: *"invalid role … must be 'premise' or 'attempt'"*. Both fixed;
both now clean:

```
$ python3 code/trustcheck.py --deployment code/clio.json --sources skip --chunks-dir skip \
    validate proofs/registry/cylindric-M-convexity.json --files-dir .
OK: proofs/registry/cylindric-M-convexity.json is valid (status: proved, deployment: clio)   [exit 0]

$ python3 code/registry_validate.py proofs/registry/cylindric-M-convexity.json --proofs-dir .
16 warning(s); OK: proofs/registry/cylindric-M-convexity.json is valid (status: proved)
```

Two standing tool defects re-confirmed, neither new:

* **`registry_validate.py` exits 0 with problems printed.** It printed `1 problem(s)` and
  returned 0. A detector that cannot fail cannot gate anything — read the output, never the
  exit code.
* **the root-dir bug is still live in `registry_validate.py`.** Default invocation gives 20
  spurious `file ... not found under /home/clio/projects/proofs`, including for my new file.
  `--proofs-dir .` is correct (`trustcheck`'s spelling is `--files-dir .`). Confirmed again
  by running both forms, not by recall.

* The 16 remaining warnings are all the known non-firing extraction gate: sources written in
  the ID + locator form my own protocol prescribes fail to resolve and are reported as soft
  `not in the sources index`. Pre-existing, untouched today, and still the case that the
  gate never actually binds.

## What is now owed

Nothing in this registry is carried by brute force alone. The honest remaining items are
unchanged from 09-25 c2 and are *not* about B-EXC:

* `polymatroid-exchange`'s re-review by Rick is still owed on the `a6c83ed` build (his
  endorsement predates the erratum to the displayed sign).
* `AffineAdditive` still *defines* the paper's Lemma 2.3 criterion rather than deriving it
  (affine symmetric group not in Mathlib) — a different registry.
