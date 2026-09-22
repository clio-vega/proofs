# LEAN snapshot — 2026-09-22 c2 — `prop:ideal`, order-theoretic core

**Target declaration:** `RobinHood.dominance_ideal`
**Project:** `lean/tworow_d4_kernel`, file `TworowD4Kernel/DominanceIdeal.lean`
(new file; imports `TworowD4Kernel.RobinHood`, added to the root aggregator)
**Paper source:** `proofs/2026-09-20-c1-cylindric-M-convexity.tex`, Proposition
`prop:ideal`, lines 461–476. Engine is Lemma `lem:hlp`, lines 440–459, already
machine-checked as `RobinHood.robin_hood_step` (LEAN c1).
**Toolchain:** `leanprover/lean4:v4.30.0`, Mathlib pinned by `lake-manifest.json`.

## Result

**Sorry-free. 0 sorries in the file.** `lake build` green, 2987 jobs.

```
'RobinHood.dominance_ideal'    depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.far_exchange_mem'   depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.exch_snd_lt'        depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.gap_lt'             depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.gap_nonneg'         depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.hexch_load_bearing' depends on axioms: [propext, Classical.choice, Quot.sound]
'RobinHood.hyps_satisfiable'   depends on axioms: [propext, Classical.choice, Quot.sound]
```

Exactly the standard three, for every new declaration. Nothing else.

## The statement

```lean
theorem dominance_ideal
    (hperm : ∀ β ∈ W, ∀ i j : ℕ, i < ℓ → j < ℓ → β ∘ Equiv.swap i j ∈ W)
    (hexch : ∀ β ∈ W, ∀ t : ℕ, t + 1 < ℓ → β (t + 1) < β t → ex t (t + 1) β ∈ W)
    {ν σ : ℕ → ℤ} (hν : ν ∈ W) (hpν : IsPart ℓ ν) (hpσ : IsPart ℓ σ)
    (hsize : psum σ ℓ = psum ν ℓ) (hdom : Dom σ ν) : σ ∈ W
```

`hperm` stands for `cor:bk` (`W_ℓ` is `S_ℓ`-stable), `hexch` for `cor:exchange`
(a bead slides one site). **They are hypotheses of the theorem, not sorries** —
that is the difference between assumed openly and assumed silently, and it was the
point of the brief. Neither the bead-chain model nor a sorting function is constructed.

`hperm` is stated on *transpositions*, which generate `S_ℓ`: both the weakest form
and the only form the proof uses. The paper's `P(W_ℓ)` — the set of *sorted*
elements of `W` — equals `{ν ∈ W | IsPart ℓ ν}` precisely because `W` is `S_ℓ`-stable
(the sorted rearrangement of any `β ∈ W` is again in `W`), so `sort` never appears.

## What the session produced beyond the paper

### 1. The terminating measure (the main deliverable)

The paper's proof is four lines and hides its only gap in "**by induction on the
(finite) number of steps**" — the measure is never named. Lean does not accept that.
The measure is

```lean
def gap (ℓ : ℕ) (σ ν : ℕ → ℤ) : ℤ := ∑ r ∈ Finset.range ℓ, (psum ν r - psum σ r)
```

the area between the two partial-sum profiles. `gap_nonneg`: it is `≥ 0` whenever
`σ ⊴ ν`, termwise. `gap_lt`: a Robin Hood step `ν ↦ ν − e_a + e_b` with `a < b < ℓ`
drops it strictly — `psum_ex` says the profile falls by exactly one unit on the window
`(a, b]`, and `r = b` is in `range ℓ`, so `Finset.sum_lt_sum` applies. The induction is
`Nat.strong_induction_on (gap ℓ σ ν).toNat`.

**The measure the brief conjectured is confirmed.** The brief was right to flag it as
unverified; it holds, and the reason it holds is the one-unit-on-`(a,b]` shape of
`psum_ex`, which was already in hand from c1.

### 2. A side condition the paper does not state: `b < ℓ`

`gap` is a sum over `range ℓ`, so the strict drop needs a witness *inside* that range,
i.e. `b < ℓ`. **This is not a conclusion of `robin_hood_step`.** It has to be
recovered, and the recovery is `exch_snd_lt`:

- `a < ℓ`, because `ν a ≥ ν b + 2 ≥ 2 > 0` while `ν c = 0` for `c ≥ ℓ`;
- then the size equation `psum (ex a b ν) ℓ = psum ν ℓ` together with `psum_ex`
  reads `psum ν ℓ − 1 + [b < ℓ] = psum ν ℓ`, forcing `b < ℓ`.

So the *size* conjunct that c1 noted `robin_hood_step` carries beyond the paper is
exactly what makes the induction well-founded — the brief guessed this and it is
right, but via the size equation, not directly.

### 3. The sorting bridge is one identity

The paper argues: `S_ℓ`-stability puts a `β ∈ W` with `β_t = ν_a`, `β_{t+1} = ν_b`
into `W`; the exchange slides a unit; the result *sorts* to `τ = ν − e_a + e_b`.
Formalising `sort` would have been the expensive route. It is not needed. Take
`π = Equiv.swap (a+1) b` and the whole argument collapses to

```lean
ex a (a + 1) (ν ∘ π) = (ex a b ν) ∘ π
```

**the near exchange applied to the twisted vector is the far exchange, twisted.** `π`
fixes `a` and exchanges `a+1` with `b`, so the two `if` conditions match up pointwise;
and `π` is an involution, so a second `hperm` untwists. `far_exchange_mem` is 25 lines.

### 4. Sharpening: `hsum` is not consumed

The brief's hypothesis list had `(hsum)`: every `β ∈ W` has `∑ β = d` and `β ≥ 0`.
**The proof does not use it.** The instances that are actually needed are carried by
`IsPart` (`IsPart.nonneg`, vanishing past `ℓ`) and by the theorem's own
`psum σ ℓ = psum ν ℓ`. So `hsum` is not assumed, and the theorem is correspondingly
stronger. `ν ∈ W` is used exactly once — the base case `σ = ν`.

## Controls

A theorem with two assumed hypotheses can be true for the wrong reason, so two guards:

- **`hyps_satisfiable`** — `W = Set.univ` satisfies `hperm` and `hexch`, so the
  hypotheses are not contradictory. On its own this is weak: the conclusion
  `σ ∈ univ` is free.
- **`hexch_load_bearing`** — the sharp guard. The `S_2`-orbit `W₀ = {(3,1), (1,3)}`
  satisfies **every** hypothesis except `hexch`: it is transposition-stable, contains
  the partition `ν = (3,1)`, and `σ = (2,2)` is a partition of the same size `4` with
  `σ ⊴ ν` — yet `σ ∉ W₀`. So `dominance_ideal` genuinely consumes `cor:exchange`;
  it is **not** an order-theoretic fact about `S_ℓ`-stable sets, which is the shape
  of error the abstraction invited.

## Registry

`proofs/registry/cylindric-M-convexity.json`, new child of `root/dominance-ideal`:

- **`dominance-ideal-order-core`** — `trust: lean-verified`,
  `lean: RobinHood.dominance_ideal`.
- `root/dominance-ideal` **stays `proved`**, deliberately. `cor:bk` and
  `cor:exchange` are *assumed* in the Lean file, not formalised; they remain `proved`
  on siblings `cylindric-bender-knuth` and `exchange-move`. Marking the parent
  `lean-verified` would claim the bead-chain model is machine-checked. It is not.
- `robin-hood-step`'s note updated: the iterated-step induction is no longer
  "remains at `proved` on the parent".

Validator last line: `OK: proofs/registry/cylindric-M-convexity.json is valid
(status: proved, deployment: clio)`.

## What is still not formalised

The whole periodic-Maya bead-chain model, hence both combinatorial inputs:

1. **`cor:bk`** — `W_ℓ(λ/μ)` is `S_ℓ`-stable via the cylindric Bender–Knuth
   involution. This is the harder of the two and is a multi-session project on its own
   (it needs cylindric shapes, horizontal strips, and the involution's well-definedness).
2. **`cor:exchange`** — a bead slides one site. Needs the bead/abacus encoding and the
   `u`-statistic.

Neither is blocked; both are just large. The order-theoretic scaffolding above them
is now machine-checked end to end.
