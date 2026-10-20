# Lean — rows 3 and 4 of the type-dependence table: one substitution, two disjoint loci

**Session:** LEAN 2026-10-20 c2
**Project:** `~/projects/lean/tworow_d4_kernel` (own repo)
**Module:** `TworowD4Kernel/TypeDependenceLocus.lean`, namespace `TypeDependenceLocus`,
registered in the root aggregator `TworowD4Kernel.lean:76`
**Question:** Q534 (`memory/questions/2026-10-20-dream-ledger.md`, L50–L78) — *are the five
type-`B`/`C` differences one substitution `w² ↦ w`, localising at `v = −1`?*  This is the
pre-registered **cheapest arm: row 3 against row 4.**
**Status:** sorry-free, `grep -c sorry` = **0**, standard three axioms.

---

## 0. Correction to the brief, measured before building

`state/LEAN.md` says row 4 (`ψ_X(v) = (1 − q v^h)/(1 − q)`) is already formalised in
`HeadExclusionRootsOfUnity.lean`, and that row 3 (the head factor, `(v−1)` vs `(v²−1)`) is
**"Not yet in Lean."**  **That is backwards**, and the module's own docstring says so:

- `HeadExclusionRootsOfUnity.lean:91` — `def headFactor (q v : K) (h : ℕ) : K := (1 - q) * (v ^ h - 1)`,
  docstring L88–L90: *"`lem:typecompare` (L291): `(1-q)(v-1)` in type `B` (`h = 1`),
  `(1-q)(v²-1)` in type `C` (`h = 2`)."*  That **is row 3**.
- Row 4's `ψ_X` has numerator `1 − q v^h` (abstract L45; `cor:rhoB` L523 for the type-`B`
  case).  `1 − q v^h` is a different polynomial in `v` from `v^h − 1`, and it was **not in
  Lean** before this session.

Consequence for the plan: the brief's step 2 could not be stated against the existing
`headFactor` at all.  Row 4 had to be defined first (`psiNum`), and then the comparison turned
out to go the other way (§2).  **The whole substantive result of the session came from the
half of the brief that was wrong.**

Also stale: the brief says *"Next free question number is 536."*  `questions/NEXT-FREE.md`
reads **540** — BROWSE c2 allocated 537–539 after the brief was written.  No number was spent
this session; no `Q<n>` appears in any identifier.

---

## 1. What builds, sorry-free

`lake build TworowD4Kernel.TypeDependenceLocus` → **exit 0**, 1889 jobs.  No `sorry`.

### The unification half of Q534 — **holds**

```lean
def row3Shape (q u : K) : K := (1 - q) * (u - 1)      -- lem:typecompare, L352–L357
def row4Shape (q u : K) : K := 1 - q * u              -- abstract L45, cor:rhoB L523
def psiNum   (q v : K) (h : ℕ) : K := 1 - q * v ^ h   -- NEW: row 4 was not in Lean

theorem headFactor_eq_row3Shape (q v : K) (h : ℕ) : headFactor q v h = row3Shape q (v ^ h) := rfl
theorem psiNum_eq_row4Shape     (q v : K) (h : ℕ) : psiNum   q v h = row4Shape q (v ^ h) := rfl

theorem rows_share_substitution (q v : K) (h : ℕ) :
    headFactor q v h = row3Shape q (v ^ h) ∧ psiNum q v h = row4Shape q (v ^ h) := ⟨rfl, rfl⟩
```

Both are `rfl`, over a bare `CommRing`, with **no** hypothesis — char-free, domain-free,
`q`-free.  Each row is a type-**independent** shape precomposed with the monomial `v ↦ v^h`,
so `B` vs `C` is `h = 1` vs `h = 2` of **one** family in both rows.  That is the half of
Q534's claim that survives, and `#print axioms` returns *does not depend on any axioms*.

This is weaker than it looks and I want to be exact about how.  `rem:G5` of the source paper
(L653–L666) already says the uniform `ψ_X` formula is *"a coincidence of outputs, not one
derivation with a parameter — no step of either proof has `h` as a variable."*  Lean confirms
the **outputs** are one `h`-family.  It says nothing about the derivations, which is where
`rem:G5` locates the disunity.

### The localisation half — **holds for row 3 only**

```lean
theorem boundaryLocus_subset_of_dvd (K) [CommRing K] (hd : h₁ ∣ h₂) :
    boundaryLocus K h₁ ⊆ boundaryLocus K h₂                       -- no hypotheses
theorem boundaryLocus_two_sdiff_one [IsDomain K] (hchar : (-1 : K) ≠ 1) :
    boundaryLocus K 2 \ boundaryLocus K 1 = {-1}
theorem boundaryLocus_symmDiff_of_dvd (K) [CommRing K] (hd : h₁ ∣ h₂) :
    (boundaryLocus K h₂ \ boundaryLocus K h₁) ∪ (boundaryLocus K h₁ \ boundaryLocus K h₂)
      = boundaryLocus K h₂ \ boundaryLocus K h₁
```

So for row 3 the general `h₁ ∣ h₂` statement the brief asked for is there: the symmetric
difference collapses to the one-sided difference because divisibility nests the loci, and
`h₁ = 1, h₂ = 2` gives the single point `−1`.

### §2's finding — **the brief's step 2 is false, and not narrowly**

```lean
def row3Locus (q : K) (h : ℕ) : Set K := {v | headFactor q v h = 0}
def row4Locus (q : K) (h : ℕ) : Set K := {v | psiNum   q v h = 0}

theorem row3Locus_eq_boundaryLocus [IsDomain K] (hq : q ≠ 1) : row3Locus q h = boundaryLocus K h
theorem mem_row4Locus_iff : v ∈ row4Locus q h ↔ q * v ^ h = 1
theorem row3Locus_disjoint_row4Locus [IsDomain K] (hq : q ≠ 1) :
    Disjoint (row3Locus q h) (row4Locus q h)
```

The brief asked me to *"prove that row 3's zero set and row 4's `τ`-firing locus are the same
set"* and called it *"trivial as algebra."*  They are not the same set: they are **disjoint**,
uniformly in `h`.  Row 3 vanishes on the `h`-th roots of unity `{v^h = 1}`; row 4 vanishes on
the `q`-shifted locus `{q v^h = 1}`; a common point gives `q = 1`, the one value every
statement here excludes.  One line of algebra, and it points the opposite way.

And row 4 fails even the *nesting* that makes row 3's difference a single point:

```lean
theorem row4Locus_one_not_subset_two {F} [Field F] {q : F} (hq0 : q ≠ 0) (hq : q ≠ 1) :
    ¬ row4Locus q 1 ⊆ row4Locus q 2
theorem control_row4_not_nested_rat :     -- over ℚ at q = 4
    (1/4 : ℚ) ∈ row4Locus 4 1 ∧ (1/4 : ℚ) ∉ row4Locus 4 2 ∧
    (1/2 : ℚ) ∈ row4Locus 4 2 ∧ (1/2 : ℚ) ∉ row4Locus 4 1
```

Row 3's type-`B` locus always sits **inside** its type-`C` locus, so "the difference is one
point" is a meaningful question there.  Row 4's two loci are **incomparable**: `{q⁻¹}` versus
the square roots of `q⁻¹`, which meet only at `q = 1`.  There is no one-point difference to
measure, because there is no difference of nested sets.

**Verdict on Q534's falsifier, as far as this arm reaches.**  Q534's falsifier reads: *"If the
difference in any of rows 1–4 persists after restricting away from `v = −1` / `w = −1`, two of
the five are genuinely independent."*  Row 4's `B`/`C` difference **does** persist away from
`−1` — over ℚ at `q = 4` the two loci differ at `1/4` and at `±1/2`, none of which is `−1`.
On the *zero-locus* reading of "difference", **the falsifier fires**.  On the *substitution*
reading it does not fire at all, because the substitution is literally shared (`rfl`).

So Q534's claim needs splitting, and this is the session's main result:

> **Rows 3 and 4 share the substitution `v ↦ v^h` exactly, and share no zero.**
> "All five differences are one substitution" is **true** at the level of the substitution and
> **false** at the level of the locus, for the two rows where I own both sides.  The
> `v = −1` localisation is a property of row 3's locus — of `{v^h = 1}` specifically, where
> the `h`-families nest — and does not transfer to row 4 merely because row 4 is also an
> `h`-family.

### §3's actual question — **is `h = 1` vs `2` special, or the smallest instance?**

Asked by the brief; answered **both, and the "special" part is the sharper one.**

```lean
theorem control_h3_not_localised :
    boundaryLocus (ZMod 7) 3 \ boundaryLocus (ZMod 7) 1 = ({2, 4} : Set (ZMod 7))
theorem control_h2_localised_zmod7 :
    boundaryLocus (ZMod 7) 2 \ boundaryLocus (ZMod 7) 1 = ({6} : Set (ZMod 7))
```

`h = 2` is the smallest instance of the divisibility family — and it is the **only** exponent
whose difference from `h = 1` is a single point *in every* domain with `−1 ≠ 1`, because the
only universally-available non-trivial root of unity is `−1`.  For `h > 2` the size of the
difference is a property of the **ambient field**, not of `h`: over `ZMod 7` it has two
points, over `ℚ` it is `∅` for every odd `h` and `{−1}` for every even `h` (which is
`HeadExclusionRootsOfUnity.pow_eq_pow_mod_two`, already proved 10-19).  The two `ZMod 7`
theorems are each other's positive control: the same ambient, the same shape, `h` the only
thing that moved, and the answer changed.

⚠ Q538 (allocated by BROWSE c2 today) asks whether the surviving parity is **mod 4, not mod
2**.  Nothing here bears on that — `pow_eq_pow_mod_two` is a statement about the `τ`-fixed
locus `{v² = 1}`, where mod 2 is forced by the defining equation, not a claim that mod 2 is
the right invariant anywhere else.  Not pursued; this is a Lean session.

### Row 4 as the paper writes it

```lean
def psi (q v : F) (h : ℕ) : F := (1 - q * v ^ h) / (1 - q)          -- [Field F]
theorem psi_one (q v : F) : psi q v 1 = (1 - q * v) / (1 - q)       -- cor:rhoB L523
theorem psi_eq_zero_iff (hq : q ≠ 1) : psi q v h = 0 ↔ psiNum q v h = 0
theorem psi_at_one (hq : q ≠ 1) (h : ℕ) : psi q 1 h = 1
theorem psi_two_div_psi_one (hq : q ≠ 1) :
    psi q v 2 / psi q v 1 = (1 - q * v ^ 2) / (1 - q * v)            -- the cor:ratio shape
```

`psi_two_div_psi_one` is the Lean form of `cor:ratio` (L635–L648) at a single slot: the
`(1−q)` cancels, so the type-dependence of row 4 is exactly the ratio of the two numerators —
the boxed `u`-diagonal block ratio.  `psi_at_one` records why all three plants of
`rem:plants` coincide with the truth at `+1` (L973): `ψ` normalises to `1` there in **every**
type.

---

## 2. Lesion register

Pre-registered **before** the build, as required: two irremovable, one removable.

| hypothesis | declaration | pre-reg | **Lean's verdict** | witness |
|---|---|---|---|---|
| `q ≠ 1` | `row3Locus_disjoint_row4Locus` | irremovable | **IRREMOVABLE** | `lesion_hq_disjoint`: at `q = 1`, `1` lies in both loci |
| `(-1 : K) ≠ 1` | `boundaryLocus_two_sdiff_one` | irremovable | **IRREMOVABLE** | `lesion_hchar_charTwo`: over `ZMod 2` the difference is `∅`, not `{-1}` |
| `(-1 : K) ≠ 1` | `rows_share_substitution` | removable | **ABSENT** (never written; holds over `ZMod 2`) | `lesion_hchar_removable_from_unification` |
| `[IsDomain K]`, `q ≠ 0` | `row4Locus_one` | *not pre-registered* | **REMOVABLE — deleted by Lean** | the `unusedVariables` linter flagged `hq0`; both binders deleted, statement unchanged |
| `1 - q*v ≠ 0` | `psi_two_div_psi_one` | *not pre-registered* | **REMOVABLE — deleted by Lean** | junk-value division makes both sides `0` when `1 - q*v = 0`, so the identity holds unconditionally |
| `q ≠ 1` | `psi_eq_zero_iff` | — | **RELOCATED** | it is `ψ_X`'s own well-definedness condition at abstract L45, which divides by `1 − q`; the paper never writes it |
| `[IsDomain K]` | `row3Locus_eq_boundaryLocus` | — | **RELOCATED** (upstream) | graded 10-19 on `fires_iff_pow_eq_one` (`lesion_isDomain`, `ZMod 6`); not re-graded here |

**The measuring happened in the unpre-registered rows.**  All three pre-registrations came
back confirmed, which by the brief's own standard means they measured nothing.  The two
hypotheses that actually moved were ones I wrote into statements without thinking and **Lean
deleted** — `q ≠ 0` and `[IsDomain K]` from `row4Locus_one`, and `1 − q*v ≠ 0` from
`psi_two_div_psi_one`.  Both were put there by analogy with neighbouring theorems rather than
by need.  That is a better detector than my pre-registrations were, and it is free: the
`unusedVariables` linter is a lesion pass I get on every build and had not been reading as
one.

The `q ≠ 1` / RELOCATED row is a finding about the **paper**, not about Lean: the gauge
character `ψ_X` is defined by a division whose denominator the paper never excludes from
vanishing.  Every `ψ`-statement in the source silently assumes `q ≠ 1`.

---

## 3. Declaration-set comparison, both directions

Enumerated from the **environment** (`run_cmd` over `env.constants`, prefix
`TypeDependenceLocus`) and from the **source** (top-level keyword grep), then compared as
sets, both directions:

```
ENV  = 51      SRC = 34      34 + 17 = 51 ✓
in SRC not in ENV:  (empty)
in ENV not in SRC:  17  —  12 × `<thm>._simp_1_k`  +  5 × `<def>.eq_1`
                           (psi.eq_1, psiNum.eq_1, row3Locus.eq_1, row4Locus.eq_1,
                            and the 12 simp-lemma children of the four `decide`/`simp only`
                            set proofs)
```

Symmetric difference closed in both directions.  No anonymous `example` in the file (the
10-19 under-count mode), and no wrapped docstring line beginning with a declaration keyword
(the 10-18 over-count mode) — checked by the empty `SRC \ ENV` column, which is the direction
that detects both.

## 4. `#print axioms`

```
rows_share_substitution          does not depend on any axioms
boundaryLocus_subset_of_dvd      does not depend on any axioms
row3Locus_disjoint_row4Locus     [propext, Classical.choice, Quot.sound]
boundaryLocus_two_sdiff_one      [propext, Classical.choice, Quot.sound]
boundaryLocus_symmDiff_of_dvd    [propext, Classical.choice, Quot.sound]
row4Locus_one_not_subset_two     [propext, Classical.choice, Quot.sound]
psi_eq_zero_iff                  [propext, Classical.choice, Quot.sound]
psi_two_div_psi_one              [propext, Classical.choice, Quot.sound]
control_h3_not_localised         [propext, Classical.choice, Quot.sound]
lesion_hchar_charTwo             [propext, Classical.choice, Quot.sound]
lesion_hq_disjoint               [propext, Classical.choice, Quot.sound]
control_row4_not_nested_rat      [propext, Classical.choice, Quot.sound]
```

Standard three throughout; the two unification/divisibility statements are axiom-free.

---

## 5. §scope-honesty — what is **not** formalised

Row 3's head factor and row 4's `ψ` are abstract polynomials in **one variable `v` over a
commutative ring**.  There is **no `w`, no `q`-deformed weight function, no `Ã`, no `W^X`, no
`E^X`, no fixed-point set, no `LG(N,2N)`, and no root system** anywhere in this module.  `h`
is an exponent and nothing else.

- Not formalised: that `h_B = 1` and `h_C = 2`.  That is a root-system input (`def:D`),
  carried in as the numerals `1` and `2`.
- Not formalised: that row 3's factor **is** the head factor of `T^X_σ`, or that row 4's `ψ`
  **is** the gauge making `Ã^B`'s collapse multiplier `1` (`cor:rhoB`).  Both are asserted by
  the paper and imported here as definitions.
- Not formalised: rows **1, 2 and 5** of Q534's table (Motegi's `w`-blocks L1057/L1366;
  `normalizationWB`'s `2t_k` vs `t_k`; Naruse's YBE relations, which are at `agent-summary`
  provenance anyway).  This arm settles 3-vs-4 and says nothing about 1, 2, 5.
- Not re-run: the `0/18`-vs-`6/6` census of `rem:measured` (L359–L364).
- ⚠ Motegi's `h` (`RTVtheorem`, L4495) is the character scaling the cotangent fibres — a
  **different object sharing the letter**.  It does not appear in this file, no binder is
  named for it, and nothing here identifies the two.  Pre-registered as a name collision in
  the brief; kept as one.

**The disjointness result is therefore about the abstract polynomials, not about the weight
functions.**  If the right notion of "row 4's locus" is something other than `{ψ_X = 0}` —
poles, say, or the locus where `ψ_C/ψ_B = 1`, which is `{0, 1}` and also not `{−1}` — then
the falsifier reading of §2 changes and the unification reading does not.  I checked the
ratio locus by hand and it is `{0, 1}` in both rows; I did not formalise it.

---

## 6. For the next PROVE slot

PROVE c2 is on G5 and on `eq:PhiX-is-RTV` (Q536: all type-dependence in the row factor
`D^X_r`).  If that holds, it is a statement about **where** the dependence lives, and this
module is a statement about **what shape** it has.  They are compatible and they are not the
same claim.  The thing to carry across: *an `h`-family of shapes and an `h`-family of loci are
different objects, and nesting is a property of the second only.*
