# LEAN 2026-10-18 c2 — the two colength counts of Theorem A, and a retraction

**Target.** The combinatorial skeleton of `lem:colenG`, `lem:colenL`, `lem:typeB` of
`proofs/2026-10-18-atype-divisors-and-the-gkm-lattice.tex` (Theorem A, `L_ij = G_ij`).

**Project.** `projects/lean/tworow_d4_kernel`, module `TworowD4Kernel/AtypeColength.lean`
(added to the root import `TworowD4Kernel.lean` at line 7).

**Status: sorry-free.** 0 `sorryAx`, 0 nonstandard axioms, over **75 environment records /
37 source-level declarations**, both enumerated from the environment.

---

## What builds

| Paper | Lean declaration | Content |
|---|---|---|
| (implicit `2 ≤ N`) | `two_le_of_ne` | `i ≠ j` in `Fin N` **forces** `2 ≤ N` — derived, not assumed. Free strengthening of the paper's implicit hypothesis. |
| `defn` L178–186 | `refl` | `s_ij M = if |M ∩ {i,j}| = 1 then M ∆ {i,j} else M`, stated in the paper's own case form. |
| — | `refl_comm`, `mem_moved_iff`, `refl_eq_of_left`, `refl_involutive`, `refl_mem_moved`, `refl_ne_self` | `s_ij` is an involution on all of `P`, preserves `P^mv`, and is fixed-point-free there. |
| **`lem:colenG`** | **`movedOrbits_card`** | `s_ij` partitions `P^mv` into exactly `2^(N-2)` orbits of size two. Counted **independently** of the label side: `|P^mv| = 2^(N-1)` (as a disjoint union of the two half-labels), then fibrewise over `M ↦ {M, s M}` with every fibre of card exactly `2` (`fiber_card`). |
| **`lem:colenL`** | **`ord_sum`** | `∑_{I ∈ P} ord Ã_{I,I} = #{I : i∉I, j∈I} = 2^(N-2)`, the label count via the bijection `I ↦ I.erase j` onto `(univ \ {i,j}).powerset`. |
| Theorem A, comparison step | **`colength_match`** | `#movedOrbits = ∑_I ord_I`. |
| `lem:typeB` | `ord_sum_of_indicator`, `ord_sum_typeB` | see TASK 2 below. |

**Not in scope, and deliberately so.** The module theory is *not* formalised:
δ-surjectivity, `ker δ = G_ij`, Smith normal form, "colength of the image of an injective
endomorphism over a DVR = valuation of the determinant", and `(M3)+(M4)`
(`det Φ = ±∏_I Ã_{I,I}`). Those are the steps that turn the two cardinalities into two
colengths. What is certified here is exactly the counting, which is what the paper reduces
both colengths *to*.

---

## TASK 1 — the pre-registered lesion fired. My banner was wrong.

`SUMMARY.md` said of Theorem A:

> *"The two counts run over **different index sets** — two-element `s_ij`-orbits on one
> side, labels with `i∉I, j∈I` on the other — and neither is derived from the other, which
> is **the only way the match carries information**."*

**`AtypeColength.orbitLabelEquiv` type-checks**, sorry-free:

```lean
noncomputable def orbitLabelEquiv {i j : Fin N} (h : i ≠ j) :
    {I : Finset (Fin N) // I ∈ labels i j} ≃
      {o : Finset (Finset (Fin N)) // o ∈ movedOrbits i j}
```

`I ↦ {I, s_ij I}`, with inverse *"take the member of the orbit that excludes `i` and
includes `j`"* — injective because `s_ij I` always **re-admits** `i`
(`mem_refl_of_mem_labels`), surjective because the `i`-present half of each moved orbit
reflects into a label (`refl_mem_labels`). The derived corollary
`labels_card_eq_movedOrbits_card` proves the two cardinalities equal **from the `Equiv`
alone** — it cites neither `labels_card` nor `movedOrbits_card`.

So the two index sets are **canonically isomorphic**, and the stronger `SUMMARY.md`
phrasing is **retracted**. Banner corrected in place.

**The honest residue.** The two *routes* to the number remain independent —
δ-surjectivity + module theory on one side, Smith normal form + factor-by-factor valuation
on the other — and that independence is where the information actually lives. But the
*index sets* are not independent. ***A set is itself a projection*** (the 10-17 scar: `x`
and `xψ` gave the same 57 **and the same set**; only the denominators differed).

**The paper was already right.** `rem` at L551 says only that the bijection "is not used
anywhere in the proof" — which is true, and is the correct weaker claim. The error was
introduced in the *summary*, not in the mathematics. That is its own lesson: **the banner
is an artifact with its own failure mode, and it can be strictly stronger than the theorem
it advertises.**

**⚠ Theorem A is untouched.** Theorem A equates **colengths of lattices**; `orbitLabelEquiv`
is about **index sets**. The theorem stands either way. Only the epistemic claim about *why*
the match is informative moved. Do not let this be miscopied as "Theorem A broken."

---

## TASK 2 — type independence, and exactly how little it shows

`ord_sum_of_indicator` isolates the real content: the sum depends on the diagonal **only**
through the indicator `[i ∉ I ∧ j ∈ I]`. `ord_sum_typeB` is then that lemma applied to a
type-`B` order function, with the paper's indicator identity as an explicit **hypothesis**
(`hordB`) rather than smuggled in — type `B` and type `C` are two instances of one
declaration, so there is nothing to re-prove.

**What it shows:** the *combinatorial skeleton* of `lem:colenL` is type-independent.
**What it does not show:** that the weight functions are. The whole type-dependence sits in
`2t_k` vs `t_k` (Motegi `normalizationWB` L1386 vs `normalizationWC` L1078,
arXiv:2610.15344; Q514, closed 10-18).

**The absence of risk is the thing to state.** A numerical engine once produced a green
type-`B` arm while normalising type `B` by `E^C` — a shared wrong premise, structurally
invisible to its control. Lean cannot repeat that here, and the reason is specific: `ord`
**contains no normaliser at all**. So this corollary is evidence about the indicator and
about nothing else. Recorded in the declaration's docstring so a later reader cannot take
it for evidence about the weight functions.

---

## Lesions — graded on the STATEMENT, not the script

Every lesion is graded by **proving the negation of the lesioned statement at a named
instance**, never by deleting a binder and watching `Unknown identifier` (which grades the
script — 5 of 6 lesions on 10-17 did exactly that).

**`i ≠ j` is irremovable from five theorems**, each with an explicit falsifier at
`N = 2, i = j = 0` (all by `decide`):

| Lesion | Falsifier |
|---|---|
| `Lesion.two_le_lesion_fails` | `Fin 1` is inhabited, so `i, j` exist with `¬(2 ≤ 1)`. |
| `Lesion.mem_moved_iff_lesion_fails` | `{0} ∈ moved 0 0` but the RHS disjunction is unsatisfiable at `i = j`. |
| `Lesion.labels_card_lesion_fails` | `labels 0 0 = ∅`, card `0 ≠ 1 = 2^0`. |
| `Lesion.ord_sum_lesion_fails` | the sum is `0 ≠ 1`. |
| `Lesion.movedOrbits_card_lesion_fails` | the count is `≠ 2^0 = 1`, proved by `decide` at `N = 2`. ⚠ My *explanation* of the mechanism — that at `i = j` the map is `M ↦ M.erase i`, injective on the moved locus, giving `2^(N-1)` orbits — is **not verified**: `#eval` on `movedOrbits` does not terminate in budget (nested `Finset (Finset (Fin N))` universes), so it is checked only at `N = 2`, where `2 ≠ 1`. The falsifier is proved; the story about *why* is a conjecture. |
| `Lesion.refl_involutive_lesion_fails` | `s_00 (s_00 {0}) = ∅ ≠ {0}`. |

**And `i ≠ j` is REMOVABLE from `refl_ne_self` — a free strengthening.**
`refl_ne_self_of_mem_moved (i j : Fin N)` drops the binder entirely: on the moved branch
`s_ij M = M ∆ {i,j}`, and `M ∆ T = M` forces `T = ∅`, while `{i,j}` is never empty. So no
point of `P^mv` is fixed, for **any** `i, j` — including `i = j`.

This is the sharp pair the discipline is for: the binder is **irremovable** from
`refl_involutive` and **removable** from `refl_ne_self`, two adjacent theorems about the
same map, and **neither half is visible from the other**. (Twin of the 10-18 `a ≠ b` pair
in `konckiFamily_ne_std` / `konckiFamily_not_diagonal`.)

---

## Controls run

- **Olean deleted, not `touch`ed.** `touch` does not force re-elaboration (lake keys on
  content hash). Removed `AtypeColength.olean` + `.ilean`; mtime moved
  `19:25:38 → 19:26:05`, 2.7 s elaboration. Green is not a replay.
- **Declarations enumerated from the ENVIRONMENT**, via `run_cmd` over `env.constants`, not
  from a source regex (on 10-18 a regex read 28 against the environment's 27: it missed two
  auto-generated `.congr_simp` lemmas and invented one from a wrapped docstring — two faults
  nearly cancelling into a plausible typo). **75 records total**; **37 source-level**;
  **38 auto-generated** (`._simp_*`, `._proof_*`, `.eq_*`). Records split as a *set*, not a
  total.
- **Positive control on the axiom scan.** Planted `theorem plantedSorry : (1:Nat) = 2 := by
  sorry` into the same namespace and re-ran the identical scanner: it reported
  `AtypeColength.plantedSorry :: [sorryAx]`, `RECORDS = 76, sorryAx-carrying = 1,
  nonstandard = [sorryAx]`. The clean arm: `RECORDS = 75, sorryAx-carrying = 0,
  nonstandard = []`. The instrument fires and spares correctly.
- **No `native_decide`** anywhere (it would add `Lean.ofReduceBool`).

## Build status — what I did and did NOT verify

**Verified:** `lake build TworowD4Kernel.AtypeColength` is green, re-elaborated from source
after deleting the olean (not `touch`ed). The module was also imported and fully elaborated
twice more, by the axiom scanner and its planted-sorry control.

**NOT verified: the whole-project build did not finish in this session.** `lake build`
(3471 jobs) reached `[3463/3471]` and was then terminated — it was still elaborating the
root module `TworowD4Kernel.lean`, which imports all 73 modules, when competing `lake`
invocations of mine contended for the build lock. `TworowD4Kernel.olean` is **absent**; the
log contains no `error` and no `Build completed`. So the honest statement is *terminated*,
not *passed* and not *failed*. The project was green at 3469 jobs before this session;
I added one module and one import line.

**What that root build could have caught, and the substitute I ran:** essentially only a
name collision at import time. Checked instead from the **environment**: all **75** records
carry the `AtypeColength.` prefix, so nothing leaks into the root namespace, and exactly
one module in the project declares `namespace AtypeColength`. (I first reached for a source
regex over declaration keywords for this and it gave `35`, a number about indentation
rather than about namespaces — the environment enumeration is the right instrument and the
regex was the wrong one. Same scar as 10-18.)

---

## `#print axioms` — collected from the environment for all 75 records

```
RECORDS = 75   sorryAx-carrying = 0   nonstandard = []
```

Main declarations, each `[propext, Classical.choice, Quot.sound]`:

```
AtypeColength.movedOrbits_card              -- lem:colenG
AtypeColength.ord_sum                       -- lem:colenL
AtypeColength.colength_match                -- Theorem A, comparison step
AtypeColength.orbitLabelEquiv               -- TASK 1 (the retraction)
AtypeColength.labels_card_eq_movedOrbits_card
AtypeColength.ord_sum_typeB                 -- lem:typeB, indicator half
AtypeColength.refl_ne_self_of_mem_moved     -- free strengthening
```

`AtypeColength.two_le_of_ne` and `AtypeColength.ord` use only `[propext, Quot.sound]`.

---

## Q allocation

**None.** No `Q<n>` appears in any identifier in this module — the lesion namespace is
`AtypeColength.Lesion`, a descriptive name. (On 10-18 I named a namespace `Q515Plant` with
515 next-free and no question opened; *a `Q<n>` in an identifier is an allocation*.)
