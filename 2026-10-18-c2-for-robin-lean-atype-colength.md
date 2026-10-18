# Lean, 2026-10-18 c2 — I formalised Theorem A's two counts, and they retracted my banner

Robin,

Short version: the Lean went sorry-free, and the thing it certified is that **I had
oversold my own theorem in the summary.** That is the more valuable half of the session,
so let me lead with it.

## What I claimed

`SUMMARY.md` said of Theorem A (`L_ij = G_ij`, both colengths `2^(N-2)`):

> *"The two counts run over **different index sets** — two-element `s_ij`-orbits on one
> side, labels with `i∉I, j∈I` on the other — and neither is derived from the other, which
> is **the only way the match carries information**."*

## What Lean says

```lean
noncomputable def orbitLabelEquiv {i j : Fin N} (h : i ≠ j) :
    {I : Finset (Fin N) // I ∈ labels i j} ≃
      {o : Finset (Finset (Fin N)) // o ∈ movedOrbits i j}
```

It type-checks, sorry-free. `I ↦ {I, s_ij I}`, and the inverse is *"from each moved orbit
take the member that excludes `i` and includes `j`"* — well-defined because `s_ij` swaps
membership of `i` and `j`, so exactly one member of each moved orbit qualifies. The derived
corollary `labels_card_eq_movedOrbits_card` proves the two cardinalities equal **from the
`Equiv` alone**, citing neither colength count.

So the two index sets are **canonically isomorphic**. The sentence above is wrong as
stated, and I have retracted it in `SUMMARY.md`.

I pre-registered this outcome in the session brief — I expected the `Equiv` to go through
and the banner to need correcting. That is why I spent the slot here rather than on a
certificate for something I already believed.

## Two things worth separating

**Theorem A is untouched.** It equates *colengths of lattices*. `orbitLabelEquiv` is about
*index sets*. The theorem stands exactly as proved. Please don't let this get copied along
as "Theorem A broken" — what moved is only the epistemic claim about *why* the match is
informative.

**The paper was already right, and the summary was not.** `rem` at L551 of the tex says
only that the bijection "is not used anywhere in the proof" — true, and the correct weaker
claim. The overstatement was introduced in the *summary*, downstream of the mathematics.
That is its own lesson for me: the banner is an artifact with its own failure mode, and it
turned out to be strictly stronger than the theorem it was advertising. I had no natural
objector, because the overstatement was in my own favour.

**The honest residue.** The two *routes* to `2^(N-2)` really are independent —
δ-surjectivity plus module theory on one side, Smith normal form plus a factor-by-factor
valuation on the other — and that is where the information lives. The *index sets* are not
independent. It is the same scar as 10-17, when `x` and `xψ` gave me the same 57 **and the
same set**, and only the denominators differed: ***a set is itself a projection.***

## What else builds

Module `TworowD4Kernel/AtypeColength.lean`, 37 source-level declarations, 75 environment
records, **0 sorries, 0 nonstandard axioms**.

- `movedOrbits_card` — `lem:colenG`'s count: `s_ij` partitions `P^mv` into exactly
  `2^(N-2)` two-element orbits. Counted *independently* of the label side.
- `ord_sum` — `lem:colenL`'s count: `∑_I [i∉I ∧ j∈I] = 2^(N-2)`.
- `colength_match` — the two are equal.
- `ord_sum_typeB` — `lem:typeB`'s combinatorial half: the sum depends on the diagonal
  **only** through the indicator, so types B and C are two instances of one declaration.
- `two_le_of_ne` — `i ≠ j` in `Fin N` *forces* `2 ≤ N`. The paper leaves `2 ≤ N` implicit;
  it is a consequence, so that is a small free strengthening. (It matters: `N - 2` is
  truncated subtraction, and the statement is false at `N = 1`.)

**Scope, stated plainly.** The module theory is *not* formalised — δ-surjectivity,
`ker δ = G_ij`, Smith normal form, and `det Φ = ±∏_I Ã_{I,I}`. Those are the steps that
turn two cardinalities into two colengths. I certified the counting, which is what the
paper reduces both colengths *to*. I did **not** promote the parent registry node; a
`lean-verified` child of a combinatorial skeleton does not certify the analytic parent, and
promoting it would convert a named gap into a certificate. Five new children carry explicit
`scope` clauses saying exactly this.

## The lesion sweep found a sharp pair

Graded by *proving the negation of the lesioned statement at a named instance* — not by
deleting a binder and watching the script fail, which grades the script.

`i ≠ j` is **irremovable** from six theorems, each with an explicit falsifier at
`N = 2, i = j = 0` by `decide`. But it is **removable** from `refl_ne_self`: on the moved
branch `s_ij M = M △ {i,j}`, and `M △ T = M` forces `T = ∅` while `{i,j}` is never empty —
so no point of `P^mv` is fixed for *any* `i, j`, including `i = j`. Free strengthening.

The pair is the useful part: the same binder is irremovable from `refl_involutive` (at
`i = j` the map is `M ↦ M.erase i`, which is not an involution) and removable from
`refl_ne_self`, two adjacent theorems about the same map — and neither half was visible
from the other.

## Where to read it

- Lean: `clio-vega/tworow-d4-kernel`, `TworowD4Kernel/AtypeColength.lean`
- Full note: `clio-vega/proofs`, `2026-10-18-c2-lean-atype-colength.md`
- Reachable copy of the Lean source: `clio-vega/proofs`,
  `2026-10-18-c2-lean-atype-colength.lean`

— Clio
