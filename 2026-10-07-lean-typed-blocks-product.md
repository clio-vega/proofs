# Lean session, 2026-10-07 c1 — the type-split product of labelled-block counts

**Project:** `/home/clio/projects/lean/tworow_d4_kernel` (Lean 4.30.0 / Mathlib v4.30.0)
**New file:** `TworowD4Kernel/TypedBlocks.lean` (248 lines, 23 declaration blocks — 22 named
plus one anonymous `Decidable` instance)
**Registry node:** `typed-blocks-product-lean`, a sibling of `labelled-blocks-multinomial-lean`
under `theorem-A-existence-trivial` in `proofs/registry/ls-corollary-4-coupling.json`
**Paper proof:** `proofs/2026-10-06-ls-corollary-4-coupling.tex`, Theorem A
(`arXiv:math/0202090`, Lenart–Sottile, *Skew Schubert polynomials*, Corollary 4 and the
remark following Theorem 2)

---

## The target declaration

```lean
theorem TworowD4Kernel.TypedBlocks.card_typedBlockFunctions_eq_prod_multinomial
    {A V T : Type*} [Fintype A] [DecidableEq A] [Fintype V] [DecidableEq V]
    [Fintype T] [DecidableEq T] {τ : A → T} {m : T → V → ℕ}
    (hm : ∀ s, ∑ v, m s v = Fintype.card (typePart τ s)) :
    #(typedBlockFunctions A τ m) = ∏ s, Nat.multinomial univ (m s)
```

where `typePart τ s = {a : A // τ a = s}`, `restrict τ f s = fun x => f x.1`, and

```lean
TypedBlockSizes τ m f  ↔  ∀ s, LabelledBlocks.BlockSizes (m s) (restrict τ f s)
```

so `typedBlockFunctions A τ m` is the finset of `f : A → V` that cut **each** type part
into blocks of the sizes prescribed for that type.

The mathematical content is that the per-type conditions are **independent** — they
constrain disjoint parts of the data. That is isolated as

```lean
def typedBlockEquivPi (τ : A → T) (m : T → V → ℕ) :
    {f : A → V // TypedBlockSizes τ m f} ≃ ∀ s, {g : typePart τ s → V // BlockSizes (m s) g}
```

## What builds sorry-free

**Everything. There are no sorries.** Full-tree `lake build` green (3200 jobs);
`lake test` exit 0.

| declaration | what it is |
|---|---|
| `card_typedBlockFunctions_eq_prod_multinomial` | the target |
| `typedBlockEquivPi` | the types are independent |
| `restrictEquiv`, `restrictEquiv_apply` | restriction to type parts, = Mathlib's `Equiv.piCongrFiberwise` |
| `card_typedBlockFunctions_eq_zero_of_sum_ne` | **necessity** of the hypothesis |
| `prod_multinomial_pos`, `typedBlockFunctions_nonempty` | existence half of Theorem A at this abstraction |
| `sum_card_fibre` | the blocks of `g : B → V` partition `B` |
| `tauEx`, `mEx`, `sum_mEx`, `card_typedBlockFunctions_ex`, `prod_multinomial_ex` | positive case, by enumeration |
| `mBadEx`, `sum_mBadEx_ne`, `card_typedBlockFunctions_badEx`, `prod_multinomial_badEx` | negative control |

## `#print axioms`

All twelve audited declarations, including the main one, report exactly the standard three:

```
'TworowD4Kernel.TypedBlocks.card_typedBlockFunctions_eq_prod_multinomial' depends on axioms:
  [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.typedBlockEquivPi'                            [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.card_typedBlockFunctions_eq_zero_of_sum_ne'   [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.typedBlockFunctions_nonempty'                 [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.sum_card_fibre'                               [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.prod_multinomial_pos'                         [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.card_typedBlockFunctions_ex'                  [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.prod_multinomial_ex'                          [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.card_typedBlockFunctions_badEx'               [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.prod_multinomial_badEx'                       [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.sum_mEx'                                      [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.TypedBlocks.sum_mBadEx_ne'                                [propext, Classical.choice, Quot.sound]
```

No `native_decide` anywhere (its `Lean.ofReduceBool` is outside the allowlist); the four
`by decide` witnesses are kernel reductions.

---

## Mathlib pre-flight

Searched by **declaration-block parser**, not line-based grep — the parser splits a file
into blocks starting at a declaration keyword and matches regexes against the whole
block, so a statement spread over several lines is found. Positively controlled on
`Multiset.countPerms_filter_ne`, which is exactly the declaration a line-based grep missed
in the 2026-10-06 session (`card` and `multinomial` sit on different lines): the parser
finds it.

**Found, and used as such — the transport is Mathlib's:**

- `Equiv.piCongrFiberwise` (`Mathlib/Logic/Equiv/Basic.lean:963`) — splits `A → V` into
  `∀ s, (typePart τ s → V)` along the type map
- `Equiv.subtypePiEquivPi` (`Mathlib/Logic/Equiv/Basic.lean:412`) — pushes a pointwise
  predicate through a `Pi`
- `Fintype.card_pi` (`Mathlib/Data/Fintype/BigOperators.lean:132`) — `Pi` of subtypes to
  a product of cardinalities

**Not found — the count is not there.** No Mathlib declaration block states `card` of a
fibrewise-constrained function set as a product. `Multiset.bell` and `Multiset.countPerms`
are still arithmetic definitions under the standing TODO in
`Mathlib/Combinatorics/Enumerative/Bell.lean` (*"Prove that it actually counts the number
of partitions as indicated"*), and `bell_mul_eq` is pure arithmetic. So the off-ramp in
the brief — *if Mathlib has the whole thing, stop and report that* — did not apply: the
multinomial half is the sibling node's, the product half is this file's, and only the
plumbing between them is Mathlib's. The file's docstring says which half is which.

**Repo pre-flight:** 0 declaration blocks of this statement shape existed
(`piCongrFiberwise` / `subtypePiEquivPi` / `card_pi` / `multinomial` / `∏`-with-`card`).

---

## Verification

### Sorry scan — and the canary killed the obvious instrument

The natural scan is to grep the build output for `declaration uses 'sorry'`. **It never
fires in this toolchain.** With a `sorry` planted in `prod_multinomial_pos` the warning
count read **0**, exactly as it did with a clean file. Had I not canary-validated, I would
have reported a green scan that could not have been anything else.

The two instruments that *do* fire, both validated in both directions over the whole tree:

| instrument | sorry planted | restored |
|---|---|---|
| `sorryAx` in full `#print axioms` build output | **2** hits | **0** |
| comment-stripped grep (`/tmp/sorryscan.py`) | **1** hit | **0** |

The `sorryAx` count of 2 is itself informative: it propagated to
`typedBlockFunctions_nonempty`, which depends on the sorried lemma. A *plain* grep is
useless here — 12 docstrings in this repo mention the word "sorry" in prose, which is why
the scan strips comments before matching.

File `sha256` after restore: `b45cd49b5664073189cee3d7533197e6c5c0f43f148f4212edfdc8a99c03e0e6`,
**byte-identical** to before planting.

### Ablation, in two strengths

Every mutation asserted its anchor present with **count == 1** before mutating (an
ablation that does not ablate is indistinguishable from an absent dependency), and read
the exit via `PIPESTATUS[0]` (`lake` is not on PATH by default; exit 127 is otherwise
invisible).

*Rename ablation* — all four give ablated exit **1**, restored exit **0**:
`Fintype.card_pi`, `Equiv.subtypePiEquivPi`, `Equiv.piCongrFiberwise`,
`LabelledBlocks.card_blockFunctions_eq_multinomial`.

But renaming only proves the name is **referenced**, not that the step is **necessary**.
So, *deletion ablation*:

| mutation | exit | error |
|---|---|---|
| remove the `Fintype.card_pi` step from the `rw` chain | 1 | `typeclass instance problem is stuck` |
| remove the transport along `typedBlockEquivPi` | 1 | `rewrite failed: Did not find an occurrence of the pattern` |
| *(control)* append a no-op `<;> skip` | **0** | — harness correctly reported "not necessary" |

The third row is the point: the harness can tell a real dependency from a mutation that
changes nothing.

### Negative control — it fires

`mBadEx` prescribes sizes `(2,0)` on a type part with **3** elements. The count collapses
to **0** while the product of multinomials stays at **2**. So the summation hypothesis is
load-bearing: without it the stated equality is **false**, not merely unproved. This is
recorded as a theorem (`card_typedBlockFunctions_eq_zero_of_sum_ne`, plus the three
`decide` witnesses), not as a shell transcript.

### Positive case by a different mechanism

`A = Fin 5` split `3 + 2` by `tauEx`, sizes `(2,1)` and `(1,1)`. Predicted
`3!/(2!·1!) · 2!/(1!·1!) = 3 · 2 = 6`. Kernel enumeration of all **32** functions
`Fin 5 → Fin 2` gives **6** — not the proof's route, which goes through `typedBlockEquivPi`
and the multinomial theorem. The product `∏ s, Nat.multinomial univ (mEx s) = 6` is
evaluated separately.

All witnesses are also `#guard` shadows in `TworowD4KernelTests.lean`, which is outside
`defaultTargets` precisely so it can fail alone.

I did **not** prove the same `Prop` twice and call it corroboration — definitional proof
irrelevance makes that comparison constant.

---

## Side repair: a docstring that promised a declaration nobody wrote

`LabelledBlocks.lean`'s module docstring read:

> *"The product formula itself is recorded as `multinomial_prod_pow`, the form in which
> Theorem A uses it: a product of independent multinomials with repeated block sizes."*

`multinomial_prod_pow` **does not exist** — not in that file, not in the repo, not in
`projects/proofs/`. The 2026-10-06 session's docstring named a declaration it never wrote,
and the brief for this session quoted that sentence back to me as though it described
something on disk. A docstring is a claim about contents. The paragraph now points at
`TypedBlocks.card_typedBlockFunctions_eq_prod_multinomial` and records the correction
in place, so the next reader sees that the pointer moved rather than wondering where it went.

---

## Scope — what is still not formalised

**Theorem A itself is not formalised.** The labelled Bruhat order, increasing chains,
types, `Γ_α`, and the numbers `I_α(u,w)`, `c^w_{u,v}` do not exist in Lean. What exists is
the abstract product-over-fibres count. The dictionary to the paper is:

> `A = Γ(u,w)`, `T` = the set of types `α`, `τ = type`, `V = ⨆_v Γ(w₀v, w₀)` over **all**
> types, `m α C = c^w_{u,v}` for `C ∈ Γ_α(w₀v, w₀)` and `m α C = 0` for any block `C`
> whose own type is not `α`.

One thing fell out of writing it that way, and it is the only new observation in the
session: condition (i) of the paper's Definition — that `f` be **type-preserving** — is
*not* a separate hypothesis in the Lean statement. It is **forced** by the zero
prescriptions, since a chain of type `α` cannot land in a block of type `≠ α` when that
fibre is prescribed empty. Blocks with `m s v = 0` contribute `0! = 1` to the denominator,
so `∏_s Nat.multinomial univ (m s)` is exactly `N(u,w)`. Type-preservation is a *shadow of
the size prescription*, not an independent axiom — which, read back into the paper, is
another way of saying what Corollary A already says: the fibre sizes are an input to the
construction, not an output.

But that dictionary is a dictionary **on paper**. It is not a Lean definition, and nothing
on the paper side is promoted by this file. The parent node `theorem-A-existence-trivial`
stays `proved`. `unproved ≠ unformalised`.

## Registry

`typed-blocks-product-lean` added with `trust: lean-verified`; all 11 names in its `lean`
field verified present in the file by declaration-match, not by eye.

`trustcheck.py ... --files-dir .` → exit 0, 0 problems.
`registry_validate.py --proofs-dir .` → the `--proofs-dir proofs` default double-prefixes
to `proofs/proofs/...` and reports **every** node of a valid registry as a broken file
pointer; with `.` it is clean.
