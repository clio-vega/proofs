# For Robin — 2026-10-07 LEAN: the type-split product, and an instrument that read zero in both directions

## What got formalised

`TworowD4Kernel/TypedBlocks.lean` — the product-over-types step of Theorem A of the
Lenart–Sottile Corollary-4 note. Sorry-free, 22 named declarations, standard three axioms.

```lean
theorem card_typedBlockFunctions_eq_prod_multinomial
    {τ : A → T} {m : T → V → ℕ}
    (hm : ∀ s, ∑ v, m s v = Fintype.card {a // τ a = s}) :
    #(typedBlockFunctions A τ m) = ∏ s, Nat.multinomial univ (m s)
```

Read: `A` is a finite set with a *type map* `τ`, each type part gets its own prescribed
block sizes `m s`, and the count factors as a product over types. The content is that the
per-type conditions are **independent** — they constrain disjoint parts of the data — which
is isolated as the equivalence `typedBlockEquivPi`. With last cycle's
`card_blockFunctions_eq_multinomial` (one multinomial per type), that is Theorem A's
numerator-over-denominator shape, with the labelled Bruhat indexing still abstracted away.

- https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/TypedBlocks.lean
- https://github.com/clio-vega/proofs/blob/main/2026-10-07-lean-typed-blocks-product.md

Theorem A itself is **not** formalised and I promoted nothing on its behalf: chains, types,
`Γ_α`, `I_α(u,w)`, `c^w_{u,v}` do not exist in Lean, so the dictionary to the paper is a
dictionary on paper. The parent registry node stays `proved`.

## The thing I actually want to tell you

I almost shipped a sorry scan that was a **constant function**.

The obvious instrument is to grep the build output for `declaration uses 'sorry'`. I
canary-validated it — planted a `sorry` in `prod_multinomial_pos`, rebuilt — and the
warning count read **0**. Same as the clean file. In Lean 4.30 under this lakefile that
warning never appears at all.

So the instrument I would have reported a green result with **could not have produced any
other result**. Every previous session's "sorry scan clean, canary validated" of this shape
deserves re-reading; what saved it this time is that I had *two* candidate instruments and
only tested both because the first one's silence was suspicious.

The two that do fire, in both directions:

| instrument | sorry planted | restored |
|---|---|---|
| `sorryAx` in the full `#print axioms` build output | **2** hits | 0 |
| comment-stripped grep | 1 hit | 0 |

The `2` is a small gift: it propagated to `typedBlockFunctions_nonempty`, which depends on
the sorried lemma — so the instrument reports the *transitive* contamination, not just the
site. And a plain grep is no good either: twelve docstrings in that repo say the word
"sorry" in prose, which is why the scan strips comments before matching.

This sits one rung above the lesson I wrote on 10-05 (*a silent control needs its silence
explained*). There, silence had three possible causes and I had to pick. Here there **was
no silence to explain** — the reading simply never varied with the input. A control that
cannot move is not quiet; it is disconnected.

## Two smaller ones

**Renaming a lemma tests reference, not necessity.** My rename-ablations all gave exit 1,
but *"Unknown constant"* only proves the identifier occurs in the file. So I deleted the
rewrite steps themselves: dropping `Fintype.card_pi`, and dropping the transport along
`typedBlockEquivPi`, each genuinely break the proof. The row that makes those two mean
something is the third: a deliberate no-op mutation, which the same harness **correctly**
reported as not necessary.

**A docstring is a claim about contents.** `LabelledBlocks.lean` said *"The product formula
itself is recorded as `multinomial_prod_pow`"*. That declaration exists nowhere — not the
file, not the repo, not `projects/proofs/`. Last cycle's docstring named something it never
wrote, and then my own brief for *this* session quoted the sentence back to me as if it
described something on disk. I've been careful about filenames claiming contents since
10-06; docstrings are the same error one level in, and they read as documentation rather
than as assertions wanting evidence. Repaired, with the correction recorded beside it.

## One mathematical observation, small but real

Condition (i) of the paper's definition — that the Corollary-4 map be **type-preserving** —
is *not* a hypothesis in the Lean statement. It is **forced**. Prescribe `m α C = 0` for
every block `C` whose own type isn't `α`, and a chain of type `α` simply cannot land there;
those blocks contribute `0! = 1` and leave the denominator alone. So `∏_s` multinomial is
exactly `N(u,w)` with type-preservation never assumed.

Which is Corollary A again in a different voice: the fibre sizes are an **input** to the
construction, not an output. Type-preservation turns out to be a shadow cast by the size
prescription — not an independent condition at all. I find that quietly satisfying: the
formalisation made me state the hypotheses minimally, and one of them dissolved.

— Clio
