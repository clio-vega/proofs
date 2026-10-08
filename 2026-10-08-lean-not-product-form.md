# LEAN 2026-10-08 c1 — `not_product_form`: a root in `(-1,0)` kills the product form

**Target:** `TworowD4Kernel.ProductForm.not_product_form`
**Project:** `/home/clio/projects/lean/tworow_d4_kernel` → `clio-vega/tworow-d4-kernel`
**New file:** `TworowD4Kernel/ProductFormObstruction.lean` (142 lines, 6 declarations)
**Paper:** Theorem D, `proofs/2026-10-07-two-part-green-polynomials.tex`

## Result: sorry-free, 6/6 declarations

```lean
theorem not_product_form (b : ℕ) (hb : 3 ≤ b) (hodd : Odd b) :
    ¬ ∃ (c : ℕ) (d : Multiset ℕ), (∀ k ∈ d, 1 ≤ k) ∧
      ∀ t : ℝ, t ^ b - t ^ (b - 1) + 1 = t ^ c * (d.map (fun k => 1 - t ^ k)).prod
```

`#print axioms`, verbatim:

```
'TworowD4Kernel.ProductForm.one_sub_pow_pos' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.one_sub_pow_ne_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.prod_one_sub_pow_ne_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_product_form' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.prod_eq_zero_of_exponent_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.ProductForm.not_prod_ne_zero_without_hd' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Standard three only. No `sorryAx`, no `native_decide`. **Zero sorries** — nothing is bookmarked.
Full root build (`lake build TworowD4Kernel`) exit 0, 3204 jobs.

## The finding: cyclotomic theory was a detour

Yesterday's `NonCyclotomicRoot.lean` left the consequence open and named the route in its own
docstring — `Polynomial.isRoot_cyclotomic_iff`, `IsPrimitiveRoot`, and the input *"every root of a
cyclotomic polynomial has modulus 1"*.

**None of it is needed.** Against the *specific* shape `t^c ∏(1 - t^{d_i})` the obstruction is an
elementary sign argument with no roots of unity anywhere in it:

> For `t ∈ (-1,0)`, `|t| < 1` so `|t^k| = |t|^k < 1` for every `k ≥ 1`, hence `1 - t^k > 0`.
> And `t ≠ 0`, so `t^c ≠ 0`. A finite product of nonzero reals is nonzero — but `D_b` *vanishes*
> there. So `D_b` is not of that shape.

Four lines of Lean per step. The cyclotomic route would prove something more general than Theorem D
needs; the paper asked for the modulus-1 theorem because it was thinking about *all* products of
cyclotomics, not about this product.

## What is NOT formalised — please read this part

1. **The `±1` exponents.** Theorem D's display allows factors in the *denominator*. I proved the
   numerator-only shape. Clearing denominators is the same argument on paper — the left side still
   vanishes at the root, the right side still doesn't — but **that restatement is not in the Lean
   file**, and `not_product_form` should not be read as covering it. This is the one real gap
   between my declaration and Theorem D's display, and it is mine, not the paper's.
2. `D_{a,b} = Y^{(a,b)}_{(a,b)} = t^b - t^{b-1} + 1` is paper-side (Theorem C at `m = b`). In Lean
   the polynomial is simply *written down*.
3. The quantification over all `ℓ(λ) ≤ 2` and all two-part `ρ`.

`Y`, Hall–Littlewood `P`, Kostka–Foulkes and charge have **no Lean definitions in this project**.
`unproved ≠ unformalised`. The registry parent `thm-D-product-form-obstruction` therefore stays
`proved`; only the new child is `lean-verified`.

## Instruments

`declaration uses 'sorry'` is **dead** in this repo (reads 0 with a sorry planted), so the grade
rests on `#print axioms`, run as a two-arm canary:

| arm | `sorryAx` count |
|---|---|
| clean | **0** |
| sorry planted in `one_sub_pow_pos` | **4** |

The 4 are exactly `one_sub_pow_pos → one_sub_pow_ne_zero → prod_one_sub_pow_ne_zero →
not_product_form`, and the 2 independent declarations were correctly spared — so the propagation
pattern doubles as a check that the dependency structure is what I think it is.

**`lake build` exited 0 *with* the sorry planted.** The exit code is not the instrument.

Also: my first overlap grep used `--include='*.lean'`, which made `ugrep` emit *"No such file or
directory"*. Counts happened to be right, but I re-ran via `find | xargs` with a positive control
(fires on 2 files) and a negative control (reads 0) before trusting any of them.

## Ablations — four, each guarded by an occurrence count before mutating

| ablation | expected | got |
|---|---|---|
| delete `hd : ∀ k ∈ d, 1 ≤ k` | fail | **exit 1**, `Unknown identifier 'hd'` |
| rename the `exists_root_Ioo` citation | fail | **exit 1** |
| **delete that step outright** | fail | **exit 1**, `Unknown identifier 'ht₀'` |
| no-op whitespace in an unrelated docstring | **build** | **exit 0** |

The third is there because renaming tests that a *name is referenced*, not that a *step is
necessary*. The fourth is there so the harness can be shown to distinguish the two — without it,
three failures prove only that I can break a build.

Better than any of them: **the hypothesis is refuted as a statement, not just as a proof.**
`not_prod_ne_zero_without_hd` shows that with `hd` dropped the conclusion is **false** — at
`d = {0}` the factor is `1 - t^0 = 0`, so the product vanishes for every `c` and every `t`, with the
interval hypothesis never consulted. A build failure is a fact about my proof; this is a fact about
the statement. It sits beside `not_exists_root_Ioo_four` from yesterday.

## The overlap guard found a real near-miss

I searched the statement *shape*, not my intended name. `PhiNonvanishing.prod_factors_ne_zero`
already proves `∏(1 - X^{e_i}) ≠ 0` over a `List`, with the **same** hypothesis `∀ m ∈ e, 1 ≤ m`,
and for the same reason (`k = 0` makes the factor vanish).

It is a **different statement and does not imply mine**: it says the *polynomial* is nonzero in
`ℤ[X]`, seen by evaluating at `0`; I need nonvanishing *at one particular real point*, and a nonzero
polynomial may certainly vanish at a point. The implication runs the other way. Worth knowing the
cousin is there — had I searched for my intended name instead of the shape, I would not have found
it, and I would not have been able to say which way the implication runs.

I used `Multiset` (exponents carry no order) and recorded in the file docstring why, given the
`List` precedent next door.

## Grades, naming the tool

- `trustcheck.py ... validate --files-dir .` → **`OK`**
- `registry_validate.py --proofs-dir .` → **2 problems**

The two disagree, as they do about `peer-claimed`. Both problems are on
`rick-two-point-formula-thm25`, which I never touched; running the *same* validator on the pre-edit
backup also reports **2**, so my edit introduced **zero** new problems.

`trustcheck` grades the tree, not the leaf — it will say `OK` with a nonexistent declaration in a
`lean` field. My evidence that `TworowD4Kernel.ProductForm.not_product_form` is a real name is that
**Lean itself accepted it** in the `#print axioms` output pasted above.

## For Robin — the one question

Is the `±1`-exponent version worth formalising, or is the numerator-only case the honest statement
of the obstruction? Clearing denominators is routine on paper, but it would roughly double the Lean
file, and I would rather you told me whether Theorem D's display is load-bearing for the shape
question before I spend a session on it.
