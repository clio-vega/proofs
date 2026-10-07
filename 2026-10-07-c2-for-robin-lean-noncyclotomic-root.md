# The root that carries Theorem D is now machine-checked

2026-10-07, LEAN cycle 2. For Robin — and relevant to Rick, since this is the statement he asked
permission to cite in an FPSAC draft.

## What is verified

```lean
theorem exists_root_Ioo (b : ℕ) (hb : 3 ≤ b) (hodd : Odd b) :
    ∃ t ∈ Set.Ioo (-1 : ℝ) 0, t ^ b - t ^ (b - 1) + 1 = 0
```

Sorry-free, in `TworowD4Kernel/NonCyclotomicRoot.lean` —
`github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/NonCyclotomicRoot.lean`
(commit `c6c5cb1`). 8 declarations, 0 sorries, axioms `[propext, Classical.choice, Quot.sound]` on
every one, no `native_decide`. `lake build` exit 0, `lake test` exit 0.

This is the one step of Theorem D with analytic content, and the only one whose failure would cost
the theorem. Lean follows the paper line for line: `Odd.neg_one_pow` and `Nat.Odd.sub_odd` give
`D_b(−1) = −1 − 1 + 1 = −1`, `zero_pow` twice gives `D_b(0) = 1`, and Mathlib's
`intermediate_value_Ioo` on `Icc (−1) 0` produces the root. No cleverness; that is the point.

## What is *not* verified, said plainly

The Lean file formalises **one arithmetic lemma, not Theorem D.** Three things stay on paper:

1. `D_{a,b} = Y^{(a,b)}_{(a,b)}` — the `m = b` case of Theorem C. In Lean, `t^b − t^{b−1} + 1` is
   simply written down.
2. Root ⇒ not a product of cyclotomics and monomials. This needs "roots of cyclotomic polynomials
   have modulus 1". It was this session's stretch goal and I did not start it.
3. The product-form null itself, quantified over all `λ` with `ℓ(λ) ≤ 2` and all two-part `ρ`.

`Y`, Hall–Littlewood `P`, Kostka–Foulkes and charge have **no Lean definitions in this project**.
A paper↔Lean dictionary is not a Lean definition. So the registry parent
`thm-D-product-form-obstruction` stays `proved`; only the new child is `lean-verified`.

**For Rick specifically:** nothing here weakens the citation, and nothing here strengthens it
beyond the paper. Theorem D needs **one** non-cyclotomic witness family, and odd `b ≥ 3` supplies
infinitely many — so the separate even-`b` observation (`gap:evenb`, *computed* for `3 ≤ b ≤ 15`,
not proved) is not load-bearing and never was.

## The hypothesis is necessary for the theorem, not just for my proof

I ran the ablation twice, because the two arms test different things.

- **Mechanical.** Delete `Odd b` and its single use: `lake build` exit 1, residual goal literally
  `⊢ Odd b`. This shows *my proof route* needs it.
- **Mathematical.** `not_exists_root_Ioo_four` proves the hypothesis-free statement is **false**:
  `b = 4` satisfies `3 ≤ b`, and `t^4 − t^3 + 1` has no root in `(−1,0)` at all, because
  `D_four_pos_of_mem_Icc` shows it is strictly positive on the whole of `[−1,0]`.

The second is the one worth having. A build failure is a fact about me; a counterexample is a fact
about the statement.

## One thing you should know about the registry

I nearly reported this grade on the validator's word. Two arms on `trustcheck … validate`:

- Replace `"lean": "TworowD4Kernel.exists_root_Ioo"` with a **nonexistent** declaration name —
  reads **`OK: … is valid`**.
- Misspell the trust value as `lean-verfied` — reads **`2 problem(s)`**, catching the enum *and*
  the parent boundary rule.

So the validator checks the **shape of the tree** and never resolves the `lean` field: it does not
confirm the declaration exists, let alone that it states what the node claims. Every
`lean-verified` grade in my registry — including the ones from earlier sessions — is corroborated
by `#print axioms` and by nothing else automated. That seems worth either fixing in
`trustcheck.py` or writing down loudly, and I have done the second.

Full write-up:
`github.com/clio-vega/proofs/blob/main/2026-10-07-c2-lean-noncyclotomic-root.md` (`71e5476`).

— Clio
