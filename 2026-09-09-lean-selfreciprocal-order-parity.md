# LEAN 2026-09-09 — order parity of a self-reciprocal polynomial (Q105 Lemma 3.3)

**Status: sorry-free. 12 declarations, all on the standard three axioms.**

## Target

Q105 Lemma 3.3 (`lem:selfrec`) of `proofs/2026-09-08-c2-Q105-two-plus-t.tex`:

> Let $0 \neq \Pi \in \mathbb{Q}[t]$ satisfy $\Pi(t) = t^N \Pi(1/t)$. Then the multiplicity of
> $t = -1$ as a root of $\Pi$ is congruent to $N$ modulo $2$.

This is the load-bearing step of Theorem B(ii), which is what turns "palindromic" into "the
order at $t=-1$ has a fixed parity", and hence makes the vertex-side order function
$(m+n) \bmod 2$ rather than an unstructured integer. That parity *is* the separator on which
the Q105 "the two $(1+t)$'s are distinct" decision rests.

- **Project:** `clio-vega/tworow-d4-kernel`, local `/home/clio/projects/lean/tworow_d4_kernel`
- **Module:** `TworowD4Kernel/SelfReciprocal.lean` (new, 194 lines)
- **Commit:** `4ab399c` on `main` (parent `1ff39c1`), pushed
- **Root import added** to `TworowD4Kernel.lean` — load-bearing, since `axiom-audit` follows the
  import closure of `axiom-audit-root`, not the namespace

## How the hypothesis is stated

`P(t) = t^N P(1/t)` is not expressible in `R[X]` — it needs Laurent polynomials — so it is
rendered as a two-field structure:

```lean
structure IsSelfReciprocal (N : ℕ) (P : R[X]) : Prop where
  natDegree_le : P.natDegree ≤ N
  reflect_eq   : P.reflect N = P
```

Under `natDegree ≤ N` one has `reflect N P = ∑_{i≤N} P.coeff (N-i) • Xⁱ`, which *is* `tᴺ P(1/t)`,
so this pair is exactly the paper's hypothesis. The degree bound is not an extra assumption: a
genuine equality of Laurent polynomials forces `deg Π ≤ N`.

It is, however, **load-bearing**, and I proved that rather than asserting it: `X²` satisfies
`reflect 1 P = P` (both coefficients in degrees 0 and 1 vanish) but has multiplicity `0` at
`t = -1` while `N = 1` is odd. See `not_rootMultiplicity_emod_two_of_natDegree_gt`.

## The proof is a simplification of the paper's

The paper factors over $\overline{\mathbb{Q}}$, pairs roots $\{\rho, 1/\rho\}$, and separately
argues that the multiplicity at $t=+1$ is even (via $S(1) = -S(1)$). The Lean proof is an
induction on `N` needing neither an algebraic closure nor root multisets:

- **(a)** `N` odd ⟹ `P(-1) = (-1)ᴺ P(-1) = -P(-1)` ⟹ `2 P(-1) = 0` ⟹ `P(-1) = 0`.
  One substitution (`Polynomial.eval₂_reflect_mul_pow` at `x = -1`), no factorisation.
- **(b)** `P = (X+1) Q` ⟹ `Q` is self-reciprocal with `N-1`, by `Polynomial.reflect_mul` and
  cancelling the nonzerodivisor `X + 1`.
- **(c)** Induct, stepping the multiplicity by `Polynomial.rootMultiplicity_mul_X_sub_C_pow`.

**The paper's preliminary reduction to $\Pi(0) \neq 0$ (writing $\Pi = t^\alpha Q$) is therefore
removable**: step (b) goes through verbatim when the constant term vanishes, and the induction
never appeals to it. This is a genuine simplification of a printed proof of mine, found by
planning the formalisation rather than by doing new mathematics.

## Generality

Stated over any `[CommRing R] [NoZeroDivisors R] [Nontrivial R]` with a hypothesis `(2 : R) ≠ 0`
— which is used in exactly one place, step (a). The paper's `ℚ` statement and the `ℤ` statement
the application wants are both corollaries:

```lean
theorem rootMultiplicity_neg_one_emod_two_rat {N : ℕ} {P : ℚ[X]} (hP : P ≠ 0)
    (h : IsSelfReciprocal N P) : P.rootMultiplicity (-1) % 2 = N % 2
theorem rootMultiplicity_neg_one_emod_two_int {N : ℕ} {P : ℤ[X]} (hP : P ≠ 0)
    (h : IsSelfReciprocal N P) : P.rootMultiplicity (-1) % 2 = N % 2
```

## Declarations (12, all sorry-free)

Main:
- `TworowD4Kernel.IsSelfReciprocal` (structure)
- `TworowD4Kernel.IsSelfReciprocal.eval_neg_one_eq_zero_of_odd` — step (a)
- `TworowD4Kernel.IsSelfReciprocal.reflect_one_X_add_one`
- `TworowD4Kernel.IsSelfReciprocal.of_mul_X_add_one` — step (b)
- `TworowD4Kernel.IsSelfReciprocal.rootMultiplicity_neg_one_emod_two` — **the lemma**
- `TworowD4Kernel.rootMultiplicity_neg_one_emod_two_rat` — the paper's statement verbatim
- `TworowD4Kernel.rootMultiplicity_neg_one_emod_two_int`

Non-vacuity and sharpness:
- `TworowD4Kernel.isSelfReciprocal_one_X_add_one`, `TworowD4Kernel.rootMultiplicity_neg_one_X_add_one`
- `TworowD4Kernel.reflect_one_X_sq`, `TworowD4Kernel.rootMultiplicity_neg_one_X_sq`
- `TworowD4Kernel.not_rootMultiplicity_emod_two_of_natDegree_gt`

## Build evidence

| check | result |
|---|---|
| `lake build` | `Build completed successfully (2981 jobs)`, exit 0 |
| `lake test` | `Built TworowD4KernelTests (4.1s)`, exit 0 |
| `sorry` count | **0** (`grep -c sorry` on the module: 0; no `declaration uses 'sorry'` warning) |
| linter | clean on this module (no long-line, no header warning) |
| `decide` / `native_decide` | **not used** — nothing here is a finite check |

## `#print axioms`

Asserted mechanically, not by eye: the axiom set of each declaration was parsed and compared as
a **set** against `{propext, Classical.choice, Quot.sound}`, with the check failing on any
unparsed output.

```
records parsed : 12
deviant        : NONE
unparsed output: ''
ASSERTION: PASS — all equal {propext, Classical.choice, Quot.sound}
```

(The first run of this checker **failed** — it double-prefixed the namespace and got 12 unknown
constants. It reported FAIL rather than vacuously passing on zero records, which is the property
the checker exists to have.)

## Scope — what this does NOT verify

This formalises **Lemma 3.3 only**. It does *not* formalise:

- **Lemma 3.2** (self-reciprocality of $\Pi_{\nu\lambda}$ itself, i.e. $h_k h_{N-k} = h_{N-k} h_k$
  paired against $s_\nu$);
- **Theorem B(i)** (the power-sum order computation, $\min_\mu \#\{\text{odd parts}\}$);
- **Theorem B(ii)** as a whole, which combines Lemmas 3.2 and 3.3.

Those need the ring of symmetric functions and plethysm, which this repository does not have.
The registry node added below is a *premise* of `Q105-vertex-order-is-parity-of-m-plus-n`, which
stays at `proved`. **This note is not evidence that Theorem B is verified.**

## Pending

CI run **`34322729509`** (head `4ab399c`, branch `main`, `axiom-audit: true`, allowlist = the
standard three) was **in progress** when this session ended — it was 1m11s in, and a successful
run takes ~45 min, so its status is genuinely unknown, not green. Local `lake build` / `lake test` are green and the axiom sets were checked
locally through the root import. CI status to be read off the run log at next WAKE, from the run
itself and not from a summary.

## Naming debt (for Robin, not acted on)

The repo is called `tworow-d4-kernel` and now holds cross-rank ribbon commutators and a general
polynomial lemma. The name no longer describes the contents. Renaming would break every commit
hash cited in every PDF sent under PROTOCOL §2.3, so this is a question for Robin, not a rename
I should perform.
