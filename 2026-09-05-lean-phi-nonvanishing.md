# Lean session 2026-09-05 (cycle 2) — the integral-domain core of Q83's main theorem

**Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**New file:** `TworowD4Kernel/PhiNonvanishing.lean`, imported from the root module
`TworowD4Kernel.lean` (line 13), so it is inside the import closure the CI `axiom-audit`
follows.

**Paper source:** `~/projects/proofs/2026-09-05-Q83-sharpness-all-k.tex`, `thm:main`
("Sharpness for all $k$"), the proof paragraph at lines 343–350.

---

## Status: both jobs land, sorry-free. Zero sorries in the file.

```
$ grep -c sorry TworowD4Kernel/PhiNonvanishing.lean
0
$ lake build
Build completed successfully (2977 jobs).
```

### `#print axioms`

```
'TworowD4Kernel.PhiNonvanishing.Phi_ne_zero'           depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.PhiNonvanishing.Phi_natDegree'         depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.PhiNonvanishing.natDegree_window'      depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.PhiNonvanishing.natDegree_prod_factors' depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.PhiNonvanishing.prod_factors_ne_zero'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.PhiNonvanishing.window_ne_zero'        depends on axioms: [propext, Classical.choice, Quot.sound]
```

Equality with the standard three, asserted as an allowlist, not as a denylist on
`Lean.ofReduceBool`. No `native_decide`, no `decide` on anything, no `sorry`.

---

## The definition

The paper writes $\Phi$ with a division. I did not formalise a quotient. For $a<b$,

$$\frac{x^{a}-x^{b}}{1-x} \;=\; x^a + x^{a+1} + \dots + x^{b-1},$$

and the geometric sum is what the paper's own derivation produces (the window
$j-e_{\max} < \Sigma_T \le j-e_{\min}$) *before* it is contracted into a fraction. So the
sum is the definition:

```lean
noncomputable def Phi (a b : ℕ) (e : List ℕ) : Polynomial ℤ :=
  (∑ i ∈ Finset.range (b - a), Polynomial.X ^ (a + i)) *
    (e.map fun m => 1 - Polynomial.X ^ m).prod
```

with $a=\min(e_{k-1},e_k)$, $b=\max(e_{k-1},e_k)$, $e=[e_1,\dots,e_{k-2}]$.

## JOB 1 — the target

```lean
theorem Phi_ne_zero {a b : ℕ} {e : List ℕ} (hab : a < b) (he : ∀ m ∈ e, 1 ≤ m) :
    Phi a b e ≠ 0
```

This is the paper's sentence *"each $1-x^{e_i}\ne0$ in the integral domain $\Z[x]$, so
$\Phi_{\vec e}\ne0$."* The two factors genuinely need **different** evaluation points, and
there is no single point at which both are visibly nonzero — that is the whole content of
the step, and it is what the prose hides:

| factor | seen nonzero at | value there | why the other point fails |
|---|---|---|---|
| $1-X^m$ | $x=0$ | $1$ | at $x=1$ it is $0$ |
| $\sum_{i<b-a}X^{a+i}$ | $x=1$ | $b-a$ | at $x=0$ it is $0$ unless $a=0$ |

`mul_ne_zero` over `IsDomain (Polynomial ℤ)` closes it.

### Finding: the paper's $e_i\ge2$ is not load-bearing here

The paper carries $e_1,\dots,e_k\ge2$ throughout `thm:main`. For **this** step only
$e_i\ge1$ is needed — $m\ge1$ is exactly what makes `eval 0 (1 - X^m) = 1`; $m\ge2$ buys
nothing. The Lean statement is therefore at `1 ≤ m`, and it is strictly stronger than the
paper's.

This is not pedantry. The $\ge2$ *is* used in the very next sentence (it is what gives
$e_{\min}\ge2$, hence the degree window), so it is a real hypothesis of the theorem —
but carrying it into the nonvanishing argument would have been a pinned parameter sitting
in a lemma that does not need it, and the harness would then have been unable to tell me
whether it was load-bearing.
→ `a-fixed-parameter-can-be-the-whole-obstruction`

## JOB 2 — the degree identity

```lean
theorem Phi_natDegree {a b : ℕ} {e : List ℕ} (hab : a < b) (he : ∀ m ∈ e, 1 ≤ m) :
    (Phi a b e).natDegree = (b - 1) + e.sum
```

The paper's next sentence: *"$\deg\Phi_{\vec e}=e_{\max}-1+\sum_{i\le k-2}e_i=E-e_{\min}-1$,
so $0\le j_0\le E-1$ and $(E-j_0,1^{j_0})$ is a genuine partition of $E$."* In the notation
of `Phi`, $e_{\max}=b$ and $\sum_{i\le k-2}e_i=$ `e.sum`, so the claim is exactly
`(b - 1) + e.sum`. **The off-by-one is where I expected trouble and it is clean:** the
window $X^a+\dots+X^{b-1}$ tops out at $b-1$, not $b$, and the identity is an equality, not
a bound — I did not need the weaker `≤` fallback the brief authorised.

Supporting lemmas, all sorry-free:

* `natDegree_factor` — $\deg(1-X^m)=m$ for $m\ge1$, via
  `natDegree_sub_eq_right_of_natDegree_lt`.
* `natDegree_prod_factors` — $\deg\prod(1-X^{e_i})=\sum e_i$, by list induction on
  `natDegree_mul`, which needs both factors nonzero and so **consumes JOB 1's
  `prod_factors_ne_zero` as a hypothesis**. JOB 2 is not independent of JOB 1; it rests on it.
* `natDegree_window` — $\deg\sum_{i<b-a}X^{a+i}=b-1$, by
  `natDegree_eq_of_le_of_coeff_ne_zero`: bound above by `natDegree_sum_le_of_forall_le`,
  and the coefficient at $b-1$ is $1$ because exactly the term $i=b-1-a$ survives.

## Cross-check against the paper's own worked example

`2026-09-05-Q83-sharpness-all-k.tex`, the example *"the excluded case of `thm:wit-cited` is
not an obstruction"*, takes $(e_1,e_2,e_3)=(4,2,3)$ and computes by hand
$\Phi_{\vec e}(x)=\frac{x^2-x^3}{1-x}(1-x^4)=x^2-x^6$. In `Phi`'s notation that is
$a=2$, $b=3$, $e=[4]$. Three `example`s in the file check it:

```lean
example : Phi 2 3 [4] = Polynomial.X ^ 2 - Polynomial.X ^ 6
example : (Phi 2 3 [4]).natDegree = 6
example : Phi 2 3 [4] ≠ 0
```

Neither number was tuned. The degree $6$ comes out of `Phi_natDegree` as $(3-1)+4$ and
independently equals the paper's $E-e_{\min}-1=9-2-1$ — two routes to the same number,
one from the Lean theorem and one from the paper's formula, agreeing.

## Coverage note — read this before inferring test-driver coverage

`Polynomial ℤ` is `noncomputable`, so **nothing in this file can be exercised by `#guard` in
`TworowD4KernelTests`**. The three `example`s above are elaboration-time checks inside the
library, not test-driver checks. Validation here is `lake build` plus the `#print axioms`
allowlist above, and nothing else. In particular the standing gap in
`CI-NEGATIVE-CONTROLS.md` still applies: a failing `test` step leaves `axiom-audit`
reporting `skipped`, and "axiom-audit skipped" can never be read as "axioms fine".
→ `inherited-validation-covers-only-executed-paths`

## What is *not* formalised

`Phi_ne_zero` and `Phi_natDegree` are the two sentences of the paragraph. The paragraph's
remaining content — that $[x^{j}]\Phi_{\vec e}$ *is* the matrix entry (eq. `phi`), and that
the $(1+t)$-adic valuation conclusion follows — is not touched, and the surrounding
`thm:main` remains `proved`, not `lean-verified`. The Lean nodes below are attached as
children, not as a re-grading of their parent.
