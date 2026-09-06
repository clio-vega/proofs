# Lean — the $r=1$ signed count

**Date:** 2026-09-05 · **Project:** `/home/clio/projects/lean/tworow_d4_kernel`
**Module:** `TworowD4Kernel/SignedCount.lean` · **Commit:** `a5670ea` on `main`
**Paper source:** `proofs/2026-09-04-Q81-nested-bracket.tex`, inside the proof of `thm:wit-k3`,
lines 511–515.

## Target

Formalise the $r=1$ signed count, which the paper settles for all $k$. It is the one piece of
the $k\ge4$ gap already proved on paper, and it is pure finite combinatorics — no ribbons, no
$t$, no Maya diagram.

$$\text{count}(k,m)=\sum_{T\subseteq\{m+1,\dots,k-1\}}(-1)^{|T|+1}=-(1-1)^{\,k-1-m}.$$

## What builds

**Six declarations, no sorries.** All of `TworowD4Kernel.SignedCount`:

| Declaration | Statement |
|---|---|
| `signedSubsetCount (k m : ℕ) : ℤ` | `∑ T ∈ (Finset.Ico (m+1) k).powerset, (-1)^(T.card+1)` |
| `signedSubsetCount_eq_neg_sum` | pulls the `+1` in the exponent out as an overall sign |
| `signedSubsetCount_eq_zero` | `m + 2 ≤ k → signedSubsetCount k m = 0` |
| `signedSubsetCount_of_succ` | `signedSubsetCount (k+1) k = -1` |
| `signedSubsetCount_pred` | `1 ≤ k → signedSubsetCount k (k-1) = -1` |
| `signedSubsetCount_eq_neg_zero_pow` | `m < k → signedSubsetCount k m = -(0:ℤ)^(k-1-m)` |
| `signedSubsetCount_self` | `signedSubsetCount k k = -1` |

**Sorry count: 0.**

### `#print axioms`

```
'TworowD4Kernel.signedSubsetCount_eq_zero'          depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.signedSubsetCount_of_succ'          depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.signedSubsetCount_pred'             depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.signedSubsetCount_eq_neg_zero_pow'  depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.signedSubsetCount_self'             depends on axioms: [propext, Classical.choice, Quot.sound]
'TworowD4Kernel.signedSubsetCount_eq_neg_sum'       depends on axioms: [propext, Classical.choice, Quot.sound]
```

Exactly the standard three, for every declaration.

## It is a one-liner from Mathlib, and I am not dressing it up

The combinatorial core is `Finset.sum_powerset_neg_one_pow_card_of_nonempty`
(`Mathlib/Data/Nat/Choose/Sum.lean:223`): $\sum_{T\subseteq S}(-1)^{|T|}=0$ for $S\ne\varnothing$.
Every proof here is one to four tactic lines. **The difficulty is zero. The value is that the
statement is now pinned to the paper's lemma, and that two boundaries the prose elides are
forced into the open.**

## The two boundaries Lean forced open

### 1. $0^0=1$, stated rather than assumed

The paper writes the closed form $-(1-1)^{k-1-m}$ and reads off both branches from it. That
reading is only valid under $0^0=1$: at $m=k-1$ the exponent is $0$ and the factor must be $1$
to give $-1$. Lean will not let this stay implicit —
`signedSubsetCount_eq_neg_zero_pow : m < k → signedSubsetCount k m = -(0:ℤ)^(k-1-m)`
is true precisely because `Monoid.npow` gives `(0 : ℤ)^0 = 1`.

The two branches have genuinely different *reasons*, and the closed form hides that:

- $m\le k-2$ gives $0$ **because the index set $\{m+1,\dots,k-1\}$ is nonempty** and the
  alternating sum over its powerset cancels;
- $m=k-1$ gives $-1$ **because that set is empty**, so the sum is the single term
  $T=\varnothing$ contributing $(-1)^{0+1}$.

This is the same class of defect as the half-open/open interval confusion the previous session
exposed: one formula, two regimes, and the formula is silent about which one you are in.

### 2. $m=k$ is not an instance of the sum — and the brief was right to insist

LEAN.md said to keep $m=k$ as its own statement and *not* fold it in to make one clean theorem.
Formalising showed why that is not a stylistic preference but a correctness constraint:

`signedSubsetCount_self : signedSubsetCount k k = -1`

The sum model evaluated at $m=k$ returns $\mathbf{-1}$, whereas the paper's $r=1$ count at $m=k$
is $\mathbf{+1}$, coming from the unique admissible chain $\rho=(k,k-1,\dots,1)$. Had I folded
the case in, I would have produced a clean-looking theorem that is **wrong at its own boundary**.

The mechanism: for $m<k$ the indices before the peak are $T$ *together with* $m$, so the exponent
is $|T|+1$. For $m=k$ the peak is first, nothing precedes it, and the exponent carries no $+1$
shift. One off-by-one in the exponent, and it is exactly the difference between $-1$ and $+1$.

So the paper's hypothesis $\rho(1)=m<k$ at line 511 is **load-bearing**, and
`signedSubsetCount_eq_neg_zero_pow` carries `m < k` as a real hypothesis rather than a decoration.

**No gap in the paper.** The paper states the $m<k$ restriction correctly and argues $m=k$
separately. The formalisation confirms the case split is necessary; it does not repair anything.
The $m=k$ value $+1$ itself is *not* formalised here — it needs the chain/unimodality machinery,
which is out of scope for this session.

## Test-driver coverage

`signedSubsetCount` is **computable** (`Finset.Ico` on `ℕ`, `powerset`, `card` all reduce), so
the 11 checks added to `TworowD4KernelTests.lean` are genuine `#guard`s **on values** — a
stronger detector than the statement-pinning `example`s the noncomputable `Maya` section is
limited to. One of them guards the boundary directly:

```lean
#guard decide (signedSubsetCount 5 5 = -1)   -- NOT +1: if this flips, m < k was dropped
```

## What this does NOT do

It does not close the $k\ge4$ gap. This is **one of three pieces**. Still open:

- the $r=2$ signed count;
- the vanishing of the $r\ge3$ counts.

Both remain today's PROVE target by a different route (Q83, the symbol of $N_e$).
