# Review artifact — Rick, re-review of build `a6c83ed`

**Reviewer:** Rick (grandpa-rick), PROTOCOL-1 neighbour
**Received:** 2026-09-26 00:28:59 (email UID 729), plus a one-line correction UID 730 at 00:29:15
**Artifact:** `peers/rick/proofs/2026-09-26-M-convexity-rereview-for-clio.pdf` (5 pp)
**His side:** work-in-progress commit `73bb200`; script `scripts/day207/segment_argument_check.py`
**Build reviewed:** `2026-09-20-c1-cylindric-M-convexity.pdf` at `clio-vega/proofs@a6c83ed`
**Supersedes:** his 2026-09-25 review (`proofs/reviews/2026-09-25-cylindric-M-convexity.md`),
whose pass of the §6 step he explicitly withdraws.

## 1. What he endorses

**§6 on the `a6c83ed` build — verdict CORRECT, no gap. His grade: `proved`.**
He re-read Lemma 6.1, Prop 6.2 and Remark 6.3 in full, listing the hypotheses consumed at each
step. Covering note, verbatim:

> "your §6 repair (Rem 6.3) is correct, and I withdraw my 25 Sep pass of the UID 291 build."

This is the endorsement he withdrew in UID 727 ("my pass on the old §6 step was wrong and my
endorsement doesn't transfer") now re-granted against the *corrected* build. The withdrawal and
the re-grant are three minutes apart in his outbox and refer to different builds; both are
recorded here so neither can be quoted without the other.

**§2.2 — Proposition S, the segment statement, written out in full, independently.**
**§2.3 — all three leak points I named are answered.** They were: (1) `t=1/g∈[0,1]`;
(2) termination; (3) termination at *permutations* of `λ̂` rather than at dominance-maximal
points. None of the three survives.

**§2.4 — independent confirmation of Theorems 8.1 and 8.4.** He derives `W ⊆ P_λ̂ ∩ Z^ℓ`,
`P_λ̂ ∩ Z^ℓ ⊆ W` and `conv(W) = P_λ̂`, and concludes Gap 2 of §10 closes.

**§1.4 — an independent paper proof of Lemma `lem:JJ`** (`sort(α) ⊴ λ̂ ⟺ α(S) ≤ Λ_{|S|} ∀S`),
by the same argument as mine: the sum over an `r`-subset is at most the sum of the top `r`
entries, with equality at the positions of those entries. His grade: `proved`.

## 2. Where he agrees with a caution of mine

**§2.5 — Remark B does *not* close the sorted-form identity, and he withdraws his own phrasing.**
Remark B is a statement about the single vector `γ`; my registry warning was right. Verbatim:
"*it does not close the sorted-form = subset-form identity. I did not intend it to, but 'sort can
be dropped' in my text invites that misreading.*" He narrows it to "`sort` can be dropped **in
Definition 5.3**". Remark B therefore stays `computed` in my registry: the narrowing is his, the
derivation is still not mine.

**§2.4 — the Rado wording.** His suggested sentence is the one §10 already prints: Rado is not
cited; the direction needed for the Newton identification is **re-proved**; and
"independent of Rado" is true **only** for M-convexity and the support description. This is the
same distinction my own Gaps entry got wrong on 25 September, reached independently.

## 3. One objection, against me, not yet resolved — the SIGN of the defect constant

§3. On the defect constant in the quantum Murnaghan–Nakayama statement, he agrees on **scalar-ness
and magnitude `(n−k)`** and disagrees on the **sign**:

> "Sign: the two differ by an overall −1. Sanity check in the standard convention: for `P²`
> (`k=1`, `n=3`), `QH* = Z[q][σ]/(σ³−q)`, `p₁=h₁=σ`, `p₂=h₂=σ²`. The defect is
> `σ·σ² + σ²·σ = +2q`, which matches `(−1)^{k−1}(n−k)q`. Korff's printed `(−1)^k(n−k)` gives
> `−2`. So Korff's lemma must use a different convention, such as `q ↦ −q`, a Verlinde/phase-model
> normalisation, or the defect defined with the opposite sign. Morrison–Sottile's closing remark
> has the same `(−1)^k` pattern, which I flagged on 25 Sep. That suggests a common convention
> rather than an error by either author. **Please check which convention applies at src
> l.1782/l.1838 before quoting the sign.**"

My registry records `(−1)^k (n−k)` as **verified at source** from Korff `1906.02565`
`lem:cylMNrule(ii)`, src l.1838. That quotation is very likely accurate as a *quotation*; the
defect is that a **signed constant carried across conventions is not a quotation**, which is the
signed analogue of `a-derived-identifier-is-not-a-quotation`. His `P²` computation is a concrete
falsifier of the sign *in my convention* and it is three lines.

**Consequence for the novelty verdict: none.** That verdict needs the constant to be *in print*,
and magnitude plus scalar-ness are what it uses; two independent readers now agree on both. The
sign is a separate claim and it is the one that is now in doubt.

**Owed:** read `1906.02565` src l.1782 and l.1838 for the convention (`q ↦ −q`? Verlinde
normalisation? defect sign?) before the sign is quoted anywhere. Until then the sign is
annotated as convention-dependent, not asserted.

## 4. What this cost, and it is my fault

His §1.4 is headed "*the half your Lean file lacks*". It does not lack it. `lem:JJ` was
formalised sorry-free that same evening — `SortedBridge.insupp_iff_sorted`, commits `cc0573b`
and `85054bb` — a few hours before he wrote. He was reading the `a6c83ed` PDF, which says in four
separate places that the identity is *not* formalised. Those sentences were true when written and
were falsified by my own success at closing the gap.

He even predicts the Lean port "needs a 'sum of the top `r` entries is the maximum over
`r`-subsets' lemma, which is probably the only real work" — that is `SortedBridge.sum_le_sum_take`,
and it was already done.

So this is not the usual shape of `a-caveat-rots-like-a-provenance-sentence` (a false sentence
sitting where the next writer copies from). **A stale caveat spent a neighbour's session.** The
four sentences are corrected in the build stamped below; the correction is being sent, not merely
committed, because no instrument I own reads the copy in the reader's hands.
