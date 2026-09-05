# `X0-closed-form-E3-zero` — convention question closed; textual correction, no demotion

**Date:** 2026-09-05
**Reviewer:** Clio Vega
**Registry:** `proofs/registry/rick-beta-prime-peer-claims.json`, node `X0-closed-form-E3-zero`
**Source reviewed:** `grandpa-rick/work-in-progress` @ `4fa7f30`; Day 158 PDF (proof source `6d48722`,
announced `1d6d480`); email UID 684 (2026-09-02), UID 691 (2026-09-04)
**Scripts:** `reviews/day158_Ebasis_diagnosis.py`, `reviews/day158_P1_chain_identity.py` (re-run today)

## What moved

**Grade: unchanged, `proved`.** This is *not* a demotion event.

**Text: corrected.** The `approach` field asserted, as of this morning:

> "CONVENTION SWAP, FLAGGED BY HIM AND NOT YET CHECKED HERE"

That sentence was false when read today. The check was performed 2026-09-03 and
sent to Rick the same morning; he conceded 2026-09-04 (UID 691). A
re-affirmation note was appended to the same field on 09-04 without correcting
the sentence above it, so the node simultaneously claimed the question was open
and recorded that it was closed.

## The question, and the answer

Rick asked (UID 684) whether his "top"/"sub-top" naming swapped against the
Day 152/154 tower, and whether the chain identity $\log W = \partial\Xi$ is
therefore false (he had checked $n=2$).

Re-verified today, independently, from the Day 152 §1 definitions:

| check | result |
|---|---|
| $\log W = \partial\Xi$, $E$-basis, $E_3$ free | holds, $n=1..6$ |
| restricted to $E_3=0$, $\partial$ applied **before** restricting | holds, $n=1..6$ |
| $\partial$ evaluated **inside** the $E_3=0$ slice | fails, every $n\ge2$ |
| dropped term $E_2\partial_{E_3}\Xi|_{E_3=0}$ accounts for the gap exactly | yes, $n=1..6$ |

$\partial = 3\partial_{E_1} + 2E_1\partial_{E_2} + E_2\partial_{E_3}$ has a third
term that does not vanish at $E_3=0$, so $\partial$ does not commute with
restriction to that slice. His $n=2$ computation was correct; the inference from
it was not. There is no convention swap.

## Why this was never a demotion risk

The brief for the session asked whether my `computed → proved` upgrade (09-04)
had quoted the identity under a superseded convention. It had not: the grade
rests on Day 158 **Theorem 1** ($[T^n]\Xi = E_2Y_n/n$) and **Theorem 2**
($X^{(0)} = \tfrac12\log W$), both of which I reproduced to $n\le10$ on an
independent pipeline, plus Prop. A verified symbolically. $\log W = \partial\Xi$
is a Day-152 identity upstream of the node, not part of the graded statement.

## Conditions carried forward (unchanged from 09-04)

The `proved` grade remains conditional on the three-line induction for
$\deg_u g_m \le m+2$ (Day 158 §3, asserted as "empirically (and provable ...)")
being adopted into the Day 158 text; the induction is supplied in
`reviews/2026-09-04-review-rick.md` §2.3. Scope is the statement only: it does
not extend to the deleted Day-158 consequence, and does not endorse C.5 or
Day 162 Theorem B.

## Memory

Third firing of `refutations-do-not-propagate-backwards`: a later session closed
the question and no back-edge updated the field that had asserted it open. The
09-04 pass appended rather than corrected. Annotate-don't-rewrite is right for
*conclusions*; it is wrong for a sentence that states an open question which is
no longer open.
