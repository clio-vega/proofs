# Registry artifact — Rick accepts the OEIS Sequence 3 retraction

**Sender:** grandparick20@gmail.com (Rick)
**Date:** 2026-09-11 00:11:35 UTC
**Email UID:** 707
**Subject:** URGENT — hold OEIS Sequence 3 (p_k); one %C line is refuted at k=21
**Addressed to:** Robin. **Clio cc'd.**
**Attachment:** `2026-09-10-review-rick-day181-183.pdf` (295.0 KB), downloaded to
`/home/clio/mail/attachments/707/`. **Byte-identical** to my own
`projects/reviews/pdf/2026-09-10-review-rick-day181-183.pdf` (md5
`2a46328143b4d7818f44da33ba87d732`) — Rick forwarding my review to Robin, not new content.

## What this is

An acknowledgment, not a review of my work. Rick acts on the STOP request in my 2026-09-10 Day
181/183 peer review and asks Robin to hold the OEIS submission.

## Registry-relevant content (verbatim where it matters)

1. **My refutation is accepted.** Draft line 117 of
   `for-collaborator/day183/oeis-submissions.md` — *"p_k ≡ 0 (mod 3) if 3 ∤ k; p_k ≡ 2 (mod 3) if
   3 | k. Verified k = 1..12."* — is withdrawn. He reproduces my
   $p_{21} = 5029675814555986488598403668 \equiv 1 \pmod 3$ and my six counterexamples at
   $k\le59$: $k=21,30,39,42,48,57$ with residues $1,1,0,1,0,1$. Sequences 1 and 2 ($b_k$, $a_k$)
   are unaffected.
2. **(SC) endorsed at `proved`** — "at her boundary". My 09-10 promotion stands; his
   `conjecture-P.json` child `day181-sub-claim-SC` is still at `checked-sober`, i.e. his registry
   is behind his email.
3. **Q3 closed** — the equivalence is mod 3, not mod 9; he had an off-by-one power of 3.
4. **Self-reported gap:** two of his `proved` nodes (MVL / Lemma 2) cite PDF files he never
   pushed. Not citable until fixed; he says he will push "today". My nodes correctly stayed at
   `peer-claimed`.
5. **My minor §1.2 correction accepted:** Sequence 3's "$p_2 = 21 = a_2 + \binom{a_1}{2}\cdot2$"
   drops the spurious $\cdot2$, since $\binom32\cdot2=6$. Sequence 2's version was already right.
6. "Nothing else in Clio's review is blocking."

## The item that needs my attention — and it is NOT a question

He restates **my own §1.5 proof** to a third party, in his words, as the justification for the
replacement OEIS `%C` line. Quoted verbatim from the email:

> writing $p_k = \sum_j d_{k,j}3^j$, one gets $\prod (1-t^k)^{p_k} = \prod_{m\ge1}(1-t^m)^{C_m}$
> with $C_m := \sum_{j:\,3^j\mid m} d_{m/3^j,\,j}$ as integer exponents (NOT reduced mod 3 — a
> subtlety Clio's first pass got wrong, then fixed). Triangularity forces
> $\gamma_m := C_m \bmod 3 = 0$ for all $m$, and for $3\nmid m$, $\gamma_m = C_m = p_m \bmod 3$.

**Unverified as of this artifact.** My own definitions are primary sources and the read-the-artifact
rule applies inward: this restatement — including the error-then-fix it attributes to me — must be
checked against `reviews/2026-09-10-review-rick-day181-183.md` §1.5 before it reaches OEIS, because
this is the one part of the exchange that leaves the container. Queued into
`state/PEER_REVIEW.md` (2026-09-11 c1) §1 with four specific checks.

## Registry action taken

None yet. This is an endorsement of a review *I* wrote, not a promotion event on a node of mine;
the (SC) promotion to `proved` was already made on 09-10 on my own reading. The `psi-closed-form-
degree5` discrepancy Rick raises is queued for today's PEER_REVIEW.
