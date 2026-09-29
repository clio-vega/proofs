# Review artifact — Rick, email UID 723, 2026-09-25 10:37

**Sender:** Rick (grandparick20@gmail.com)
**Subject:** reply re Theorem B
**Attachment:** `peers/rick/proofs/2026-09-25-reply-clio-theoremB.pdf` (5 pp, 258.5 KB)
**Registered:** 2026-09-29 WAKE (three days late — container down 09-26 08:00 → 09-29 12:10,
API weekly quota; see that day's digest to Robin).

## The objection, as he states it

> Theorem B holds for all `n` — but **do not claim it as new**. `R_e(−1)` is multiplication by
> `ψ(p_e)` in `QH*(Gr_{k,n})`, i.e. classical Newton pushed through Bertram's `ψ`. A quantum
> Murnaghan–Nakayama rule already exists: **Morrison–Sottile, arXiv:1507.06569, Thm 2**. Clio's
> `R_5(−1)` reproduces their example (4.5). The `r = n` defect is the scalar `(n−k)(−1)^{k−1}q`.
> Verified independently for `n ≤ 9`.

## Classification: this is a NOVELTY event, not a correctness event

**No trust downgrade is warranted and none is applied.** He affirms the mathematics ("holds for all
`n`", independently verified `n ≤ 9`). What he attacks is an attribution, and per the Q254
precedent a registry node is **not demoted on novelty**. `thmB-t-equals-minus-one` and
`propMN-consequence` stay `proved`.

## What I checked before recording this — the overclaim is NOT in the papers

Grepped every `.tex` in `proofs/` for novelty language around Theorem B and the quantum MN rule:

- `2026-09-23-c2-ribbon-height-is-a-heap-orientation.tex:450` states `prop:MN` as a **consequence**,
  explicitly *"Assume Theorem~\ref{thm:B}(ii)"*, and gap (G1) says it *"derives the quantum
  Murnaghan–Nakayama rule **from** (ii); it does not prove (ii)."* No novelty is claimed.
- `2026-09-24-c1-cylindric-newton-identity.tex:97` already says **"Theorem~\ref{thm:newton} is their
  corollary and belongs to Morrison–Sottile"**, and l.635/666 record *"I have **not** opened
  Morrison–Sottile"* / *"**Not read by me**"*.

So my own papers had already conceded the attribution, in writing, the day before Rick raised it.
**Rick's email confirms a concession I had already made rather than exposing one I had not.** That
is worth stating plainly instead of filing this as a near-miss.

## What is genuinely outstanding

1. **Morrison–Sottile `1507.06569` is `agent-summary` and unread by me.** Both my own paper and
   Rick's email are second-hand about what its Thm 2 says. Rick's claim that my `R_5(−1)`
   reproduces their Example 4.5 is *his* computation, not mine.
   **Owed: read `1507.06569` at source and check Thm 2 and Ex. 4.5 myself.**
   Until then the prior-art attribution rests on two `agent-summary` reports that agree — which is
   better than one, and is still not a reading.
2. **Theorem 5 (ribbon height = heap orientation) has had no novelty check at all.** Rick names this
   explicitly. It is the one genuinely structural claim in that paper and nobody has searched
   against it.
3. **The sign of the defect is disputed** (his UID 729): he and I agree on the scalar and on the
   magnitude `(n−k)`, and **disagree on the sign**. He asks that Korff's convention be checked at
   the cited lines before either of us quotes a sign. Note my own 09-23 paper carries both
   `(−1)^{k−1}(n−k)q` (Thm B(iii)) and a full-wrap value `(−1)^{k−1}kq` — those are different
   scalars and the relation between them should be stated, not assumed.

## Action taken

- Artifact filed (this document).
- Prior-art note added to `propMN-consequence` in
  `proofs/registry/ribbon-height-heap-orientation.json`; trust unchanged.
- Items 1–3 carried into `state/PEER_REVIEW.md` for the review slot.
