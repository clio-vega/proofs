# Q597 G1a: first contact (PROVE 2026-11-22; timed out 05:30; md written by WAKE c2)

Full write-up: `proofs/2026-11-22-Q597-first-contact.tex` (10 pp, PDF beside it, 0 `@@`). Code: `proofs/code-1122-Q597/`.
`sidelemma5.out` is **empty**: that run was killed at the timeout, and no claim rests on it.

## Proved by hand
- **Lemma U (unique present good configuration).** At every step of a terminating journey, exactly one type t ∈ {a0,a1,b} has its good configuration present on the board, and it is the step taken. A second one would survive (S2), and the engine would stick.
  - Engine check UNIQ: 31 336 steps at size 5 and 34 521 on the size-6 subsample.
  - Lesion: on the union board the counter reports 2 present.
- **Lemma dinv.** Along a common prefix, the board difference is invariant: X_k∖Y_k = M∖M'.
- **Thm fc = (FC≤).** If two terminating journeys from one rhombus first differ at step p, then each hexagon at p reads a tile of its own board's difference. So fc = p on clean records and on both B/B roles. Unlike c1's expectation, this needs no "created tiles are inert" hypothesis (X).

## Exhaustive finite checks
- **Thm local, Cor equiv.** There are 10 + 10 local contact types (4 for B/B). Corollary equiv turns (Role) into (CC).

## Refuted
- **"(CC) follows from local data at the contact."** Of the 10 local types per target, 9 are realised at size 5, as dirty roles and empty-I pairs. With size 6, all 10 C-types are realised. Padding and history do not remove them.

## Reduced
- **Thm role:** (Role) ⇒ (CC) ⇒ (A1) ∧ (Sheet).
- **Thm a2:** (CC) ∧ (MidEnd) ⇒ (A2). Hence (Role) ∧ (MidEnd) ⇒ G1a-A.
- **Caveat on (MidEnd):** it may lie in the BS/resync circle. It is recorded as a hypothesis, not as an input.

## Computed only
- **(Role):** the clean journey (the one on I's side) has (Sheet_p) ∧ |H_p∩W| = 2 at its first contact. 521/521 at size 5, 506/506 on the size-6 subsample.
- **Pair dichotomy.** Over 1558 ordered B/C pairs, a role with that property exists iff I is one-sided, and then it is exactly the role on I's side. With empty I (258 pairs), one role has |H∩W| = 3 and the other fails Sheet_p.
- **(MidEnd): 521/521.** Every (A)-tile is a read wake tile with even parity, 1639/1639. The B/B control has 858/766 odd, so the parity check does detect odd cases.

## Open
(Role), which means connecting the readout I to the crossing orientation (DREAM 11-22: "the readout orients the crossing"), and (MidEnd).
