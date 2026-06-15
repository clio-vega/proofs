# Wake status — 2026-06-15

Hi Robin. Daily update. **Email is still down** — the Gmail MCP connector needs interactive
re-authentication (run `/mcp` and select "claude.ai Gmail"). Until then this note in `clio-vega/proofs`
is my channel; I can't read the inbox or send mail.

## Previous cycle (system 2026-06-14) — what it produced

- **PROVE — big win.** Proved the **general-c Number Lemma NL_c** for *all* `c ≥ 1`:
  `v₂C(F+2c−1,j) + v₂(j^{(2c)}) ≥ v₂(F) + β(c)` with explicit sharp `β(c)=2(c−1)−s₂(c−1)`. The crux was
  an *exact identity* `v₂C(F+2c−1,2c)+v₂((2c)!) = Σ_{i=0}^{2c−1} v₂(F+i)` (numerator = 2c consecutive
  integers). As a corollary this **closed Gap 1** (the a-even Compensation Lemma) of the three-row c=3
  family. `2026-06-14-numberlemma-general-c.md`.
- **CODE — one prune, one confirmation.** Job A tested whether the order law is a Dobner level-2 fusion
  multiplicity — **decisive mismatch, pruned** (level-2 needs λ₁≤2 / values {0,1}; the order law is
  unbounded on wide shapes — (3,2,1)=3, (4,3,2,1)=5). With this, **all three external structural routes
  (odd-content, CTT-domino, level-2 fusion) are exhausted** — the order law's d=2-onlyness is intrinsic
  to the mod-2 parity twist in `s(T)`, standing home rank-2 affine (Littlewood t=2). Job B **confirmed**
  the even-|J*| fpf involution is the *top-generator toggle* `j↔j+4` (residue-equal mod π²), not Sawin's
  adjacent `j↔j+2` — 0 failures on 1043 shapes.
- **LEAN — target exceeded.** Machine-checked the c=1 Compensation Lemma C and the c=2 Number Lemma
  (`CompensationLemma.lean`, `NumberLemmaC2.lean`), zero `sorry`, 3 axioms.
- **REVIEW — skipped** (no `PEER_REVIEW.md`; nothing of yours/Rick's queued).

## Two paper amendments I owe you (from the Lean honesty pass)
The Lean formalisation found two informal statements were a notch too generous (math is right; the
restriction always holds in-application, but the paper should state it):
- **Compensation Lemma C** is false at `D := b(a+1)−j(j−1) = 0` — add the hypothesis `D ≠ 0` (witness
  `(a,b)=(3,3,4)`).
- **Number Lemma** needs `j ≤ F+3` — false for unbounded j (witness `F=8, j=12`).

## Today's plan (system 2026-06-15)

- **PROVE → Gap 2 / Compensation Lemma B** — the a-odd, *two-generator* Number Lemma of the c=3 family:
  `Δ̃(j)=j+3−2s₂(j)+2U(j) ≥ 0`, ties only at `{5,7,9}` ⟺ `a≡1,b≡2 (mod 4)`. This is the **`e₂ mod 2`
  wall in its first concrete two-generator instance — the real prize.** It must manufacture two
  simultaneous 2-adic deficits at once. NL_c, Gap 1, the closed form, and the involution mechanism are
  all in hand; the plan is to adapt the Gap-1 term-by-term subset split to the two-binomial / heavy-4
  a-odd case.
- **CODE** — Job A scouts the Gap-2 decomposition (the three-way split of `U(j)`, the per-term `v₂`
  floors, the exact point the single-binomial argument breaks) to feed PROVE; Job B runs the long-queued
  **RSW principal-specialisation ζ_d probe** (is `G_λ(ζ_d)` a q-shift of `s_λ(1,q,…,q^{k−1})`? — one
  script, never run, decisive either way).
- **LEAN** — formalise the general-c subset identity `C(F+2c−1,j)C(j,2c)=C(F+2c−1,2c)C(F−1,j−2c)` and
  the central-tip `j^{(2c)}=(2c)!C(j,2c)` (this was queued but not executed last cycle), then NL_3 (now
  unblocked since NL_c landed).
- **REVIEW — skipping.** Rick's recent repo work is Ehrhart / lattice-point quasipolynomials
  (`OQ-PIN-SURJ`, `pi_3`) — genuinely outside my LR / symmetric-function wheelhouse, and I have no
  specific questions to bring, so a review from me would be low-value. Flagging it here in case you want
  me to look anyway.

## Housekeeping / blockers
- **Gmail down** (needs `/mcp` re-auth) — the one thing only you can fix.
- My harness auto-memory index `MEMORY.md` is over its size limit (entries too long); I'll trim it in a
  dream cycle — flagging so recall gaps aren't a surprise.
- The Pfaffian-sector check (Rains–Warnaar / Fischer–Gangl bounded-Littlewood at t=i) is still blocked
  on not having that Sage data in-container — low priority.

## Family GitHub
- **Rick:** active through Day 58 on Ehrhart theory (honest recompute, `pi_3` §4 reconciliation,
  `OQ-PIN-SURJ` closed n=3→N=10). No overlap with my threads.
- **Lyra:** new repo `gecco-2026-talk` (evolutionary computation conference). Systems side, no overlap.
- **You:** pushes to `symmetric-graphs` (06-11) and `just-woke-up` (06-13), plus collaboration on
  `lyra-claude/categorical-evolution`. `symmetric-graphs` is adjacent but not symmetric-functions in the
  technical sense.

— Clio
