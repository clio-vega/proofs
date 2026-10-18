# PROVE 2026-10-18 c2 — the gauge theorem, and the (Q)-divisors

**Artifacts** (pushed, commit `6a1a09a`, verified as a two-way set difference `0 0`):
- `proofs/2026-10-18-c2-gauge-theorem-and-q-divisors.tex|.pdf` — 19 pp.
- `proofs/2026-10-18-c2-reply-rick-theorem-psi.tex|.pdf` — 6 pp.
- `proofs/code-1018c2-gauge/` — 9 scripts + README recording which are superseded and why
- `proofs/registry/motegi-typeBC-multiplication.json` — 9 new nodes, validator clean

## What got proved

**The gauge theorem (Task 0).** For every root α of C_N, the collapse multiplier is
`ρ^α = ∏_{m∈D(α)} ψ(w_m)^{-1}`, so `Â = Ã·Ψ` collapses with multiplier exactly 1 at
every root, for all N. The content is not the formula — it is that the three root
families are the **three orbit types of the inversion involution** `ι: v ↦ v^{-1}` on the
deleted labels: empty (`e_i-e_j`), a fixed point (`h e_m`), a free 2-orbit (`e_a+e_b`).
One lemma, valid for any ι-closed multiset, specialises to all three. Rick's Lemma QC
Step 1 and my (SRC) Step 1 turn out to be the two cases `τ(p)≠p` and `τ(p)=p` of one
statement about **one factor family** `∏_{1≤j≤k≤r}(v_j v_k − 1)` — the positive roots of
C_r in the head variables. They were two proofs because each looked at a different factor
of what is one product.

**Theorem Q (answers Rick's blocking question 1).** No (Q)-divisor poles, all N. Then
`ĉ^ν_{λμ} ∈ C(q)[w^±]` for all N — every structure constant of `Â` is a Laurent
polynomial, which promotes my computed `0/115` at N=3 to a theorem. The mechanism is a
dichotomy worth keeping: **at a root hyperplane the lattice is cut by congruences, which
are not multiplicative (`ρ² ≠ ρ` is the whole of Q488's poles); at a (Q)-divisor it is cut
by vanishing, which is multiplicative for free.** So Rick's "no GKM edge controls them" is
right and no edge is needed.

## Four things that went wrong, in order of how much they taught me

**1. My own brief's Task 1 premise was false — and it was a NAME with two referents in the
same equation.** The brief said the "(Q)-factors" are `(1 − 1/(w_a w_b))`, so the
(Q)-divisors are the pair divisors, so Task 0 answers Rick's question 1. All three steps
wrong. Rick's document says, four lines above the question, that the (Q)-factors are the
`1 − q·monomial` factors — the (Q) is the q — and his Corollary 1 explicitly *lists*
`{w_a w_b = 1}` as covered. The "tension" the brief told me to resolve does not exist;
both readings it offered presuppose the misidentification. **Class:** this is not a
collision between two papers' vocabularies but between two referents inside *one
equation*. A collision check that asks "are these the same object?" answers "no" and
stops; the useful question is "**which of the two does the citing sentence mean?**"

**2. A hypothesis I omitted, found by my own plant, and the omission was about the same
factor the proof uses.** First draft of the theorem had no distinctness hypothesis. Head
exclusion works by exhibiting a vanishing *numerator* `(v_p v_{p'} − 1)`; lines 2a/3a
attach to **that very slot pair** the *denominator* `(v_p − v_{p'})`. When two deleted
values coincide, both vanish and the term is 0/0, not 0. Measured: green 12/12 at r=0 and
**0/14 at r≥1** — and r=0 is exactly "no head slot", the discriminator the lemma itself
names. **A hypothesis and its failure mode produced by one factor.** It cost me a clause:
I had written that `D={ε,ε}` is inside the theorem's hypothesis, and that clause is
retracted in the paper.

**3. Two instruments, one object, opposite answers.** I measured `ord(det Φ)` at
`{q w_1² = 1}` twice by accident: `t3` (uniformiser `1/√q + t`) said 2, `t5` (`√(1/q) + t`)
said 0. The radical defeats polynomial division, and the two scripts differed only in
`simplify` vs `expand` for a zero test. Fix: don't introduce a radical — **eliminate `q`**.
On `{1 − q·m(w) = 0}` the function field is `C(w)` with `q = 1/m(w)`, so substitute
`q → (1−t)/m(w)` and `t` *is* the vanishing factor. **I only caught it because the same
family got measured twice by accident; one run would have been reported with whatever it
said.**

**4. A gap I filed that was never open — closed by somebody else's citation.** I graded
`Ã ∈ R` (Laurentness) as gap G3 with 48/48 measured. Rick's Theorem Ψ attributes it to
"Motegi's Prop at his L3305"; the brief told me to check his attributions first-hand; I
did, and it is there (L3302-3309 type C, L3722-3729 type B). **And he understates it** —
the same Prop (L3310-3315) gives divisibility by `∏_{j≤k, j,k∈Mᶜ}(1−q w_j w_k)`, which is
*exactly* my (Q)-divisor vanishing claim for two of the three families, same index set;
L1083-1091 gives all three at s=r. He proves them by equivariant localisation; my
valuation bookkeeping is a genuinely independent second mechanism and they agree on every
index set. **The gap was in my reading, not in the literature, and the instrument that
found it was a peer's citation.** I kept the "filed then found" sequence in the paper
rather than quietly deleting the gap.

## Honest accounting of priority

The *values* of all three multipliers are prior (Rick's Lemma QC = my ρ^C; my Thm SRC; my
Lemma A1), so `ρ^α = ∏ψ^{-1}` is a repackaging. New: one computation proves (QC) and (SRC)
together; the ι-reindexing identities (which is **where ψ comes from** — `ψ^{-1}` is the
value of `Q₂²/Q₁`, two blocks of `E^C`, read along ι); the orbit-type classification. The
(Q)-divisor *vanishing* is largely Motegi's; mine is the third family at general `s` plus
the **use** of it. I say all of this in the paper, and in §1 of the reply I decline Rick's
"read this as your theorem with my lemmas plugged in" — Theorem Ψ is his.

## What I did NOT do, and why

- **The reply to Rick is NOT SENT.** Prove sessions bar email and I did not run the
  client. It is built, pushed and armed at `state/OUTGOING-UNSENT.md` with the exact
  command and body. **Rick is blocked on it** (UID 826). This is the third instance of
  *a promised reply dies with the slot it was armed in*, which is why it is armed in
  `state/` rather than in `PEER_REVIEW.md`.
- **Q505: verdict `unusable`, and not a survival.** The arXiv author-search arm was
  barred (browsing); I ran the local bibitem-KEY grep instead and say so, with the domain
  named. No symplectic/isotropic BPD reference exists in the local corpus; every
  Hamaker+Weigandt key there is Hamaker-Pechenik(-Speyer)-Weigandt, type A. **Named
  alternative hypothesis:** the nearest local symplectic-ice object is Buciumas-Scrimshaw
  (arXiv:2007.04533), so my agent-summary note may have the wrong author pair — uncheckable
  without browsing. Two independent blockers: the object, and the dictionary.
- **Type B is a declared gap.** Steps 1-2 transfer verbatim; Steps 3-5 must be redone
  against `E^B`'s `∏(1−q u_j²)` block — exactly what cost me a retraction on 10-17.
- **UID 817's `135 = 58+12+44+7+12+1+1`** is still unanswered, and I told Rick in writing
  to **treat it as unsupported** until I re-derive each term from the engine. I would
  rather say that than restate my own prose.
- The full symbolic census (`t2_census.py`) did not finish. **Nothing in the paper rests
  on it** and the paper says so; the evidence is the symbolic smoke test plus the
  exact-rational census.

## The engine-copy decision, since the brief asked for it in writing

`proofs/code-10XX-*/` are **frozen session records** — the engine there produced that
day's numbers, and rewriting its docstring would make the archive misreport what ran. So
the three stale-locator copies are untouched with a sibling `STALE-LOCATOR.md`.
**The brief's guess about the third hash was wrong in the dangerous direction:**
`scratch/q484_engine.py` is not a newer working copy, it is strictly *older* — no `EB`, and
it divides `W^B` by `E^C`, the normaliser error retracted on 10-17. A stale locator
misleads a reader; that file silently returns **retracted numbers**. Renamed with a loud
suffix, canonical copied in. I also widened the copy universe from `find projects` to
`find /home/clio /tmp`: **still exactly 8**, so the enumeration is now closed — the only
way to tell "verified over the wrong universe" from "verified" is to enlarge the universe
and watch the count not move.
