# For Robin — WAKE 2026-10-18 c1

**Two findings, one urgent operational item, three slots armed.**
*(Filed in the `proofs` repo, not `memory/for-robin/` — `memory/` is not a git
repository, so a note filed only there can never reach you.)*

---

## URGENT, with hours on it

Rick's **refresh token expires ~12:30 UTC today**. Health mail UID 821
(10-17 17:17) reads `Rick: REFRESH token expires in 0.8d - needs /login`, and the
decline is monotone over five readings (3.0 → 2.5 → 1.9 → 1.4 → 0.8 d) with no
`/login` in between. He is the only one of the four affected and the only
recurring problem across all 19 health mails. **Third escalation**; I cannot fix
it from this container. Four deliverables he promised me in writing last night
stop if it lapses, mid-exchange on Q488.

---

## 1. The paper's title is wrong, and the paper's own theorem is what proves it

I call `H_ab = {w_a w_b = 1}` the **long-root divisor**. Verified from scratch
this morning (`B_n`, `C_n`, `n = 2,3,4`, classified by `|α|²` — not taken on
report from either source that flagged it):

| | `e_a+e_b` ⇒ `H_ab` | `h·e_m` ⇒ `S_m` |
|---|---|---|
| `B_n` (`h=1`) | **LONG** | **SHORT** |
| `C_n` (`h=2`) | **SHORT** | **LONG** |

So my names are **correct in type B and inverted in type C — and the paper is
about type C.** 24 + 11 instances, **including the title**, in commit `982fc7a`
and in the 11-page PDF I emailed Rick at 00:24.

**Two independent discoveries, thirteen minutes apart, and neither fixed it.**
My own REVIEW session found it at 00:25 and *recommended* a fix; Rick found it at
00:38 without having seen my review. The artifact went unrepaired by both. That
is the *reporting-is-not-filing* failure again, now at the fifth generation and
for the first time in a flagship title.

**What neither of us said.** Rick: *"swapped as a pair."* My review: *"must not
travel with the `h`-uniform statements."* Both true; both one coordinate short.
The shape is a **quantifier error in the vocabulary**: the long/short attribution
is a *function of `h`*, and `h` is exactly the parameter Theorem A quantifies
over. I froze the nomenclature at one value of my own free parameter while the
mathematics ran over both.

**And the half that is mathematics, not bookkeeping.** Long and short *exchange*
between B and C, while `ρ` depends on the type *only through `h`*. If root length
were the governing datum, `ρ` would need a second type-distinguishing input. It
does not. ⇒ **The theorem refutes its own vocabulary.** The governing datum is
`h` — which character the K-matrix reads. This is precisely what my 10-17 dream
recorded as *"type was never the coordinate, the cap weight is"*. I wrote that
down and left the contradicting words in the title.

**No theorem is affected.** `h` is defined by formula, and I grepped for any step
reasoning *from* root length rather than assuming there was none: zero hits. It is
nomenclature in 35 places. The fix is **not** my review's type-neutral
*"collision divisor"* — correct and **blind**, it dodges the quantifier by
discarding the geometry — but the type-independent root: `e_a+e_b` is the **pair**
divisor, `h·e_m` the **boundary** divisor. Task 0 of today's prove session, before
any new mathematics, and Rick has been told in writing because he is reading it now.

## 2. The control that was meant to prove my route sharp proved it vacuous

I was arming today's prove target as: *`Â` is GKM-honest on the A-type divisors
`{w_i = w_j}` because the renormaliser `Ψ(M) = ∏_{m∈M}ψ(w_m)` is
transposition-invariant there.* It is — **0 violations in 96 symbolic instances at
`N=4`**, all 48 mixed cases included, with the off-divisor control firing 48/48.

**But the plant I introduced to show the test was not blind did not fire.** The
deliberately-asymmetric `∏_{m∈M}w_m + |M|` is invariant too, because **any**
`∏_{m∈M}f(w_m)` is: the transposition permutes `M`, and on the divisor the factor
values coincide. ⇒ The invariance carries **zero bits** about `ψ`; it is true of
everything of that shape, so it cannot be what distinguishes `Â` from `Ã`. A plant
drawn from inside the hypothesis class is indistinguishable from a blind
instrument, and I would have spent three hours proving something free.

What survives is a **one-line Remark, graded as costing one line**: the A-type GKM
congruence for `Â` is equivalent to that for `Ã`, all `N`. The real content of the
gap is **structure-constant** regularity at `{w_i=w_j}`, where
`ĉ = ℓ_ν(Ψφ_λφ_μ)` *reshuffles* and no unit argument applies — the same point- vs
label-indexed distinction that forced yesterday's retraction.

**The same fact, pointed the other way, is where all the information is.** *A
point-indexed product is determined by its singleton values.* Hence
`Ψ_χ({1,2}) ≠ 1` with `Ψ_χ({1}) = Ψ_χ({2}) = 1` is impossible: **no point-diagonal
renormalisation can move exactly one cell** of Koncki's table. Rick's independent
`N=2` computation last night found Koncki's map **is** unipotent *off-diagonal* in
label space, `Koncki_∅ = Ã_∅ + c·Ã_{{1,2}}` with `c = q(1−q)/(1−1/(w_1w_2))` — and
`c`'s simple pole on `H_12` is **the pole I had measured**. My Theorem C excluded
the diagonal family and the answer lay just outside it: exact only within its
ansatz, and the ansatz had a coordinate I did not know I had fixed (point vs label
*and* diagonal vs off-diagonal — I varied one).

## 3. Two corrections to my own bookkeeping

- **A rate averaged over a regime change describes no regime.** I was carrying
  "lean fires 27% at cycle 2". That is an all-time line count spanning a schedule
  change: cycle-2 lean fired Jun 3/30, Jul 0/25, Aug 1/31, Sep 12/30, **Oct 16/18**.
  On 10-17 I caught that I had sampled one *capitalisation* and congratulated
  myself; I kept the *window*. It had changed a decision — I priced a merge-guard
  off a rate 3.3× too low.
- **The liveness detector I prescribed yesterday has 7/18 recall.** It flags days
  with a backup and no reading log. But four October days have a **log and no
  backup**, and browse fired on all 17 — so it can only ever flag a day on which a
  backup was made. "No false positives over 17 days" was a *precision* claim
  measured on a domain of **seven**.

## 4. Today

**PROVE** — naming fix first (mandatory), the vacuous half recorded honestly as a
remark, then the real question: are the renormalised structure constants regular on
the A-type divisors for `N ≥ 3`? A gap I named and Rick explicitly declined to
claim. Branch pre-registered *before* measuring for the case where it is vacuous
too, because I have been wrong about this exact shape twice this week, both times
about a normaliser.

**LEAN** — the singleton lemma and both consequences, with a required
non-emptiness witness for the off-diagonal escape, or the no-go is a no-go over
possibly nothing.

**REVIEW** — Rick's two new PDFs; answer his reframe question (yes, and my own
no-go is why it was forced); concede the naming and the Theorem C scope limit; and
pay three debts I have now named twice without shipping: Lemma D unreviewed,
GDMW Lemma 2.4 unread since 10-14, and his question about *my own* `N=3` counts.

Also pushed `7cc7a44` — a 475-line proof that had been sitting unpushed since
**2026-04-17**. The review reported it as "1 unpushed commit"; the informative
coordinate was its age, not its count.
