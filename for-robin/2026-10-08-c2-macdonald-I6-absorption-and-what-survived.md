# M^(L) is Macdonald I.6 — and the one step his Example 10 does not take

**2026-10-08, cycle 2 (PROVE).** Paper: `2026-10-08-c2-exact-denominators.tex` (7pp, compiles clean).
Code: `2026-10-08-c2-exact-denominators-code/`.
URL: https://github.com/clio-vega/proofs/blob/main/2026-10-08-c2-exact-denominators.tex

## The short version

Yesterday I wrote a paper identifying `M_{ρν} = ⟨p_ρ,h_ν⟩`. Today's brief sent me to prove a
corollary of it behind a novelty gate. Three things happened, in increasing order of how much they
should change what I do next.

1. **The corollary was already proved** — in yesterday's paper, three lines below the theorem the
   brief named as "the likely route" to it.
2. **The gate answers *equal*, and not by computation.** `⟨h_ν,m_σ⟩=δ` gives
   `⟨p_ρ,h_ν⟩ = [m_ν]p_ρ` outright. So `M` **is** Macdonald's `L = M(p,m)`, I.6 (6.8). That was
   already Lemma 1 of yesterday's paper. The gate was answered inside the document the gate was
   written to test.
3. **Essentially the whole cluster is Macdonald I.6**, including the strengthening I had planned as
   today's actual contribution.

## The absorption, node by node

Verified first-hand from **rendered page images** at 200 dpi (PDF pp. 114, 121 = book pp. 103, 110),
not from the text layer — my own `sources.json` entry for this scan records that its OCR "yields
garbage that looks like content", and I had violated that instruction on the first pass.

| my node | Macdonald I.6 |
|---|---|
| `lem-three-faces` | (6.9) |
| `thm-structure-diagonal-blocks` — billed "THE RESULT" yesterday | (6.9)+(6.10)+**Remark p.103** |
| diagonal entry `d_ρ` | Ex. 10, parenthetical — *his* proof is the orbit–stabiliser argument I gave |
| `cor-det-M-L` | two lines from Ex. 10 |
| `cor-exact-denominators` | two lines from Ex. 10 |
| **today's intended strengthening** | **Ex. 10, final sentence, verbatim** |

That last one, quoted exactly: *"Since LU⁻¹ is unitriangular, (m̃_μ) is a **Z**-basis of the subring
**Z**[p₁,p₂,…] of Λ."* I had planned to prove that the lattice spanned by the power sums equals the
lattice spanned by the augmented monomials, hence `Λ_n^Z / P_n ≅ ⊕_ρ Z/d_ρ`. I verified it
computationally (SNF, all 36 pairs `(n,L)`, `n≤8`, 20 with non-cyclic cokernel, 0 mismatches)
*before* opening the page. The verification stands; the novelty does not.

**This is the fifth absorption in 48 hours and the first that is not a closed form.** The previous
four fit the standing lesson — a closed form is a coordinate on a densely sampled orbit. This one is
a structure theorem, absorbed by a 1995 textbook that was already in `projects/library/`, in the
chapter the gate itself named. The lesson needs widening: **a novelty gate that prescribes a
computation when the source is on the shelf tests the wrong thing.** I ran an entrywise numerical
comparison when the honest instrument was `pdftoppm`.

## A correction to the brief, in case the framing recurs

The brief calls this a matched pair of bounds and assigns the difficulty backwards:

> "'Divides `d_ν`' is the easy half and is a consequence of Cramer plus `det`. 'Is not a proper
> divisor' is the content."

Cramer bounds denominators by `det M^(L) = ∏_{ℓ(ρ)≤L} d_ρ`, not by `d_ν`. At `n=4` that is **96**,
while the claim for `ν=(4)` is that the denominator is **1**. Cramer cannot deliver the divisibility
half at all. The divisibility half is the substantive one; the sharpness half is one line, because
the diagonal entry `1/d_ν` lies in row `ν`.

## What survived, and why it is not in Example 10

Macdonald gets `u_μ | L_{λμ}` because `Aut(μ)` acts on `E_{λμ}` **freely** (every such `f` is
surjective). On the **inverse** side the same kind of group acts on `Π_{ℓ(λ)}` **non-freely**, and
the fixed points are exactly where the denominators live. That is the gap his argument leaves.

**Theorem.** `Aut(λ)` acts on `Π_k` preserving each fibre `F_{λν} = {π : λ^π = ν}` and preserving
`μ(0̂,π)`. Hence

> `(M⁻¹)_{λν} = Σ_O μ_O / |Stab(O)|` over `Aut(λ)`-orbits `O ⊆ F_{λν}`,

and `den` divides `lcm_O (|Stab(O)| / gcd(|Stab(O)|, μ_O))`.

Consequences, all exact:
- a **single-orbit fibre gives the denominator exactly** — 351 of 396 support entries for `n≤8`;
- diagonal: `F_{λλ}={0̂}`, entry exactly `1/d_λ`;
- first column: `F_{λ,(n)}={1̂}`, entry exactly `(-1)^{k-1}(k-1)!/d_λ`, denominator exactly
  `d_λ/gcd(d_λ,(k-1)!)` — a **proper** divisor of `d_λ` in 38 of 66 rows. At `λ=(1^n)` this is the
  classical `[p_{(n)}]e_n = (-1)^{n-1}/n`.

So the exact-denominator theorem reads like a statement about a row but is carried by single
entries: **267 of 396 support entries have a denominator strictly below `d_λ`.** It is not an
entrywise statement, and I had assumed the diagonal was the unique attaining entry — it is not, in
24 of the 42 non-vacuous rows.

## Honesty on the instruments

- **Falsifiable set counted before the verdict.** Tightness of the bound is an *identity* on
  single-orbit fibres, so only **45** of 396 entries can falsify it, not 396. On those 45: 0
  violations, 39 tight, 6 loose. The headline is 39/45, not 390/396.
- **Planted control fires.** Substituting the off-by-duality quantity `|O|` for `|Stab(O)|` gives
  **248** violations against **0** for the true stabiliser. The test reads the stabiliser order.
- **An instrument fault I nearly banked.** My `Aut(λ)` generator appended `tuple(cur)` where `cur`
  was a `dict` — that yields the dict's *keys*, so every "group element" was the identity and every
  orbit had size 1. My guard `assert len(G)==d_λ` passed the whole time: **the cardinality was right
  while the content was wrong.** First run reported the bound loose in 95 of 204 cases; corrected,
  6 of 396. Fixed with `len(set(G))` plus an entrywise check that each `w` fixes the part sequence.
- The LaTeX error counter was itself tested against a planted bad macro (1 error, 1 undefined ref)
  before its zero reading on the real document was believed.

## Open, and stated precisely

All six loose instances have **exactly two orbits** and miss by a factor of **exactly 2**:
`(2,2,1,1,1)→(3,2,2)`, `(2,2,2,1,1)→(6,2)`, `(2,2,1,1,1,1)→(3,2,2,1)`, `(2,1⁶)→(6,2)`,
`(2,1⁶)→(4,2,2)`, `(2,1⁶)→(3,3,2)`. Every loose `λ` contains a part 2 beside 1s. The uniform factor
2 is unexplained and six instances is too few to lean on. I did not chase a closed form for
`gcd(d_λ, Σ_{π∈F} μ(π))`, because by this week's evidence a closed form is a coordinate.

## Not attempted, with the reason

Q394 (is `M^(L)` a signed count over orientations of a Hasse diagram, à la `2610.08500`?). That
paper sits at **`agent-summary` with zero locators**, this was a no-browsing session, and the brief
itself carries the caveat "I have not verified their path graph is the Hasse diagram of anything."
Any claim resting on it today would have been unsourced. It needs a deep-read first.

## What I'd ask you

The absorption rate is the signal, not the individual results. Five in 48 hours, and today's was
sitting on my own disk. I think the fix is a cheap discipline rather than a better gate: **before
any novelty claim in classical symmetric-function territory, open Macdonald at the relevant chapter
and read the Examples, not just the numbered results.** Four of today's six absorbed items are in
*Examples*, which is exactly where I never look — and my own source entry already warned that
example numbering is edition-risky, which I had read as "cite equations instead" rather than
"the Examples carry content."
