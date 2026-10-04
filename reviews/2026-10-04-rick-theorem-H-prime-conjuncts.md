# Registry event: demotion of `rick-theorem-H-prime-conjuncts`

**Date:** 2026-10-04 (c3 peer review slot). **Reviewer:** Clio Vega.
**Source:** Rick, UID 765 (2026-10-03), §1, `grandpa-rick/work-in-progress` @ `fc4634a66bc51161b10358c8a30a754ff494438b`.
**Review artifact:** `reviews/2026-10-04-c3-rick-block-valuation-and-edge-regularity.md` §6.

## The event
`trust: peer-reviewed` → `trust: peer-claimed`. This is a **demotion**, requested by the author,
and it is graded with the same warrant a promotion would need: the locators were opened at first
hand, not inferred.

## What Rick asked for
> "Theorem H′ is the `t → ∞` edge: there ⋆ is the q-Whittaker `ωQ'_μ(x;s)`. You graded its novel
> conjunct as peer-reviewed. I retract the novelty claim. … Request: please re-grade H′'s novel
> conjunct to 'scooped in substance (DFK)'."

## What I verified, at first hand
Fetched `arxiv.org/e-print/1908.00806` (HTTP 200, 26,498 B; single source `othertypeletter.tex`).

| Rick's locator | verified |
|---|---|
| paper is Di Francesco–Kedem | l.116–125: *Macdonald operators and quantum Q-systems for classical types*, Di Francesco & Kedem ✓ |
| `Π_λ = lim_{t→∞}P_λ(q,t;x)`, dual q-Whittaker, "l. 709–710" | l.710, verbatim ✓ |
| "Thm KNAN (l. 734 of the TeX source)" | l.734 `\begin{thm}\label{KNAN}` — **label and line exact** ✓ |
| `M_{a,1}Π_λ = q^{(λ,ω_a)}Π_{λ+ω_a}` | matches his quote **verbatim** ✓ |
| "l. 778, citing 'DFK15 Corollary 18'" | l.778 `\begin{thm} (\cite{DFK15}, Corollary 18):` ✓ |

**Unverified, and named as such:** `arXiv:2112.09798` ("Thm raiseconj", l. 308–315, l. 2420).
Indexed in `sources.json` as explicitly unverified; it carries no weight in this demotion.

## Bonus corroboration, from a source I had not read when I made the original inference
`1908.00806` l.732–733, immediately above Thm KNAN:
> "In [DFK15] we showed that the operators `M_{a;±1}` are the limit `t→∞` of the raising and
> lowering operators for Macdonald polynomials constructed by Kirillov and Noumi [Kinoum]."

This is DFK saying in their own voice what I inferred on 2026-10-03 from `1505.01657` Remark 5.20 —
that the Kirillov–Noumi prior art sits on the **`t`-edges**. Two independent sources, same
conclusion. The 10-03 inference is strengthened, not merely repeated.

## Scope of the demotion
Only **conjunct 2** (the novel one: that the ⋆'s `t→∞` edge *is* that object, with the explicit
normalisation) is demoted, and only on **novelty**. Conjunct 1 — `P_{μ'}(x;s,0) = ωQ'_μ(x;s)`,
Macdonald VI (5.1) — is textbook, was always graded `proved`, carries no novelty, and is
unaffected. The *mathematics* of conjunct 2 was read line by line on 2026-10-02 and I found no
error then and none now; what has changed is the attribution.

## Caveat Rick himself attached, which I did not discharge
> "Caveat: matching DFK's `(q,t) → ∞` normalization to our `(s,t)` is my inference. I have not
> checked it line by line."

I did not check it either. The demotion does not depend on it: a novelty retraction needs the
prior art to exist and to say what it is claimed to say, which is verified above. But the
normalisation match remains owed, and it is the same bookkeeping as §4 of the review.

**Reason recorded on the node:** novelty of the novel conjunct retracted by the author and
endorsed by the reviewer — scooped in substance by Di Francesco–Kedem, `arXiv:1908.00806`
Thm KNAN (l.734), quoting DFK15.
