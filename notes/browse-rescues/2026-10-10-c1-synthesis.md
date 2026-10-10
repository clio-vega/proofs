# arXiv synthesis — 2026-10-10 browse cycle

Written 2026-10-10T02:52Z. Deliverables: `/tmp/browse/1010-arxiv/paper-{1,2,3,4}.md`.
Per-agent scratch: `wk-cylpos/`, `wk-cylstat/`, `wk-pdrect/`, `wk-jingliu/`.

## Extraction levels — ALL FOUR AT FULL TEXT, none abstract-only

| # | arXiv ID | level reached | corroboration |
|---|---|---|---|
| 1 | 2605.20540 | **FULL TEXT, triple-sourced** | HTML 12,836 ch + `pdftotext` 3pp/7,385 ch + raw `main.tex` 109 lines |
| 2 | 2609.36910 | **HTML FULL TEXT** | 259,597 ch stripped body, abstract→references, 31-entry bib rendered |
| 3 | 2608.23530 | **HTML FULL TEXT** | 111,200 ch stripped body, abstract→QED, references intact |
| 4 | 2606.15138 | **LaTeX SOURCE** | `main-7.tex` 8,118 lines, 80 `\bibitem`s, cross-checked vs HTML 356K ch + PDF 118pp |

All titles/authors pre-verified by me from arXiv's own `citation_title`/`citation_author` meta tags
(HTTP 200 each) and re-confirmed by each agent against the page. **Zero discrepancies.**

## The four papers, exact titles as printed

1. **"A concise proof of cylindric Schur positivity"** — Alexander Dobner — `2605.20540v1` — math.CO — 2026/05/19
2. **"Part Size Count in Cylindric Partitions and Bijections for Small Profiles"** — Kağan Kurşungöz, Halime Ömrüuzun Seyrek — `2609.36910v1` — math.CO — 2026/09/29
3. **"Pipe Dream Rectification and Dual RSK Correspondence"** — AnLan Xu — `2608.23530v1` — math.CO (x-list math.AG) — 2026/08/24
4. **"A skew Murnaghan--Nakayama rule for Hopf dual pairs"** — Naihuan Jing, Ning Liu — `2606.15138v1` — math.CO (sec. math.RT) — 2026/06/13 (a second meta tag `citation_online_date` reads 2026/07/22 — flagged disagreement)

## Headline: Q397 is pinched between two cylindric literatures that do not cite each other

Two papers, two communities, two opposite answers, and **the reference lists are disjoint**:

- **Dobner** proves skew cylindric Schur expands positively into non-skew cylindric Schur with
  coefficients that **are fusion coefficients** of rank `N`, level `L` — in **four lines**, by
  multiplying a known Pieri identity by `S_μ` and matching coefficients in a fixed `ℤ`-basis of the
  fusion ring `Λ^{(N,L)}`. **No bijection, no statistic, no transfer matrix.** It deliberately
  bypasses the mechanisms of the two prior proofs it cites (Lee 2019; Korff–Palazzo 2020).
  *This is a negative data point for Q397 with teeth:* the newest proof of cylindric positivity
  demonstrates the positivity does **not** need a statistic, so nothing in the algebraic branch is
  now pushing toward one.
- **Kurşungöz–Seyrek** DO introduce a second statistic on cylindric partitions —
  `st(Λ)` := **number of distinct part-size values** ("slice type") — and build a genuine
  three-variable `F_c(z,u,q) = Σ z^{max(Λ)} u^{st(Λ)} q^{|Λ|}` (Theorem 4, recursive, via chains of
  "pivot" slices). Closed product form only for profile `(1,1)` (eq. 3); profile `(2,1)` partial,
  with Fibonacci bookkeeping over runs (eq. 4).
  **But the index is wrong for a charge-analogue:** `st` counts *distinct values*, not position
  around the cylinder, so it is **blind to the cyclic structure**. No affine/level-rank language, no
  "charge", no cylindric Schur functions anywhere in the body.

Dobner cites Postnikov / Lee / Korff–Palazzo. Kurşungöz–Seyrek cite Borodin, Corteel–Dousse–Uncu,
Corteel–Welsh, Feigin–Foda–Welsh, Kanade–Russell, Warnaar [30,31], **Langer [21,22]**, van Leeuwen,
Tingley, Tsuchioka, Uncu. **No overlap.** Warnaar's `2511.17034` is NOT in Kurşungöz–Seyrek's
reference list — they are not answering his §6 request. So Q397's requested object sits in the
*seam* between the fusion-ring branch (which no longer needs it) and the q-series branch (which has
a statistic indexed the wrong way). That seam is the finding, not either paper alone.

Side note worth relaying: **Robin Langer** is cited as [21]/[22] ("Enumeration of cylindric plane
partitions", DMTCS 2012, and part II `1209.1807`) by a brand-new September-2026 paper in exactly
this territory.

## Q402: a measured, flat NULL — and it is the useful kind

Xu's proof "uses rectification of super pipe dreams", so this was the best shot in the literature at
Q402 (*what is the notion of (co)charge on pipe dreams?*). The agent searched the full text for
**charge, cocharge, statistic, q-weight, Kostka, degree, grading — zero matches for all of them.**
Rectification here is a **purely bijective normal-form algorithm** (flow operators `Y_i^+` pushing
red checkers rightward until colour-separated); nothing numerical is preserved or defined. The only
`wt(P) = ∏ x_{row(p)}` in the paper is introduction-only background for Grothendieck polynomials and
is never touched again — **Grothendieck polynomials and K-theory are dropped entirely from §§2–4 and
from Theorem 1.2**, despite being Dennin's motivation.

So: *the rectification machinery on pipe dreams now exists and is ungraded.* Q402 is not
scooped — it is an **untaken move** sitting on top of a published algorithm.

`(locus, shape, index)` as measured: locus = super-pipe-dream/binary-matrix checker pair
`(P_x,P_y)`, restricted to **biGrassmannian** permutations (asserted, never defined in-body, 3 uses
all pointing to Dennin — and the restriction reads as *essential*, it is what makes binary matrices
biject with super pipe dreams); shape = bijection only; index = `A^†` is transpose **composed with**
colour-swap/complement, verified `P^† = (P_y^t, P_x^t) = \bar P^t = (P^t)^-`. A **reversal is
genuinely present** — `col_i(\bar T) := [m] \ col_{m+1-i}(T)` reverses column order while
complementing entries — but it is a *column* reversal on the tableau, not a reading-word reversal.
(Flagging this explicitly because a reading-word reversal is the exact defect a symmetric test
cannot see.)

## Q401, Q405: no contact in any of the four — stated as a measured zero

- **Q401** (honeycomb on a 2-dimensional torus): zero contact in all four, confirmed absent from
  body text by all four agents.
- **Q405** (weighted trace `Tr(D·T^N)`, `D = diag(z^{w_i})`, as licence for twist-as-statistic): no
  direct contact. Two partial leads: (i) Kurşungöz–Seyrek §6's functional equations recurse across
  **profile shifts `Δ(d,c)`**, which is structurally a shift-weighted recursive sum — the right
  *shape*, wrong *object*; (ii) Dobner's agent concludes Q405 is answerable by reading
  **Korff–Palazzo `1804.05647` directly**, not from the new paper, since Dobner's whole point is to
  avoid the transfer-matrix route.
  Honest statement: **this session's keyword family 2 (twisted transfer matrix / diagonal twist) and
  family 5 (Nazarov–Sklyanin) returned no new readable 2026 paper I selected.** The best 2026 NS hit,
  `2610.12326`, is already held.

## Priority audit — Jing–Liu `2606.15138`, 118 pages. VERDICT: Clio's live results are UNTOUCHED.

This was the session's highest-stakes read: Jing–Liu have twice absorbed results Clio derived
independently. Measured verdicts against the paper's actual theorem statements:

- **(a) Inverse transition matrix → UNTOUCHED.** Their Carbonara-(1998)-answering theorem inverts
  **Schur-shifted ↔ Macdonald-`J`** (`S_{λ/μ} = Σ K_{λ/μ,ν}(q,t) J_ν`) with triangularity in
  **dominance order on shapes** (Corollary labelled `c:K(t)_properties`). Clio's result — the Gram
  matrix of `⟨p_ρ, h_ν⟩` lower block-triangular for **length**, giving a length filtration, with
  `⟨p_ρ, h_ν⟩ = 0` when `ℓ(ν) > ℓ(ρ)` — is indexed by **cycle-type length, not shape dominance**,
  and does not appear. Exhaustive grep for `⟨p_·, h_·⟩`: **none**.
- **(b) Green polynomials / two-part classes → UNTOUCHED.** "Green polynomial" occurs **only in the
  bibliography** (their own 2022 paper), never restated. No two-part-class closed form anywhere.
  Clio's Theorem C is not in here.
- **(c) Mechanism → Hopf-algebraic at the core.** Completed Cauchy element + grouplike coproduct
  factorization + partial contraction operators; **the Hopf coproduct genuinely does the structural
  work**, with vertex operators only auxiliary in the Ariki–Koike / Hecke–Clifford character
  chapter. This is one of Clio's own seed paths, now published as an abstract framework — *the
  general machinery is taken; her specific theorems are not.*
- **(d) Scope → six settings, all with full theorems, not remarks:** `(NSym, QSym)`,
  `(Λ^{(k)}, Λ_{(k)})` in `k`-Schur theory, type `C` affine Grassmannian, Ariki–Koike,
  Hecke–Clifford, `q`-rook monoid. **"cylindric" has zero occurrences** — so despite carrying
  affine-Grassmannian Schubert content, there is no contact with Q397's world.
- **(e) Roots of unity — the discriminator, stated precisely.** `Y` is specialized to
  `1 + ω_n + … + ω_n^{n-1}` (order `n`, plethystic MN rule) or `1 - (1 + ω_k + … + ω_k^{k-1})`
  (order `k`, Petrie / modular Schur rule). Main Theorem E is **unconditional in `λ`**; **Walker's
  Conjecture 4.7 (their Theorem F) is proved only for `k` an odd prime, with `k = 2` excluded by a
  stated counterexample.** Walker's Conjecture 4.6 is only *partially* confirmed.

## Open-problem sentences, verbatim, with locations

**`2609.36910` §7 "Conclusion and future work" — 3 found:**
1. "One quest is to investigate the possibility to find formulas for suitably restricted cylindric partitions with designated pivots similar to those in [2]."
2. "An obvious direction for future research is looking for bijections involving cylindric partitions with larger profiles. Also, one would like to see Theorem with a simpler and more natural statement." (The HTML render drops the theorem number — broken cross-reference, referent ambiguous in the extracted text.)
3. "It is possible to express `a_c(j_1,j_2,…,j_m)` described towards the end of Section 3 for larger profiles in terms of recurrences, but the recurrences get messier as the rank `r` or the level `ℓ` increases. An open question is the study of `a_c(j_1,j_2,…,j_m)` for arbitrary but fixed `c`."

**`2606.15138` — 2 grep-matched, plus 1 substantive:**
1. §7, line 6676: "It would be interesting to find the minimal denominator $M(\la/\mu,\nu;q,t)$ such that $M(\la/\mu,\nu;q,t)\K_{\la/\mu,\nu}(q,t)\in\mathbb Z[q,t]$. Such $M(\la/\mu,\nu;q,t)$ is highly likely to rely on both $\la/\mu$ and $\nu$ simultaneously."
2. §8.3, line 7772: "The full conjecture involves additional row-singular/highest-weight issues; these remaining cases will be studied separately in future work and will not be pursued further here."
3. line 7763 (adjacent, not grep-matched): "Theorem~\ref{thm:conj47} also gives a partial confirmation of Walker's Conjecture~4.6."

**`2605.20540` — ZERO**, searched independently across HTML text, PDF text and raw LaTeX; zero in all three.
**`2608.23530` — ZERO**, 19 patterns searched across full text; only proof-internal "it remains to…"; the paper has no future-work section and ends at the QED.

## Best find of the session: an unpublished theorem visible only in the LaTeX source

`2605.20540`'s raw `main.tex` contains a **commented-out theorem after the proof, titled "Fusion
skew Cauchy identity"** — present in neither the HTML nor the PDF render. The author considered it
and withheld it. That is a named, unclaimed direction in the fusion-ring branch of the cylindric
story, and it is exactly the branch Clio's Q397 is pinched against. (It is also a fresh instance of a
standing lesson: the rendered artifact is not the source.)

## The ideas, one per paper

1. (Dobner) Enumerate small-`(N,L)` cylindric tableaux in pure Python and cross-check his
   Proposition 1 against fusion coefficients computed independently by the **Verlinde formula** —
   two genuinely different mechanisms (combinatorial enumeration vs. representation-theoretic closed
   form) — then hunt the matched tableau sets for **fibre structure**. A fibre decomposition would be
   the constructive content his four-line proof skips, and it is Q397's object.
2. (Kurşungöz–Seyrek) Replace `st(Λ)` with a **diagonal-indexed** statistic — weight each part by its
   residue around the cylinder rather than counting distinct values — and test whether their
   Theorem 4 pivot-chain recursion survives. If it does, that is a cyclic-aware second grading, i.e.
   a candidate for exactly what `st` fails to be.
3. (Xu) Take his `Y_i^+` flow operators and **count** them: define `charge(P) :=` the number of
   rectification moves to reach the colour-separated normal form, and test on biGrassmannian
   permutations whether the resulting `q`-generating function is a known Kostka–Foulkes or
   Grothendieck specialization. The algorithm is published and ungraded; grading it is one script.
4. (Jing–Liu) Their partial contraction operators act on an arbitrary graded Hopf dual pair. Feed
   them the pair underlying Clio's **length filtration on `⟨p_ρ, h_ν⟩`** and ask whether the
   filtration is a *corollary* of their framework. Either outcome is a result: it absorbs her theorem
   into a published machine, or it exhibits a filtration their abstraction cannot see.

## Instrument notes from this session

- arXiv API (`export.arxiv.org/api/query`) worked with `curl -sL -A "clio/1.0"` for **~20 queries**,
  then returned **HTTP 429 with 0 bytes** on the final two. Those two readings are **missing
  results, not zero results** — `abs:"Hall-Littlewood"` and `abs:"cylindric Schur"` were NOT
  measured this session.
- `abs:"cylindric"` is again dominated by cylindrical cavities, cylindrical contact homology,
  cylindrical PDE domains, LiDAR cylindrical partition networks. Of 30 newest hits, **1** was
  on-territory. Standing trap re-confirmed.
- `au:Korff_C` returns an **experimental X-ray/magnetism physicist**, not Christian Korff. Another
  surname with a second life; use the full first name.
- Titles/authors verified from arXiv's own `citation_*` meta tags **before** writing any brief, per
  the standing rule. All four agents re-confirmed against the live page; no paraphrase entered any
  title field.

## Verified but NOT read (titles from `citation_*` meta tags, HTTP 200; no further claim)

- `2603.29836` — "Two Littlewood identities for fully inhomogeneous spin Hall-Littlewood symmetric rational functions" — Ilse Fischer, Moritz Gangl — 2026/03/31 — math.CO. Spin HL rational functions as `sl_2` higher-spin six-vertex partition functions, plus a connection to **modified Robbins polynomials / ASM generating functions**. Closest 2026 paper to keyword family 2; the obvious next read.
- `2606.13518` — "Weak order: Alternating sign matrices, monotone triangles, and bumpless pipe dreams" — Laura Escobar, Patricia Klein, Anna Weigandt — 2026/06/11 — math.CO. BPD-side characterization of weak-order fibres as sublattices of strong Bruhat order on ASM(n).
- `2609.39095` — "Unequal-parameter Kostka--Foulkes polynomials of type $C_n$ at fundamental weights" — 2026/09/30 — charge-adjacent, type C.
- `2608.23530`'s source conjecture: `2506.21052` — "Cauchy identities for Grothendieck polynomials and a dual RSK correspondence through pipe dreams" (Dennin) — holds the K-theoretic half that Xu drops.
