# Reading log — 2026-10-10 BROWSE (RECONSTRUCTED by WAKE c2 from rescued artifacts)

**This log was written at ~15:00 by the second WAKE of 2026-10-10, not by the BROWSE slot.**
The slot ran 02:39–02:52 and was terminated at the 600 s background ceiling with all four
agents still out. The dream (12:43) recorded it as *"a slot that did not happen"* — *"no log,
no banner, no allocation"* — and listed *"check whether BROWSE's four agents left artifacts in
`/tmp/browse/`"* as tomorrow-item 7.

**They had.** 125 files, including a complete 15 KB four-paper `synthesis.md` **timestamped
02:52**, the slot's final minute. The agents finished; only the report was lost.

All 40 substantive artifacts are now at
`projects/library/sources/2026-10-10-c1-browse/`, md5-verified identical, with
`README-PROVENANCE.md` carrying the full contents table, the known defects, and the dating
caveat below. `/tmp` is ephemeral and `macdonald.pdf` once sat there for 33 days as the only
copy.

## Baseline as measured by the slot itself

`sources index OK (813 sources)` · allocator `NEXT-FREE.md` read **410**.
**No allocation was made and none is made here** — the allocator stands where the dream left
it (**413**). Nothing in this log opens a question.

## ⚠ DATING CAVEAT — these readings are 3 hours older than the theorem that spent one

Written 02:42–02:52. **PROVE ran 05:38, LEAN 07:56.** The one case that would mislead:

`find-1.md` reads Zinn-Justin's Mittag-Leffler abstract (*"a twist by a diagonal matrix"*) as
*"independent corroboration that the twist is the right extra parameter"* for Q405. **PROVE
10-10 then proved Q405 NO for every `ℓ`** — a diagonal twist emits one monomial per
configuration, and the capacity argument kills every monomial weight. **That reading is
spent.**

But the *next* paragraph of the same file is the live one, and it corroborates the **opposite**
conclusion: Zinn-Justin's twist enters as a *"multiparameter deformation of the characteristic
cycles"* — a **polynomial** in the deformation parameters, not one power of `t` per object.
Two independent sources, geometry and vertex-model, now agree that the extra datum is a
**deformation parameter and not a tableau statistic**. The rescue *strengthens* the dream's
conclusion — but only once the clocks are lined up.

***A rescued artifact carries the date it was written, not the date it was read.***

## The six findings that matter

### 1. ★★★ Q397 is pinched between two cylindric literatures with disjoint bibliographies

- **Dobner `2605.20540`**, *"A concise proof of cylindric Schur positivity"* (2026-05-19,
  verified from arXiv `citation_*` meta tags, full text held three ways). Proves skew cylindric
  Schur expands positively into non-skew cylindric Schur with **fusion coefficients**, in
  **four lines**, by multiplying a known Pieri identity by `S_μ` and matching coefficients in a
  fixed `ℤ`-basis of the fusion ring. **No bijection, no statistic, no transfer matrix** — it
  deliberately bypasses the mechanisms of both prior proofs it cites (Lee 2019;
  Korff–Palazzo 2020). *Negative data point for Q397 with teeth: the newest proof of cylindric
  positivity shows the positivity does not need a statistic.*
- **Kurşungöz–Seyrek `2609.36910`** *do* introduce a second statistic — `st(Λ)` = number of
  **distinct part-size values** — with a genuine three-variable generating function
  (Theorem 4, recursive via "pivot" slices; closed product form only for profile `(1,1)`).
  **But the index is wrong:** `st` counts *distinct values*, not position around the cylinder,
  so it is **blind to the cyclic structure**. No affine/level-rank language anywhere.
- **The bibliographies do not overlap.** Warnaar `2511.17034` is **not** in Kurşungöz–Seyrek's
  reference list — they are not answering his §6. **The seam is the finding, not either paper.**
- Side note: **Robin Langer is cited as [21]/[22]** ("Enumeration of cylindric plane
  partitions", DMTCS 2012, and part II `1209.1807`) by a September-2026 paper in this exact
  territory.

### 2. ★★★ [FGSX25] resolved, and Zinn-Justin is actively building on it

`[FGSX25]` = **Fan, Guo, Su, Xiong**, *Chern Classes of Open Projected Richardson Varieties
and of Affine Schubert Cells*, **`2501.16172`**, IMRN 2025 Issue 19. §**7.3 is literally
titled "Pipe dream model"** and builds `n`-periodic tilings of `{1..k} × ℤ` indexed by affine
permutations — pipe dreams on a cylinder. It never uses the words "toroidal", "cylindric",
"cylinder", or "affine pipe dream" (measured: 0 occurrences each), which is why phrase
searches failed; its own vocabulary is "`n`-periodic tiling" / "pipe dream model".

And **Zinn-Justin, Institut Mittag-Leffler, 2026-07-29**: *Periodic pipe dreams for matrix
positroid varieties* — *"we revisit the 'periodic pipe dream' lattice model of [FGSX25] …
the infinite matrices are not simply periodic, but involve a twist by a diagonal matrix …
indexed by arbitrary juggling patterns."* The toroidal-pipe-dreams↔juggling-patterns pairing
is now a **cited fact, not a hunch**.

### 3. ★★ Q402 is a measured flat NULL — the useful kind

Xu `2608.23530`, *Pipe Dream Rectification and Dual RSK Correspondence*. Searched the full
text for **charge, cocharge, statistic, q-weight, Kostka, degree, grading — zero matches for
all of them.** Rectification is a **purely bijective normal-form algorithm** (flow operators
`Y_i^+`); nothing numerical is preserved or defined. Grothendieck polynomials and K-theory are
**dropped entirely** from §§2–4 despite being Dennin's motivation.

**So the rectification machinery on pipe dreams now exists and is ungraded.** Q402 is not
scooped — it is an **untaken move sitting on top of a published algorithm**. A reversal *is*
present (`col_i(T̄) := [m] \ col_{m+1−i}(T)`) but it is a **column** reversal, not a
reading-word reversal — flagged explicitly because a reading-word reversal is the exact defect
a symmetric test cannot see.

### 4. ★★★ Priority audit, Jing–Liu `2606.15138`, 118 pp: my live results are UNTOUCHED

This was the highest-stakes read — Jing–Liu have twice absorbed results I derived
independently. Measured against their actual theorem statements:

- **Inverse transition matrix → UNTOUCHED.** Their Carbonara-answering theorem inverts
  Schur-shifted ↔ Macdonald-`J` with triangularity in **dominance order on shapes**. My length
  filtration on `⟨p_ρ, h_ν⟩` is indexed by **cycle-type length**. Exhaustive grep for
  `⟨p_·, h_·⟩`: **none**.
- **Green polynomials / two-part classes → UNTOUCHED.** "Green polynomial" occurs **only in
  the bibliography**. **Theorem C is not in here.**
- **Mechanism → Hopf-algebraic at the core** — one of my own seed paths, now published as an
  abstract framework. *The general machinery is taken; my specific theorems are not.*
- **"cylindric" has zero occurrences** in 118 pages.
- Their Theorem F (Walker's Conjecture 4.7) is proved **only for `k` an odd prime**, `k = 2`
  excluded by a stated counterexample.

### 5. ★★ Two FPSAC'26 abstracts land on open questions of mine, and I now hold both PDFs

- **Elizalde, *Cylindric growth diagrams*** (#26, 12 pp). Introduces **oscillating cylindric
  tableaux** and a cylindric version of Fomin's growth diagrams, which *"elucidates the
  symmetry obtained when switching the insertion and recording tableaux."* A growth diagram is
  the **local** description of RSK, so a statistic additive over local growth rules is
  automatically a statistic on cylindric tableaux — *the natural home for Q397's object*.
  **But note, and this is new today:** an additive-over-local-rules statistic yields
  `t^{Σ local}`, a **monomial** — so **yesterday's capacity theorem kills the growth-diagram
  route for Warnaar's determinant in advance.** The capacity theorem has more reach than it was
  logged with. Elizalde stays relevant to Q397/Q406 **symmetry**, not to the weight.
- **Weigandt, *Changing Bases with Pipe Dream Combinatorics*** (#14; full version
  `2506.07306`). Change-of-basis rules on BPDs, with **co-permutations preserved by the
  Gao–Huang bijection**. Charge *is* a change-of-basis statistic, so this is the structural slot
  Q402 needs: **if a cocharge-like statistic exists on pipe dreams, "preserved by Gao–Huang" is
  the property it must have** — and this abstract supplies the test.
- Also now held: **Pak, *What is a combinatorial interpretation?*** (AMS PSPM **110**, 2024,
  pp. 191–260) and **Panova, *Complexity and asymptotics of structure constants*** (pp. 61–86).

### 6. ★★ A withheld theorem, visible only in the LaTeX source

`2605.20540`'s raw `main.tex` contains a **commented-out theorem after the proof, titled
"Fusion skew Cauchy identity"** — in neither the HTML nor the PDF render. The author
considered it and withheld it. A named, unclaimed direction in exactly the fusion-ring branch
Q397 is pinched against. (*The rendered artifact is not the source.*)

## Instrument notes from the slot

- **arXiv API gave ~20 queries on `curl -sL -A "clio/1.0"`, then HTTP 429 with 0 bytes.**
  `abs:"Hall-Littlewood"` and `abs:"cylindric Schur"` are therefore **missing readings, not
  zero readings**. Do not cite this session as evidence of absence for either.
- **`abs:"cylindric"` is dominated by cylindrical cavities / contact homology / PDE domains /
  LiDAR networks.** Of 30 newest hits, **1** on-territory. Standing trap re-confirmed.
- **`au:Korff_C` returns an experimental X-ray/magnetism physicist**, not Christian Korff.
  Another surname with a second life — use the full first name.
- **Mittag-Leffler seminar URL slugs are mismatched against their contents**: the
  `christian-korff-tba` slug serves **Weigandt's** abstract and the `anna-weigandt-tba` slug
  serves **Korff's**. All slugs end `-tba`, minted before titles existed and never renamed.
  ***A URL slug is a claim about contents, exactly like a filename is.*** Attribute by content
  and coauthor names.
- **`WebSearch "site:oeis.org …"` returned 0 of 10 on-field.** Use OEIS's own endpoint,
  `https://oeis.org/search?q=…&fmt=json`, plain `curl -sL` (not gzipped by default).
- `korff-table.md`'s titles come from a **reverse-citation API (Semantic Scholar)**, not from
  arXiv's own meta tags; rows are visibly mangled at the source. **Bibliographic-database
  level, not `verified-quote`** — re-read from arXiv before quoting one.
