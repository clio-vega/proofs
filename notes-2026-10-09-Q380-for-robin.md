# Q380 closed: charge is never constant on a nontrivial Kostka–Foulkes fibre

**Clio, 9 October 2026 — prove session.**
Paper: `https://github.com/clio-vega/proofs/blob/main/2026-10-09-Q380-monomial-kostka-foulkes.tex`
(8pp, compiles clean; code alongside in `code-1009-q380/`.)

## The result, in one line

Q380 asked: for which `(ν,μ)` with `#SSYT(ν,μ) ≥ 2` is `K_{ν,μ}(t)` a monomial?

**There are none.** Macdonald III (6.5)(ii) says `K_{λμ}(t)` is *monic* of degree
`n(μ)−n(λ)`. A monomial therefore has coefficient **1**; and at `t = 1` that coefficient
is the Kostka number. So

> `K_{λμ}(t)` is a monomial ⟺ `K_{λμ}(t) = t^{n(μ)−n(λ)}` ⟺ `#SSYT(λ,μ) = 1`.

The content form, which is the one worth remembering: **charge takes at least two values
on every fibre with at least two tableaux**, because the maximal-charge tableau is unique.

I make no priority claim — this is a two-line deduction from a 1995 textbook and I expect
it is folklore. I didn't look. The brief told me the value here was use-value, and it was
right.

## Why it mattered (this is the part I'd like you to check)

It is the amendment owed by **[A3]** of `2026-10-05-c3-migration-length-grading.tex`.
That note proved: a monotone potential makes the fibre generating function a monomial —
and then observed, correctly, that at `c = 1` the conclusion is vacuous, with 360 of 361
swept fibres at `c = 1`.

For charge on Kostka–Foulkes fibres the situation is worse than "the sample was bad":
**the conclusion is *equivalent* to `c = 1`.** The no-go has no non-vacuous instance there
at all. Contrapositive, and this is the usable statement:

> whenever `#SSYT(λ,μ) ≥ 2`, charge is **not** a potential difference for any monotone
> process on that fibre.

Q370 asked *why* charge evades the potential no-go — by reversibility, by branching, or
because `K_{λμ}(t)` isn't a monotone-process fibre generating function at all. **The third
horn, and it is now a theorem rather than an assertion.** Q370 was gated on Q380 precisely
for the instances; the instances turn out not to exist, which answers it.

## The boundary, which is where the open problem now sits

Since monomial ⟺ Kostka number 1, Q380's boundary *is* the boundary of `{K_{λμ} = 1}`.
I did not close that, but I pinned it down:

- For fixed `λ`, `M(λ) = {μ : K_{λμ} = 1}` is an **order filter** in dominance, so it is
  described by its **minimal elements** — **at most two** for every `λ` with `|λ| ≤ 10`.
- Two exact reductions, both proved: **R1** (`μ₁ = λ₁` ⟹ strip row 1 of `λ` and the largest
  part of `μ`) and **R2** (`ℓ(μ) = ℓ(λ)` ⟹ strip column 1 of `λ` and one from every part of
  `μ`). Since `λ ⊵ μ` forces `μ₁ ≤ λ₁` and `ℓ(μ) ≥ ℓ(λ)`, these are the two opposite
  extremes. They leave **10 of 1044** irreducible pairs with `N ≤ 10`; `GL_r`
  complementation kills 4 more, leaving **6**.
- Corollaries of R1+R2: single-row `λ`; the whole two-row/two-row slice (so the brief's
  suspicion about `lem:mono` of the 10-07 two-part paper was right — it is a singleton
  fibre, vacuous); and every **adjacent transfer** `μ = λ − r e_k + r e_{k+1}`.
- **Two proved no-gos on the shape of any answer**, which is the part I'd most like a second
  opinion on: *no criterion depending only on `λ−μ`* can work (`(1,1,−2)` gives `K=1` at
  `((3,3),(2,2,2))` and `K=2` at `((4,3),(3,2,2))`), and *none depending only on the
  dominance-slack vector* (`(2,1,0)` gives `K=1` at `((4,1,1),(2,2,2))` and `K=2` at
  `((5,1),(3,2,1))`). Those were the two shapes of answer I reached for first.

**Question for you:** is there a known characterisation of `K_{λμ} = 1`? If it exists I'd
rather cite it than keep measuring. The filter structure plus "at most two minimal
elements for `|λ| ≤ 10`" smells like someone has done this.

## One gap I am citing rather than proving

`Prop. 8` (monotonicity of Kostka numbers in dominance, `ν ⊵ μ ⟹ K_{λν} ≤ K_{λμ}`) is
classical and I only *sketch* the `sl₂`-string argument. Verified on 6750 triples, `N ≤ 9`,
0 violations. It is used for the boundary structure only — **never** for the main theorem.

## Instrument note, because it nearly cost the paper

The OCR text layer of Macdonald p.239 column-shifts the `n = 5` `K(t)` matrix and yields
`K_{(2,2,1),(2,1,1,1)} = t` — a plausible monomial on a two-element fibre, i.e. a
counterexample to my own theorem, in the right notation, from the right page. The real
entry is `t + t²`. I caught it only because I had computed that fibre by hand before
trusting the page; confirmed by rendering at 200 dpi. An OCR failure is not noise; it is a
well-formed false statement.

Verification: two mechanism-disjoint engines (Lascoux–Schützenberger charge; Weyl
alternating sum over the `t`-analogue of Kostant's partition function — no tableau, no
charge in it) agree on all 53 dominance pairs `N ≤ 5`, and both agree with Macdonald's
printed tables read off page images. Census: all 1719 nonempty dominance pairs with
`N ≤ 10`, of which **1245 are non-vacuous**; 0 monicity failures, 0 degree failures,
0 pairs in the target set. Planted control: collapsing charge to a constant makes the
detector report a monomial on a 3-element fibre, so the empty target set is a measurement
and not a tautology.
