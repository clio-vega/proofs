# For Robin — PROVE 2026-10-21 c2: G3 closed, and I retracted a hypothesis of my own

**Artifact:** `proofs/2026-10-21-c2-G3-boundary-diagonal-typeB.tex` (8 pp, commit `1b39268`).
**Registry:** `motegi-typeBC-multiplication.json` → `typeB-boundary-diagonal-valuation`
promoted `computed` → `proved`, with five children carrying the lemmas separately.
**Trajectory:** `state/trajectory/clio-prove-G3-20261021c2.jsonl`, 18/18 events, 0 violations.

## The result

`ord_S Ã^B_{I,I} = [m ∈ I]` on `S = {w_m=1}`, all `N`, all `I` — type B, which had stood at
*measured 34/34, not proved* since 10-19. So §§5–7 of the 10-19 S/A paper are now **proved in
both types**. What I actually proved is type-uniform and stronger:

> `ord_S Ã^X_{I,I} = [m ∈ I] · #{a>0 : a·t_m ∈ Δ⁺(X_N) and ε^a = 1}`

Three lemmas: a unit criterion for factors `1 − c·monomial`, then "the ascent block is all
units" (every `T`-factor carries `q`), then "the inversion block has at most one non-unit".

## The thing I think is actually interesting

**"Order exactly one" is reducedness of the ambient root system.** The order counts positive
roots proportional to `t_m` vanishing on `S`. `B_N` and `C_N` are reduced, so that is ≤ 1. I ran
the pre-registered falsifier: in the **non-reduced** `BC_N`, where `t_m` *and* `2t_m` are roots,
the same computation gives **2** at `ε=1`. Since this order feeds a *colength*, a non-reduced
ambient system doubles my boundary colength. The discriminator — `α` and `2α` both roots — is
the same one Q538 is built on, which is why I flagged it rather than just noting it.

It also derives **G4** instead of measuring it: type B's boundary root is `t_m`, whose zero locus
is `{w_m=1}` alone, while `2t_m` gives `{w_m=1} ∪ {w_m=−1}`. And the failure is exact — in type B
at `ε=−1` the order is 0 for every `I`, so (BD) fails on **precisely** the `N·2^{N−1}` = 1,4,12
pairs with `m ∈ I`, each by exactly 1. (The 10-19 G4 measurement was of **(BC)**; this is the
**(BD)** half. One mechanism, both halves.)

## I retracted a hypothesis of my own — and it had no objector

My 10-19 (BD) read "at a point of `S` generic in the sense of `rem:generic`", and that remark
defined genericity as avoiding the hypersurfaces `w_k = ±q^{-1}`. **That is not a hypothesis of
the theorem.** `ord_S` is the valuation at the *generic point* of the divisor, residue field
`k = C(q)(w_j : j≠m)`; those hypersurfaces are codimension 1 *inside* `S`, hence proper closed
subsets of it, hence invisible to a divisorial valuation. I had taken a limitation of my
**instrument** — a numerical sample that landed on `w_3 = 1/q` and reported `ord_S det Φ = 8`
against the correct `4` — and written it into the **theorem**. The diagnosis in 10-19 was right;
the generalisation to a hypothesis was not.

The honest shape of this error: *an over-cautious hypothesis attached to my own statement has no
natural objector either.* I have the mirror of this logged (over-claims in my favour have none);
this is the same asymmetry pointing the other way, and it is harder to notice because it reads
as rigour.

Sharp replacement: `q ≠ ε^{-h}`, needed **only** on the `m ∉ I` half.

## P1 came back SPLIT — third consecutive session, and the sharpest form yet

I pre-registered SPLIT and got it. The two halves of (BD) carry **complementary** hypotheses:

| half | needs | does not need |
|---|---|---|
| `m ∈ I` (value 1) | `h` invertible in the base field | any condition on `q` |
| `m ∉ I` (value 0) | `q ≠ ε^{-h}` | any condition on `h` |

Because `h·t_m` lies in **exactly one** of `Inv(λ)`, `T(λ)`, neither hypothesis is visible from
the other half. Shown as a **lesion**, not argued: setting `q := 1` leaves the `m ∈ I` half inert
(1→1, 129/129) and breaks the `m ∉ I` half (0→1, 129/129), both types, `N ≤ 5`. Downstream
consequence worth having: **a `q → 1` degeneration destroys (BD)'s dependence on `I` entirely**,
turning the boundary colength `2^{N−1}` into `2^N`.

Also a by-product: the type-C half of the unit criterion needs `char ≠ 2` (in char 2,
`1−v^{-2} = (1−v^{-1})²`), while type B is characteristic-free. **The type that was already
proved is the fragile one.**

## Two process findings, both the same shape

1. **The gap was a reading gap, and the locator said so.** My 10-19 note cited **L1386** for
   "the type-B diagonal". L1386 is the *unexpanded* formula — **verbatim identical to type C's at
   L1078**, which is exactly why it looked sufficient. The explicit type-B root products are at
   **L1305–1335**, eighty lines earlier, inside the subsection Motegi opens by omitting all
   type-B proofs. I had carried a parameterised *guess* instead ("the boundary term `(1−w_i^{-h})`
   occurs") and **the guess was right** — which is precisely why it survived a year of use
   unexamined. A correct guess in a cited slot is invisible to every check that asks whether the
   conclusion is true.
2. **`tracecheck` caught a real circular verification, again correctly.** I offered the
   transcription-vs-engine comparison as evidence *for the transcription*; the comparison
   descends from it. The fix also relocated the plant census: it certifies the **instrument**,
   not the transcription, because the plants are built *from* the transcription. Fourth instance
   of this class, and the first where the right encoding changed what the evidence was *about*.

## Honest limits

- I did **not** re-derive Motegi's type-B diagonal. He omits that proof ("essentially the same as
  the corresponding type C arguments… we omit these proofs", L1281) and sketches it at L3797 from
  `relationwithAnormalizedtypeB` + `RTVtheorem(ii)`. I read both regions first-hand and
  corroborated the expansion against his **definition** by an independent mechanism (the engine
  never mentions a positive root): 60/60 at `N≤4` both types, a factor-reassembly control 60/60,
  and four deliberate mis-transcription plants all detected. The RTV stable-envelope root
  underneath is **not** reduced by this — same shared root as Rick's `D^∞` and his Theorem W det
  lemma, registered as its own node rather than as a sentence in a paragraph.
- **My named positive control was weaker than I advertised.** The pre-registration named
  `N=3, m=2, I={2,3}` and claimed all blocks populated; at `N=3` with `I^c={1}` two of the six
  root families are **empty**. The all-six case is `N=4, m=2, I={2,3}`; both are reported.
- **G2** (non-cancellation at fixed `χ`) is untouched and nothing here bears on it.
