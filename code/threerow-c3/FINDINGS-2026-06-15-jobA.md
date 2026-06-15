# FINDINGS — Job A (2026-06-15 code session): Gap-2 decomposition scout + wide sweep

**Verdict.** The decomposition that closes **Compensation Lemma B** (Gap 2 of the
`(a,b,3)` family, `a` odd, `b` even, the two-generator wall) is the **two-binomial
term-by-term split** `Q₃ = (a−1)(b−2)·H − 720·C(j,6)`, run through the
*heavy/tip* product expansion. The scout confirms it term-by-term, pins the exact
obstruction (the single-binomial Gap-1 argument provably **cannot** carry this
branch), and stretches the verification clean to **m ≤ 120 (3363 shapes, 0
failures)**. Gap 2 was proved in full this cycle by the parallel PROVE session
(`proofs/2026-06-15-compensation-lemma-B-proved.md`); this scout is the
independent computational confirmation of every link, cross-checked two ways.

Script: `threerow-c3/c3gap2_scout_consolidated.py` (self-contained; every `v₂`
computed both directly and via Kummer/Legendre — **0 mismatches**). The end-to-end
proof chain `threerow-c3/c3gap2_FULLCHAIN.py` also re-runs at **0 fails** on all 9
links (LemA, LemC, heart, NL1, (i), B, star2, star1, MAIN).

---

## 1. Three-way split of `U(j)` by `v₂(j)` and `j mod 4`

`U(j) = v₂C(a+4,j) + v₂C(b+3,j) + [v₂Q₃(j) − v₂Q₃(0)]` — call the summands
`cA, cB, dQ`. Tabulating the **minimum of each summand** over the family by
`(v₂(j), j mod 4)`, split by `(a mod 4, b mod 4)`:

**Live class `a≡1, b≡2 (mod 4)`** (the only class with a tie / full box):

| (v₂j, j%4) | min cA | min cB | min dQ | min U |
|-----------:|-------:|-------:|-------:|------:|
| (0,1)  j=5,13,… | 0 | 0 | −10 | **−6** |
| (0,3)  j=7,15,… | 1 | 1 | −12 | −4 |
| (1,2)  j=6,10,… | 1 | 1 | −12 | **−5** |
| (2,0)  j=4,12,… | 0 | 0 | −10 | **−5** |
| (3,0)  j=8,24,… | 0 | 0 | −10 | −4 |
| (4,0)  j=16,… | 0 | 0 | −9 | −5 |

The **per-`v₂(j)` minima of `U`** are `{0:−6, 1:−5, 2:−5, 3:−4, 4:−5, …}` — the two
distinct deficit depths (`−6` at `v₂(j)=0`, then a `−5/−4` plateau) are the
**signature of the box `⟨2,4⟩`**: the depth-`6` deficit is generator-`2` (carried by
`j=5`), the depth-`5` plateau is generator-`4` (carried by `j=7`).

**Structural reading (this is the scout's main message to PROVE):** *the entire
deficit lives in `dQ`* — the binomials `cA, cB` are non-negative everywhere and, in
the live class, contribute exactly the **+2 of compensation** (both `cA=cB=1` on the
odd-side rows `(0,3)` and `(1,2)`). So the proof must split `dQ` against the
binomial slack — which is precisely the **Kummer relocation**: `v₂C(a+4,j) +
v₂C(b+3,j)` absorbs `2s₂(j)` plus the boundary digit-carries, collapsing `Δ̃` onto
`j+3 + 2·dQ + 2·[carries]`. Confirmed: in the three off-classes `(1,0),(3,0),(3,2)`
the same `dQ` deficit appears but `min U = −1` (no box) — the binomials there pay it
back to within one.

## 2. Per-term `v₂` floors of the `C(a+4,j)·Q₃` / `C(b+3,j)·Q₃` expansion

Split `Q₃ = (a−1)(b−2)H − 720C(j,6)` and form the two products
`P₁ = C(a+4,j)C(b+3,j)·H` (**heavy**) and `P₂ = 720·C(a+4,j)C(b+3,j)·C(j,6)`
(**tip**), so the object is `P = (a−1)(b−2)P₁ − P₂` with `U = v₂(P) − v₂Q₃(0)`. The
candidate per-term floors (the proof's `star1`, `star2`):

> `2 v₂(P₁) ≥ 2 v₂H(0) + (2s₂(j) − (j+3))`  (**heavy, star1**)
> `2 v₂(P₂) ≥ 2 v₂Q₃(0) + (2s₂(j) − (j+3))`  (**tip, star2**)

**Both hold, and both are TIGHT (min slack = 0)** across `m ≤ 80`. So a clean
term-by-term split with uniform per-term floors **does exist** — the two-generator
inequality is *not* forced off the term-by-term track. The mechanism: the tip
absorbs `C(j,6)` on the `b`-side via `C(b+3,j)·720C(j,6) = 720·(b+3)⁽⁶⁾·C(b−3,j−6)`
(reduces to NL₁ + Lemma B / `heart`), and the heavy product uses the rewrite
`H = (a+3)(b+2)G − 6E`, `E = C(j,2)Φ` with `Φ` odd, routed through `a`-absorption +
inequality `(i)`.

## 3. Where the single-binomial Gap-1 argument BREAKS (the obstruction)

The `a`-even Gap-1 branch closed with a **single-binomial** bound `v₂(dQ) ≥ 1 −
v₂(j)`. In the `a`-odd branch this naive bound **fails 28402 times** (`m ≤ 80`) —
e.g. `(a,b,j)=(9,6,6)` gives `dQ = −6 ≪ 1 − v₂(6) = 0`. **This failure is the
two-generator coupling.** It is decisive: the proof must use *both* binomials
simultaneously (each pays one of the two deficits) — a single binomial cannot, by
itself, cover the depth-`6` deficit. The term-by-term split survives only because
`cA` *and* `cB` are both available to absorb `s₂(j)` (Kummer), and the heavy factor
`(a−1)(b−2)` supplies the extra `+2` (see §4).

## 4. Equality locus + heavy-factor toggle

- **Equality locus pinned exactly.** Tie set `{j : Δ̃(j)=0} = {5,7,9}` **iff**
  `a≡1, b≡2 (mod 4)` (box interior `b≥9`), and **empty** in every other `a`-odd
  class — verified `True`, 0 exceptions, `m ≤ 80`.
- **Heavy-factor toggle (the `+2` that funds both deficits).** Min `Δ̃` over the box
  `{5,7,9}`, by class:

  | (a%4, b%4) | min Δ̃ on {5,7,9} |
  |:----------:|:----------------:|
  | **(1,2)**  | **0**  (full box `|J*|=4`) |
  | (1,0) | 2 |
  | (3,0) | 2 |
  | (3,2) | 2 |

  Toggling the congruence away from `(1,2)` jumps `Δ̃` from `0` to `≥2` — the box
  collapses. The supplier is the heavy factor: `v₂((a−1)(b−2)) = 0` for **every**
  `a`-even shape (odd heavy factor, no box) versus `≥2` for **every** `a`-odd shape
  (`4 | (a−1)(b−2)`). Census `m ≤ 60`: `a`-even → `{v₂=0 : 810}`;
  `a`-odd → `{2:196, 3:209, 4:154, 5:101, …}`, i.e. always `≥2`. This `+2` at
  offset `3` is exactly what lets the `v₂(j)=1` (gen-2, `j=5`) **and** `v₂(j)=2`
  (gen-4, `j=7`) deficits both be paid, with their XOR `j=9` for free — the
  `j ↔ j+4` residue-mod-`π²` involution closes the box.

## 5. Item 4 — wide sweep `m ≤ 120`

`Δ̃(j) ≥ 0` for all `a`-odd `λ=(a,b,3)`, `4 ≤ j ≤ b`: **3363 shapes, 0 failures**
(was `m≤79` in the proof note; now hardened to `m≤120`). Min `Δ̃` by `v₂(j)`:
`{0:0, 1:1, 2:1, 3:1, 4:7, 5:21, 6:51}` — equality (`Δ̃=0`) occurs **only at
`v₂(j)=0`** (i.e. odd `j`, the `{5,7,9}` carriers), strictly positive and growing
elsewhere. Box-tie `{5,7,9}` recorded **756 times**, every one in class `(1,2)` and
every one the *full* set — the all-or-nothing lock confirmed at scale.

---

### Cross-checks
- `v₂` direct vs Kummer/Legendre: **0 mismatches** on all binomials used.
- `c3gap2_FULLCHAIN.py`: all 9 proof links **0 fails**.
- Per-term floors (`star1`, `star2`): both `≥0`, both tight (slack 0).
- Equality locus, tie count, heavy-factor census: all as the Theorem predicts.

**Hand-off to PROVE:** the decomposition is the two-binomial heavy/tip split with
Kummer relocation; the obstruction is that the single-binomial Gap-1 bound is
provably insufficient (28402 counterexamples) and *both* binomials must fire; the
`+2` funding both deficits is the `4|(a−1)(b−2)` heavy factor at offset `3`. All of
this is now embodied in the completed proof `2026-06-15-compensation-lemma-B-proved.md`.
