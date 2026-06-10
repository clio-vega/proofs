# FINDINGS — Job A: a 2-adic handle on `v₂(M_j)` (arms PROVE)

**Date:** 2026-06-10 (code session)
**Engine:** `job1_tie_census.py` (MN → `M_j`), `core_quotient.py` (new), `jobA_v2Mj.py`, `jobA_probe.py`, `jobA_probe2.py`.
**Recall:** `G_λ(i)=Σ_j C(m,j) iʲ(1+i)ʲ M_j`, `M_j=⟨s_λ,h₁^{2(m−j)}e₂ʲ⟩∈ℤ_{≥0}`, `M₀=f^λ`;
`val(j)=j+2v₂(C(m,j)M_j)`, `J*=argmin val`, `μ=min val`.

## Headline

PROVE wanted a clean Kummer/Lucas closed form for `v₂(M_j)`. **There isn't one** — and the data
says clearly *why not*. But the object PROVE actually needs — the **step law on J\*** — is exact,
and its two pieces split cleanly into a **provable Kummer binomial part** and a **forced M-part**.
That split, not a global `v₂(M_j)` formula, is the lever.

## 1. The table

`results/v2Mj-table.csv` — one row per `(λ,j)` with `M_j≠0`, all `λ⊢2m`, `m≤12` (45 566 rows).
Columns: `m, λ, j, M_j, v2Mj, v2binom, val, f, v2f, jbits, mbits, in_Jstar, core4, quot4, quotsize`.

## 2. No clean global closed form for `v₂(M_j)` (negative, but sharp)

Tested on all 45 566 data points / 4 114 shapes (`m≤12`):

| Candidate | Fit |
|---|---|
| `v₂(M_j) ≥ v₂(f^λ)` (monotone floor) | 17 413/45 566 (≈38%) — **FALSE** |
| per-shape Kummer shift `v₂(M_j)=carries(j,g)` | 92/4 114 shapes — **FALSE** |
| additive over set-bits of `j` (base `v₂(M_0)`) | 156/4 114 shapes — **FALSE** |
| `{j: M_j odd}` = submasks of one mask | 65/1 136 shapes — **FALSE** |

So `v₂(M_j)` is neither monotone in `j`, nor a Kummer carry-count, nor bit-additive, nor a clean
Lucas parity. A global formula is the wrong target.

### Why: the 4-core matters; the 4-quotient is not enough

Determination tests (a key = single-valued ⟺ those invariants pin `v₂(M_j)`):

| Determining set | single-valued keys |
|---|---|
| `(m,j)` | 3/65 — hopeless |
| `(m,j, 4-quotient)` | 6 317/8 460 (**≈75%**) — **insufficient** |
| `(m,j, 4-core)` | 400/819 — insufficient |
| `(m,j, v₂(f))` | 236/580 — insufficient |
| `(m,j, 2-quotient)` | 14 108/14 108 — *vacuous*¹ |
| `(m,j, 4-core, 4-quotient)` | 14 108/14 108 — *vacuous*¹ |

¹ Vacuous because core+quotient ⟺ λ (the abacus bijection); and 2-cores are staircases (unique per
triangular size), so `(m, 2-quotient) ⟹ |core| ⟹` the unique 2-core `⟹ λ`. These rows only confirm
the bijection, not a reduction.

**The non-vacuous signal:** the **4-quotient alone determines `v₂(M_j)` only 75% of the time** —
the 4-core carries the rest. This is exactly consistent with the spectral memory
`v_π(G)=v₂f^λ + base(4-core) + residual`: the valuation genuinely lives on the 4-core, and **does
not factor through the 4-quotient**. Anyone hunting a `v₂(M_j)` formula must keep the 4-core in.

## 3. The lever that DOES work: the step-law split (exact to m=14)

On every fixed-point-free involution pair `{j, j+2ᵃ}` of every tie's `J*` (a∈S), with
`Δv₂(M)=v₂(M_{j+2ᵃ})−v₂(M_j)` and `Δv₂(bin)=v₂(C(m,j+2ᵃ))−v₂(C(m,j))`:

```
val constant on J*  ⟺  2ᵃ + 2·Δv₂(bin) + 2·Δv₂(M) = 0
                    ⟺  Δv₂(M) = −2^{a−1} − Δv₂(bin).        (★)
```

**(★) holds 6685/6685 pairs (m≤14), 2398/2398 (m≤12).** The two pieces:

- **Binomial part** `Δv₂(bin)` is a **Kummer carry difference** — `v₂(C(m,j))=` #carries of
  `j+(m−j)` base 2 (Kummer's theorem), a closed form. *Provable outright.*
- **M-part** `Δv₂(M)` is then **forced** by (★). It is *not* independently clean: its distribution
  spreads (`Δv₂(M)∈{−3,−2,−1,0,1}` for a=1) precisely tracking `−Δv₂(bin)`.

**Clean sub-structure (the gift):**
- **1010/2398 (m≤12), 3744/6685 (m≤14) pairs have `Δv₂(bin)=0`** — the binomial valuation is flat
  across the pair, and then **`Δv₂(M)=−2^{a−1}` exactly** (3744/3744). These are the "M carries
  everything" toggles.
- **Generator `a=3` is ALWAYS a pure-M toggle:** every `a=3` pair has `Δv₂(bin)=0` and
  `Δv₂(M)=−4=−2²` (168/168 at m≤14). The big generators are pure M-jumps.

> **Sign note for PROVE.md:** CODE.md wrote `v₂(M_{j+2ᵃ})−v₂(M_j)=2^{a−1}−Δv₂(bin)`; the verified
> sign is **`−2^{a−1}`** (eq. ★) — `2ᵃ>0` must be *cancelled*, so the M-part drops. Fix the sign.

## 4. The one clean closed form found: `v₂(f^λ)` via the 2-quotient

For every `λ` with empty 2-core (`m≤10`, 1 214 shapes, **0 cex**):

```
v₂(f^λ) = v₂( multinomial(2m; 2|q₀|, 2|q₁|) ) + v₂(f^{q₀}) + v₂(f^{q₁}),
```

where `(q₀,q₁)` is the 2-quotient. This is the classical James–Kerber/Macdonald `d`-quotient
dimension factorisation. It is a **provable anchor for `M₀=f^λ`** — the base of every `val(j)`
ladder — and the natural seed for an induction up the `j`-ladder via (★).

## What this hands PROVE

1. **Don't chase a global `v₂(M_j)` formula** — it provably has none (§2); the 4-core obstructs any
   quotient-only form.
2. **Target the step law (★) directly.** Its binomial half is Kummer (closed form, provable). The
   theorem reduces to the M-half on `J*` pairs: `Δv₂(M)=−2^{a−1}−Δv₂(bin)`.
3. **Easiest entry: the pure-M toggles** (`Δv₂(bin)=0`, all `a=3`, and 56% of all pairs):
   prove `v₂(M_{j+2ᵃ})−v₂(M_j)=−2^{a−1}` when the binomial valuation is flat across the pair.
4. **Anchor `M₀` with §4** and climb the `j`-ladder.
5. `S⊆{1,2,3}` (max generator = 3) **re-confirmed to m=14** (no generator ≥ 16; the e₂-degree
   ceiling is never approached). The PROVE target keeps its 3-generator shape.
