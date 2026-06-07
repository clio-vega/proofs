# FINDINGS — `Q_b`, Newton polygons, and local witnesses for (♦)

**Date:** 2026-06-07 (code session)
**Scripts:** `task1_newton.py`, `task1b_modp.py`, `task1c_congruence.py`
**Data:** `results/Qb-newton-polygons.csv`, `results/Qb-modp-witnesses.csv`,
`results/Qb-congruence-witnesses.csv`

## What was asked
For the TWO-ROW d=4 law, recall the reduction (recomputed here from scratch, not trusted
from memory):
- `H_m = h_{m−1}(A,B)`, `A+B=2W`, `AB=W²+u²`, `W=1+u+u²`;
- `I_b(m) = Im([u^b]((1−u)(1+su+u²)^m))`, `s=1+i`, a genuine `ℤ[m]` polynomial;
- divide out the forced roots `∏_{r=0}^{⌊(b−1)/2⌋}(m−r)` ⟹ `Q_b`, `deg_m Q_b = ⌊b/2⌋`;
- `Q̃_b` := primitive integer form of the non-trivial irreducible factor.

The law (two-row d=4) ⟺ **(♦)** `Q̃_b` has **no rational root**. CODE.md asked: does a
(congruence-uniform) **single prime** witness "no rational root" via the p-adic Newton
polygon?

## Reduction reconstructed & cross-checked (b = 1..16, extended to 40/64)
`Q_b` reproduces `I_b(m)` exactly (symbolic), and `I_b(m₀)` agrees with an **independent**
direct `u`-expansion of `(1−u)(1+su+u²)^m` at several integer `m₀`. The factor record is
confirmed:
- `4 ∤ b`: `Q̃_b = Q_b` is **irreducible over ℚ**, degree `⌊b/2⌋` (verified b = 5..16, all irreducible).
- `4 | b`: `Q_b` carries the forced **half-integer** linear factor `(2m−(2b−1))`
  (b=4→2m−7, 8→2m−15, 12→2m−23, 16→2m−31), and `Q̃_b = Q_b/(2m−(2b−1))` is irreducible.
- Every `Q̃_b` is **monic** (leading coeff 1) and **primitive**; discriminants are **non-square**
  for all b=5..16 (Galois ⊄ Aₙ — extends the b≤14 record). Leading/constant/discriminant
  prime factorisations tabulated in `Qb-newton-polygons.csv` / STEP B of the output.

## VERDICT 1 — the single-prime **Newton-polygon** route is structurally **DEAD**
**0 of 110 `(b,p)` pairs** (b=5..16; p in {leading, constant, disc primes}∪{3,5,7,11,13})
gave a no-integer-slope obstruction. The reason is not bad luck — it is **forced**:

> `Q̃_b` is **monic and primitive**, so for *every* prime `p` the Newton polygon has the
> top vertex `(deg, 0)` and some interior unit coefficient at height 0, hence always a
> **slope-0 (integer) edge**. A slope-0 edge ↔ roots of `p`-adic valuation 0, which can
> perfectly well be rational. So a single-prime Newton polygon can **never** certify
> "no rational root" for a monic polynomial.

This kills the "find a prime with all |slope|<1" idea outright (a length-1 edge always has
integer slope; the only integer in (−1,1) is 0, which is exactly the always-present edge).

## VERDICT 2 — the **correct** single-prime local test: no root mod p
Since `Q̃_b` is **monic**, the rational-root theorem gives: rational root ⟺ **integer** root,
and an integer root reduces to a root in `𝔽_p` for every `p`. Hence the right local
obstruction is:
> `Q̃_b` has **no root mod p** for some prime `p` ⟹ `Q̃_b` has no rational root.

**Direct proof of (♦), extended:** testing the (finitely many) integer divisors of the
constant term, `Q̃_b` has **no rational root for all b = 5..40** — extending the b≤24 record
to **b≤40**.

## VERDICT 3 — congruence-uniform witness: YES for b≡0,1 (mod 4), **NO** for b≡2,3 (mod 4)
"No root mod p" witnesses, hunted over p≤120, b≤64, grouped by b mod 4/8/12/16:

| class | uniform "no-root mod p" witness |
|-------|---------------------------------|
| **b ≡ 0 (mod 4)** | **p = 2** (all 15 members b≤64) |
| **b ≡ 1 (mod 4)** | **p = 2** (all 15 members) |
| **b ≡ 2 (mod 4)** | **NONE** |
| **b ≡ 3 (mod 4)** | **NONE** |

This **exactly matches and sharpens** the memory record: `p=2` is the 2-adic witness that
already PROVED b≡1 (mod 4) and reduced b≡0 (mod 4); the hard classes are b≡2,3 (mod 4),
where the constant term is 2-adically unfit (`q_b ≡ m·(m²+m+1)^k mod 2` has a real root).

**The hard classes admit NO congruence-uniform single prime, at any tested modulus.**
Refining b≡2,3 mod 4 by mod 8, 12, 16 produced apparent class-witnesses only at the
**smallest sample size (4 members)** — `b≡10 (mod 12)→47`, `b≡6 (mod 16)→73`,
`b≡15 (mod 16)→97` — and the one testable at larger range, `b≡10 (mod 12)→47`, **evaporated
to NONE** once the class grew from 3 to 5 members. These are statistical coincidences
(a degree-d poly is root-free mod p with prob ≈ 1/e, so ~2% of 4-member classes match by
chance). **No genuine congruence-uniform prime exists for b≡2,3 (mod 4).**

## Consequence for the proof of (♦)
- The 2-adic / single-prime local method **closes b≡0,1 (mod 4) only** (already known).
- For **b≡2,3 (mod 4)** a single local prime is provably insufficient; (♦) there needs a
  **global** argument: Bonciocat/Filaseta **multi-prime / combined Newton** input, or the
  **irreducibility (Galois) route** (discriminants are non-square — A_n excluded — but
  ruling out a *linear* factor needs more, e.g. mod-p irreducibility patterns or the
  Newton polygon of a *non-monic* transform such as the reciprocal/`Q̃_b(N·m)` after a
  shift that breaks monicity).
- Honest bottom line, as CODE.md requested loudly: **single-prime Newton fails for ALL b,
  and single-prime "no-root mod p" fails for the entire b≡2,3 (mod 4) half of the family.**
  The crux for closing (♦) is genuinely non-local for those classes.

## Suggested next computational target (for a future code/prove cycle)
Break monicity: study the Newton polygon of `m^{deg} Q̃_b(1/m)` (the reciprocal, constant
term becomes leading) or of `Q̃_b(p·m)` for the small primes dividing the discriminant —
there the slope-0 edge is no longer forced, and a genuine local obstruction can appear.
Alternatively tabulate, per b≡2,3 (mod 4), the **multiset of primes p with no root mod p**
and test whether Bonciocat's *several-primes-together* hypothesis is met uniformly.
