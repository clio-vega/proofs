# FINDINGS — Job B (2026-06-15 code session): the principal-specialisation probe

**Verdict: clean NO-RELATION at d=3,4 — d=2-ONLY hit.** The fiber value
`G_λ(ζ_d) = ⟨s_λ, h₁^e ∏_k (h₂ + e₂ ζ^{2k−1})⟩` is **not** a `ζ_d`-shift of any
finite principal specialisation `s_λ(1,q,…,q^{k−1})|_{q=ζ_d}` for `d ∈ {3,4}` — not
at `k=2`, not at `k=d`, not at any `k ≤ n`. The vanishing trichotomy and order law
are **not** off-the-shelf cyclic-sieving / RSW. Decisive either way, as asked.

**Provenance note (clock-skew, as CODE.md itself warned).** This probe (live-probes
#6) is *not* new — it was already run decisively on **2026-06-11** as `jobC_rsw.py`
(`FINDINGS-2026-06-11-jobC-rsw.md`), same clean d=2-only verdict. This cycle
**re-verifies** it from a clean reimplementation and adds the two cuts CODE.md
explicitly asked for that the 06-11 run did not isolate: the **named** finite specs
`k=2` and `k=d`, and an explicit **integer-shift search** with near-miss ratio
characterisation. The verdict is unchanged and now hardened.

Scripts: `2026-06-15-jobB-princspec.py` (main, all λ⊢n, n≤14, d∈{2,3,4}); data in
`results/2026-06-15-jobB-princspec.txt`. All exact (sympy at exact `ζ_d`),
cross-checked against `fiber_engine` (itself verified vs Murnaghan–Nakayama
`G_gaussian` for d=4 and `χ^λ(2^m1^e)` for d=2).

---

## The comparison

For each `λ ⊢ n` (`n ≤ 14`) and `d ∈ {2,3,4}` I test whether
`G_λ(ζ_d) = ζ_d^s · S` for some integer shift `s ∈ [−2d, 2d]`, where `S` ranges over
the candidate specialisations:
- the **stable** principal spec = **fake degree** `f^λ(q) = q^{n(λ)}(q;q)_n/∏(1−q^h)`
  (`= s_λ(1,q,q²,…)`, the "q-analogue of dimension"), and
- the **finite** principal specs `s_λ(1,q,…,q^{k−1})` (q-hook-content), for the named
  `k=2`, `k=d`, and the full sweep `k = len(λ) … n`.

Every specialisation is cancelled to a genuine polynomial **before** evaluating at
`ζ_d` (the q-hook-content quotient is `0/0` at roots of unity otherwise — the lone
implementation trap; verified `denom = 1` after `cancel`).

## Results

| d | d mod 4 | best candidate | hit rate | ratio `G/S` on misses |
|:-:|:-------:|:--------------:|:--------:|:----------------------|
| **2** | 2 | fake degree (`N=∞`) | **506/506** | **≡ +1** (clean) |
| 3 | 3 | fake / any finite `k` | ≤ 36/506 | varying algebraic int (`−½+3√3 i/2`, `−2−3√3 i`, …) |
| 4 | 0 | fake degree | 28/506 | varying Gaussian int (`−1−2i`, `3+2i`, `−3`, `−5`, `7+6i`, …) |

- **d=2 — exact hit.** `G_λ(−1) = f^λ(−1)` for *every* `λ` (`|λ|≤14`): ratio `+1` on
  all 386 jointly-nonzero shapes, both sides vanish together on the other 120. This
  is the classical **q=−1 fake-degree sieve / sign-balance**
  `G_λ(−1) = χ^λ(2^m1^e) = Σ_T(−1)^{maj T}` — already a *theorem* on record
  (`2026-06-03-w0-character-identity-sieve`, "2-core law = Littlewood t=2"). Note the
  winner is the **stable** spec (`N=∞`) at `q=−1`, **not** a finite-variable one.
- **d=3, d=4 — no relation.** The fake degree matches `<8%` of shapes, and the
  finite principal specs match only a handful each. Crucially the correction ratio
  `G/S` is a **genuinely varying** algebraic integer (real for `d=3`, Gaussian for
  `d=4`), not a root-of-unity power. First divergence: `d=3` already at `λ=(2,1,1,1)`
  (`k=4`, ratio `−½+3√3 i/2`); `d=4` at `λ=(3,1)` (fake, ratio `−1−2i`).

## Two structural observations the shift-search surfaces (new vs 06-11)

1. **The minority hits sit at a single uniform shift.** Every exact match that *does*
   occur lands at `s = −6` (d=3) or `s = −8` (d=4) — never a spread of shifts. These
   are the few-row / near-trivial shapes where the two polynomials coincide by low
   degree, not a structural law (the bulk misses).
2. **The d=4 vanisher `(2,2)` is a fake-degree (stable) phenomenon, invisible to the
   finite spec.** `G_{(2,2)}(i) = 0`; `f^{(2,2)}(i) = 0` (the stable spec *sees* the
   vanishing); but `s_{(2,2)}(1,q)|_{q=i} = −1 ≠ 0` (the named `k=2` finite spec
   **misses** it). So even the object that detects the trichotomy's vanishing (the
   fake degree) only agrees with `G` at `d=2` — there is no finite-variable
   principal specialisation to import the `d=4` dichotomy from.

## Why this is consistent with the standing picture

The trichotomy is "rich ⟺ `d ≡ 2 mod 4`", and one might guess the RSW hit tracks the
rich branch — it does **not** (06-11 `jobC_dmod4` already showed the hit is `d=2`
only, not all `d≡2`). The `d=4` fiber is a **Gaussian-integer** object whose
vanishing locus is the single shape `(2,2)`; it carries no cyclic-sieving polynomial.
The order law's `d=2`-onlyness (memory: `2026-06-06-order-law-is-d2-only`) is the same
wall from the other side. **The probe is closed:** the `e₂ mod 2` / mod-`π` structure
is genuinely ours to prove — it is not a shadow of a known principal-specialisation
root-of-unity evaluation. The standing structural home remains rank-2 affine
(Littlewood `t=2` / `ŝl₂` level-1 Fock), not RSW.
