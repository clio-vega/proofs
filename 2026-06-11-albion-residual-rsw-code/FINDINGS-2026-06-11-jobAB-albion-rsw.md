# FINDINGS — 2026-06-11 code session: Albion z-asymmetric residual (Job A) + graded RSW (Job B)

**Scripts:** `jobA_albion.py`, `jobB_graded_rsw.py` (reuse `fiber_engine.py`, the proved master
identity, and the `psi_power_schur` engine).
**Data:** `results/jobA-albion-charge.csv`, `results/jobB-graded-rsw.csv`.
**Both probes were decisive. Both are clean NEGATIVES — and each one sharpens the gap.**

---

## Job A — Albion's z-asymmetric 4-core/4-quotient does NOT close the residual

### The object (Albion 2501.18520, definitions extracted from the paper)
The Littlewood decomposition `λ ↔ (t-core, t-quotient)` via the t-Maya diagram. A t-**core** `μ` is
encoded by its **charge vector**
`κ_t(μ) = (c_0,…,c_{t-1})`, `c_r = (#beads on runner r) − N/t`, `Σ c_r = 0`
(the runner imbalance). z-asymmetric partitions `λ=(u+z|u)` (Ayyer–Kumari's Frobenius form) are
exactly those whose `(κ, quotient)` obey the conjugacy relations `c_r + c_{z−r−1}=0` etc. (Thm 2.3).
Character factorisation (Thm 1.3): `φ_t s_λ = sgn_t(λ)·s_{λ^(0)}···s_{λ^(t-1)}` when `t-core(λ)=∅`,
where `φ_t` is the Verschiebung adjoint to plethysm by `p_t`.

**So the "z-asymmetric 4-quotient datum" = `(charge vector κ_4(4-core), ordered 4-quotient)`.**

### Result (all λ ⊢ 2m, m ≤ 9, n ≤ 18; 851 shapes with finite residual)

**Q0 — the datum is the COMPLETE invariant (tautology trap, now made explicit).**
`κ_4` is a *bijective* invariant of the 4-core (34 cores ↔ 34 charge vectors, verified). Hence
`(κ_4, ordered quotient) ↔ (4-core, ordered quotient) ↔ λ`. Testing "residual = f(κ_4, quotient)"
gives **0/851 non-constant** — but *only vacuously*, exactly as the 06-07 K4 complete-invariant
sanity row did. **This also corrects the 06-07 appendix framing:** its head-line
"residual = f(t_j(core), ordered quotient), 0/2375" was likewise the complete-invariant tautology
(`t_j(core)` canonical ⟺ `κ_4` ⟺ core), **not** a genuine reduction.

**Q1 — it separates `(1^6)` vs `(4,4,2)`, but only because their cores differ.**
`(1^6)`: core `(1,1)`, `κ=(1,0,−1,0)`, residual 0. `(4,4,2)`: core `(4,1,1)`, `κ=(0,−1,0,1)`,
residual 1. Different charge ⇒ separated — which the complete invariant already knew.

**Q2 — the decisive ladder: does the charge let us forget the quotient's fine structure? NO.**

| datum | #groups | non-constant |
|---|---|---|
| ordered quotient alone | 164 | 58 |
| charge alone | 33 | 24 |
| **(charge, quotient SIZES)** | 473 | **72** |
| (charge, quotient size multiset) | 113 | 16 |
| (charge, per-runner `v₂(f^{λ^(r)})`) | 57 | 36 |
| (charge, quotient PARTS per runner) = complete inv | 851 | **0** |

Residual is **not** a function of `(charge, sizes)`, nor `(charge, per-runner v₂f)`. It collapses to
0 non-constant **only** at the full quotient partitions — i.e. the complete invariant. The
z-asymmetric/charge bookkeeping is the *wrong kind of refinement*.

**The killer witness (trivial charge `(0,0,0,0)`, 3 boxes all on runner 2):**
`λ^(2)=(3,)` → residual **0**, `(1,1,1)` → **0**, but `(2,1)` → **1**.
Same charge, same runner, same size — residual is reading the *internal shape* of the quotient
partition. Note `f^{(2,1)}` has `v₂=1` while `f^{(3,)}, f^{(1,1,1)}` have `v₂=0`: the residual sees
**multiplicative 2-adic hook data inside the quotient**, not abacus/runner data.

**Q3 — residual is NOT additive over quotient boxes.** Fitting
`residual ≈ base(charge) + Σ_r w(charge,r)·|λ^(r)|` gives only 260/688 exact; the interaction
`Δ = residual − additive` ranges over `{−7,…,+3}` (distribution: `0:260, ±1:238, ±2:151, 3:26,
−7:4,…`). There is a genuine charge-mediated **cross-box interaction term** (e.g. two single boxes on
runners 2,3 at charge `(0,1,0,−1)` each give increment 0 but together give residual 1).

### Verdict A (clean negative, with the gap re-pointed)
Albion's z-asymmetric framework gives the right **name** for the core summary (the charge vector
`κ_4`) but supplies **no closed form** for the residual. The missing coordinate is *not* runner/charge
data — it is (i) the multiplicative `v₂`-of-hooks fine-structure **inside** each quotient partition,
plus (ii) a charge-mediated cross-box interaction. This is precisely what the **dimension-drop
verdict** predicts a real-valued abacus classification cannot see. The residual closed form, if it
exists, lives in a **hook-length / 2-adic formula for `f^λ` against its quotient pieces**, not in
abacus runner combinatorics. (Caveat from CODE.md honoured: verified on ALL n≤18, not a sample.)

---

## Job B — the graded `G_λ(q)=Σ_T q^{s(T)}` is NOT a q-shift of a principal specialisation
### …but the representation cross-check is a genuine WIN

`G_λ(q) = Σ_{T∈SYT(λ)} q^{s(T)}`, `s(T)=Σ_{i∈Des(T)} w_i`, `w_i = 2i−1` if `n−i` odd else `0`
(the parity-twisted branch-exponent statistic).

**STEP 0 (new, the valuable result): `Σ_T q^{s(T)}` equals the master-identity polynomial
`G_poly(λ)(t) = ⟨s_λ, h₁^e ∏_{k∈A}(h₂+e₂t^{2k−1})⟩` — exactly, 271/271 for all λ⊢n, n≤12, 0
mismatches.** The SYT descent-statistic generating function and the symmetric-function fiber
polynomial are **one and the same object**. This certifies that the (already-run) `jobC_rsw` verdict,
built via the master identity, transfers verbatim to the graded SYT object — the two threads
`{d_j}={s(T)}` and `G_λ(ζ_d)=⟨s_λ,ψ^m⟩` are literally the same polynomial.

**STEP 1–2 (principal specialisation `s_λ(1,q,…,q^{N-1})`, menu of N, d=2,3,4,5,6):**
no uniform q-shift. Best non-trivial is d=2 at **19/27** (still a miss); d≥3 only the trivial
N=1 one-row hits. Consistent with `jobC_rsw`'s richer numeric sweep (n≤14): the principal
specialisation does **not** import the trichotomy. (Note: d=2 *does* match the **fake degree**
`f^λ(q)` 507/507 — the classical sign-balance sieve — but that is a *different* polynomial from the
principal specialisation, and is d=2-only.)

**STEP 3 — the d=4 vanisher `(2,2)`, transparent at last:**
`G_(2,2)(q) = q⁶ + 1`. It vanishes at `q=ζ_4` because `i⁶ = −1`. The principal specialisation
`s_(2,2)(1,q,q²) = q⁶+q⁵+2q⁴+q³+q²` *also* vanishes at `ζ_4`, but it is a **different polynomial**
— the coincidence is not a principal-specialisation identity. Likewise the d=3 vanishers are
transparent: `G_(3,1)=q⁵+q+1` and `G_(2,1,1)=q⁶+q⁵+q` each vanish at `ζ_3` (`1+ζ_3+ζ_3²=0`).

### Verdict B
`G_λ(q)` is **not** a q-shift of any principal specialisation `s_λ(1,q,…,q^{k−1})` for d≥3 — the
fiber/grade is genuinely the seed's own object, not an off-the-shelf RSW/cyclic-sieving polynomial.
The escape hatch "maybe it's all principal-specialisation RSW" is **closed**. The lasting deliverable
is STEP 0: a certified identity `Σ_T q^{s(T)} ≡ ⟨s_λ, h₁^e ∏(h₂+e₂t^{2k−1})⟩`, n≤12.

---

## Bottom line for the programme
- **Job A:** Albion's z-asymmetric refinement is the wrong handle for the residual; the residual reads
  `v₂`-hook structure *inside* the quotient + a cross-box interaction, re-confirming the dimension-drop
  verdict and pointing the closed-form search at a 2-adic `f^λ`-vs-quotient hook formula.
- **Job B:** graded SYT `G_λ(q)` = master-identity fiber polynomial (certified n≤12); it is not a
  principal specialisation for d≥3; the d=3,4 vanishers are transparently cyclotomic
  (`q⁶+1` at ζ_4; `q⁵+q+1`, `q⁶+q⁵+q` at ζ_3).
