# `threerow_c4_boundary` is machine-checked — c=1,2,3,4 boundaries now all `sorry`-free

**Date:** 2026-06-19 (lean session)
**File:** `TworowD4Kernel/ThreeRowC4Boundary.lean` (snapshot:
`projects/proofs/2026-06-19-lean-c4-boundary/`)
**Status:** `lake build` 0 errors / 0 warnings; `#print axioms threerow_c4_boundary` =
`[propext, Classical.choice, Quot.sound]` (standard only). 493 lines.

## What landed

The three-row `c = 4` boundary lemma `λ = (a,b,4)`, end-to-end. For the box-interior regime
`a = 2P+b+4`, `b ≥ 6`, all four boundary indices `j ∈ {b+1,b+2,b+3,b+4}` strictly exceed the
threshold `V = val(0)` (`θ = 0` for `c` even, both parities):

> `V < valb1 ∧ V < valb2 ∧ V < valb3 ∧ V < valb4`.

This joins `threerow_c1/c2/c3_boundary`. **All four three-row `d=4` boundary families are now
machine-checked.** The discharge is the `c = 2` (also even) template lifted to four indices:

- **`b` odd** (`v₂(a−2)=0`): subadditivity + Lemma A (`two_mul_sum_digits_le`) + content
  `v₂(N₂)≥2`, `v₂(N₁)≥2`. No `lemma_F` needed.
- **`b` even** (`v₂(a−2)=1+v₂(P+b/2+1)`): bridge `vz_prod_Icc` + `lemma_F` (indices i=4,2) or
  plain divisibility (i=3,1) + content `v₂(N₃)≥1`, `v₂(N₂)≥2`, `v₂(N₁)≥3`.

As the brief predicted: **`c=4` is even, so `a−b+1 = 2P+5` is odd and `lemma_F2` is NOT needed**
(confirmed). The lone `v₂(a−2)` deficit is absorbed by the factor `P+b/2+1` sitting inside the
consecutive-product run. Content bounds use the explicit even decompositions (`N₂ = 2²·Q`,
`N₃ = 2·Q` (b even), `N₁ = 2³·Q` (b even) / `2²·Q` (b odd)) — **no factorisation of `N_i`**, so the
`c ≥ 4` irreducibility wall never appears.

## What rests on MN verification (unchanged scope)

The closed-form `Δ`-hypotheses `hΔ4..hΔ1` and the `N_i` definitions are taken as given — they are
the MN-verified §1/§2 facts about the symmetric-function object `G_λ`, out of Mathlib's reach. Same
scoping as c=1,2,3. Everything 2-adic downstream (including the `N_i` content) is proved from scratch.

## ⚠️ A verification wrinkle worth your eyes (Infoview not needed — a maths/tooling note)

I tried to *independently* re-verify the `c=4` closed forms this session by recomputing `G_j` from
scratch via Murnaghan–Nakayama (my MN character engine passes all orthogonality / sum-of-squares
checks — it is correct). **I could not reproduce the closed forms** — and, alarmingly, I could not
reproduce even the *already-accepted* `c=3` `M_j` values (e.g. the c=3 note's `M_{b+2}=555` at
`(a,b)=(14,7)`) with any of the natural candidates I tried:
`⟨s_λ, e₂^j h₁^{2(m−j)}⟩`, `C(m,j)⟨s_λ, p₂^{m−j} e₂^j⟩`, `⟨s_λ, h₂^{m−j} e₂^j⟩`.

So my naive symmetric-function bookkeeping for the *physical* `G_j` is wrong — I don't currently have
the correct expression for the Baxterised coefficient (no Sage in-container this session, only my
hand-rolled power-sum + MN). This means **I could not give an independent MN cross-check of the
hypotheses** — exactly the same epistemic gap the `c=3` lean session had (it too took the paper's
MN-verified forms on trust).

What I **could** verify, and did (pure 2-adic number theory, 0 mismatches over 3400 `(P,b,i)`):
the digit-sum form I used for `hΔ4..hΔ1` is *exactly* `v₂` of the paper's §1 explicit-product ratio
`R_i` (`2026-06-16-c4-boundary-complete.md`, §1). So the Lean hypotheses faithfully encode the
paper's §1 ratio claim; the only unverified-by-me link is "§1 ratio `R_i` = true `G_{b+i}/G_0`",
which the prove session asserts MN-verified `m ≤ 41`.

**Recommendation:** when you next have Sage/the prove-session `mn.py` handy, it would be worth a
5-minute confirmation that the `c=4` `Δ` digit-sum forms in this file match a real MN run — to retire
the lingering doubt. The number theory is airtight; the symmetric-function input is trusted, not
re-derived. (If the §1 ratio is right, the file is correct; if it's off, the file faithfully proves a
statement about the wrong object — so the cross-check is genuinely the thing to do.)

## Next

`threerow_c5_boundary` — c=5 is **odd** (θ can be 3, `lemma_F2` returns, plus `(a−b+1)|N_2`
divisibility from `2026-06-17-generalc-master-and-c5.md`). Not started; flagged as the natural next
LEAN rung.
