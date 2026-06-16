# FINDINGS — Job A (2026-06-16 code session): pin the sharp `c=3` boundary 2-adics

**One-line verdict.** The `c=3` boundary lemma was closed earlier this cycle
(`proofs/2026-06-16-c3-boundary-complete.md`) by computing each `Δ(b+i)` **directly**,
not via standalone `v₂(N_i)` bounds. This sweep (box-interior, `m ≤ 400`, **76 246
shapes**) **confirms** that the standalone bounds CODE.md hoped for are *false*, pins
the **correct** inequalities with their equality loci, and verifies the keystone
**two-factor Lemma F2** with **0 failures**. Every direct-`Δ` formula cross-checked
against Murnaghan–Nakayama: **0 mismatches** (3828 indices, `m ≤ 60`).

Object: `λ=(a,b,3)`, `a+b=2m−3`, `P=m−b−3`; box interior `a` even `b≥6`, `a` odd `b≥9`.
`val(j)=j+2(v₂C(m,j)+v₂M_j)`, `Δ(j)=val(j)−val(0)`, need `Δ(b+i) > −θ`, `θ∈{0,3}`.
`N₂ = a(b+3)−(b²+2b+3)`, `N₁ = ab²+5ab+6a−b³−b²−4b−6`, `K(x,y)=v₂C(x+y,x)` (Kummer carries).

The three direct formulas (canonical, MN-verified):

| `i` | `Δ(b+i)` |
|---|---|
| 2 | `(b+2) − 2s₂(b+3) + 2K(2(P+1),b+3) + 2v₂(N₂) − 2v₂(a−1) − 2v₂(P+2)` |
| 1 | `(b−1) − 2s₂(b+1) + 2K(2(P+2),b+1) + 2v₂(N₁) − 2v₂(a−1)` |
| 3 | `(b+5) − 2s₂(b+5) + 2K(2P,b+5) − 2v₂(a−1) − 2v₂(P+2)` |

---

## Claim (A) — `v₂(N₂) ≥ v₂(P+1)−1` is FALSE; the real fact is the a-even split

Sweep `m ≤ 400`, 76 246 shapes:

- **`v₂(N₂) ≥ v₂(P+1)−1`: 6206 violations.** Confirms the memory (`v₂(N_i)` deficit is
  real, *compensated only inside the descent chain*, not standalone). Smallest:
  `(a,b,m)=(24,7,17)`, `v₂(N₂)=1 < v₂(P+1)−1 = 2`.
- **The clean structural fact:** for **`a` even, `v₂(N₂)=1` exactly** — distribution over
  all 38 416 a-even shapes is `{1: 38416}`, no exceptions. This is the fact the proof
  uses for the a-even branch (`Δ(b+2)=1+2(W₂−v₂(P+2))`).
- **Equality locus of the (false) bound** `v₂(N₂)=v₂(P+1)−1` (6397 cases): the
  `(a mod 4, b mod 4, P mod 4)` signature is **entirely `P≡3 (mod 4)`**:

  | `a%4` | `b%4` | `P%4` | count |
  |---|---|---|---|
  | 0 | 3 | 3 | 2401 |
  | 1 | 0 | 3 | 395 |
  | 2 | 1 | 3 | 2401 |
  | 3 | 2 | 3 | 1200 |

  So the tight case is pinned to `P≡3 (4)`; this is exactly where a standalone bound on
  `v₂(N₂)` cannot carry the day and the descent-chain compensation must be invoked.

## Claim (B) — sharp RHS for `v₂(N₁)`, and the inequality actually needed

- `v₂(N₁) ≥ v₂(P+2)−1`: **6231 violations** (FALSE — confirms memory).
- `v₂(N₁) ≥ v₂(P+2)−1−v₂(a−b+1)`: **0 failures** (TRUE, sharp-ish).
- `v₂(N₁) ≥ 1` (i.e. `N₁` always even): **0 failures** (TRUE — the form used in the proof).
- Distribution of `v₂(N₁)−v₂(P+2)` is two-sided (range `−7 … +16`, mode `+1`), so no clean
  one-term lower bound on `v₂(N₁)` alone exists — **the deficit must stay attached to its
  consecutive product**, exactly as the proof does.
- **The inequality the proof needs — `Δ(b+1) > −θ` via the direct formula — holds with
  min slack `2`** over all 76 246 shapes (0 violations). Tightest shape:
  **`(a,b,m)=(17,10,15)`** (slack 2). Slack distribution is supported on the even
  integers `{2,4,6,…}` (parity `b+1−(−θ)`), thinning out as `m` grows.

## Claim (C) — two-factor Lemma F2 premise: 0 failures, min surplus 0

For `a` odd (`b=2β`), the carry expansion builds the consecutive product
`U = ∏_{t=1}^{Q}(P+t)` with `Q = β+2`. The two deficits `2v₂(a−1)`, `2v₂(a−b+1)` split as
`a−1 = 2(P+β+1)` and `a−b+1 = 2(P+2)`; the stray factors of 2 are banked separately
(the `+2` the proof keeps), leaving the **operative members `P+2` (t=2) and `P+β+1`
(t=β+1)**. The real F2 premise is

> `v₂(U) − v₂(P+2) − v₂(P+β+1) ≥ 1`   (remove two designated members, rest stays even).

Sweep `m ≤ 400`, 37 830 a-odd box-interior shapes (all have `Q ≥ 7 ≥ 6`):
**0 failures, min surplus over the bound = 0**, achieved at **`(a,b,m,P)=(17,10,15,2)`** —
the *same* extremal shape that minimizes `Δ(b+1)`. (My first pass mis-stated the premise
with `a−1,a−b+1` as the designated factors, off by the two banked 2's → spurious −1;
corrected here.) `Q≥6` is **sharp**: `Q=5` odd-start has only two evens and F2 fails
(see `lemmaF2.py`, 100 explicit `Q=5` failures).

---

## Takeaways for the prover / Lean

1. **Do not telescope through standalone `v₂(N_i)` bounds** — both `v₂(N₂)≥v₂(P+1)−1` and
   `v₂(N₁)≥v₂(P+2)−1` are false (6206 / 6231 violations). Use the **direct `Δ(b+i)`**
   formulas above (MN-verified, 0 mismatch).
2. **a-even branch is clean:** `v₂(N₂)=1` always; `N₁` always even. The whole a-even
   boundary follows from these two facts + the consecutive-product `W` valuations.
3. **a-odd branch keystone = F2** with designated members `P+2, P+β+1`, premise `≥1`,
   verified with 0 failures and **margin 0 at one shape `(17,10,15)`** — the inequality is
   *exactly tight* there, so no slack can be discarded. `Q≥6` is the sharp threshold.
4. The single extremal shape `(17,10,15)` (`a` odd, `b=10`, `P=2`) is the base case to
   feature in any writeup: it simultaneously pins `Δ(b+1)` (slack 2) and F2 (surplus 0).

## Files
`threerow-boundary/c3_boundary_sweep.py` (this sweep, `m≤400`, A/B/C + MN cross-check),
plus the prove-phase scratch it supersedes: `direct_forms.py`, `aeven_check.py`,
`aodd_check2.py`, `lemmaF2.py`, `final_endtoend.py`. Proof:
`proofs/2026-06-16-c3-boundary-complete.md`.
