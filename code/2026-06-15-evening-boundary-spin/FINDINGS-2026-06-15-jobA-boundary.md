# FINDINGS — Job A (2026-06-15 evening): Boundary residual, structural scout

**One-line verdict for the prover.** The Boundary Lemma (Gap 3) is *true and strict
with margin ≥ 2*, verified for **c = 1, 2, 3 up to m ≤ 200** (≈ 59 000 shapes, **0
ties in-gate**) and confirmed for **c = 4, 5** (margin ≥ 8). The only boundary ties
anywhere are a **finite, explicit list** — `c=1: b∈{1,2,4}`, `c=2: b=2`, `c=3: b=6` —
which are *exactly* the per-family boundary-tie shapes the interior notes already
excise. **No tie contradicts any interior theorem.** And: the natural "v₂(M_{b+i})
is bounded" hypothesis in CODE.md Task 3 is **FALSE**; the real structural reason is
the proof's own factor-in-product mechanism (Lemma F), reconfirmed below.

Scripts: `threerow-boundary/jobA_widesweep.py` (c≤3, m≤200, pure-Python closed
forms), `threerow-boundary/jobA_c45.py` (c=4,5, exact Jacobi–Trudi determinant).
Every closed form cross-checked against Murnaghan–Nakayama (`mn.py`): **0 mismatches,
6741 cases**.

Object (from the proof, §0): `G_j = ⟨s_λ, e₂^j h₁^{2(m−j)}⟩ = C(m,j)·M_j`,
`val(j) = j + 2·v₂(G_j) = j + 2(v₂C(m,j) + v₂(M_j))`, `V = min_j val(j)`,
`J* = {j : val(j)=V}`. `M_j` supported on `0 ≤ j ≤ b+c`. Boundary indices = `{b+1,…,b+c}`.

---

## 1. Closed boundary forms `M_{b+i}` (Task 1) — confirmed c = 1…5

Fitted exactly (`boundary_forms.py`) and re-derived; here factored over ℤ:

| c | `M_{b+1}` | `M_{b+2}` | `M_{b+3}` | … | top `M_{b+c}` |
|---|---|---|---|---|---|
| 1 | `b` | — | — | | `b` |
| 2 | `(b−1)(a(b+2)−b(b+1))/2` | `(b−1)(b+2)/2` | — | | `(b−1)(b+2)/2` |
| 3 | `(b−2)(a−b+1)(ab²+5ab+6a−b³−b²−4b−6)/12` | `(b−2)(b+2)(a(b+3)−(b²+2b+3))/6` | `(b−2)(b+2)(b+3)/6` | | `(b−2)(b+2)(b+3)/6` |
| 4 | (deg-7 in a,b) | (deg) | `(b−3)(b+2)(b+3)(a(b+4)−(b²+3b+8))/24` | | `(b−3)(b+2)(b+3)(b+4)/24` |
| 5 | … | … | … | | `(b−4)(b+2)(b+3)(b+4)(b+5)/120` |

**The uniform-c pattern (the route to "all c").** The **top** term is `a`-FREE and is a
single hook:
> `M_{b+c} = (b−c+1)·∏_{t=2}^{c}(b+t) / c! = (b−c+1)(b+c)! / [(b+1)! c!]`  (Lemma T).

Each step down from the top *multiplies in one extra `a`-linear factor* `N_i`
(e.g. `c=3`: `N₂ = a(b+3)−(b²+2b+3)`, `N₁ = ab²+5ab+6a−b³−b²−4b−6`). Which `D_j(p,q,r)`
of the signed 3-row Jacobi–Trudi survive at `j = b+i`: at the **top** `j=b+c` only the
identity term survives (giving the `a`-free hook); each lower `i` revives one more
permutation, contributing the `a`-linear factor. This is the surviving-`D_j` structure
that makes Lemma T uniform in `c` and isolates exactly `c−1` `a`-dependent factors per
family. (Cross-checked: Lemma T + Δ(b+c) hold `c≤5, m≤40`, 1520+2995 cases, 0 bad.)

---

## 2. The valuation gap (Task 2) — wide sweep to m ≤ 200

For each shape compute `j₀` = smallest interior (`j≤b`) argmin of `val`, gate on
**`j₀ + 2c ≤ b`** (the lemma's hypothesis: the interior 2-adic box `{j₀,…,j₀+2c}`
fits inside `[0,b]`), and compare `min_{i} val(b+i)` to the interior `V`.

| c | shapes (m≤200) | **IN-GATE min Δval = val(bdry)−V** | in-gate ties | sub-gate (excluded) b-values |
|---|---|---|---|---|
| 1 | 19 900 | **2** (at (4,3,1)) | **0** | b ∈ {1,2,4} |
| 2 | 19 701 | **2** (at (5,5,2)) | **0** | b ∈ {2,3} |
| 3 | 19 306 | **2** (at (17,10,3)) | **0** | b ∈ {3,4,5,6,8} |

c=4 (m≤40): min Δ = **8**, 0 ties. c=5 (m≤34): min Δ = **8**, 0 ties.

So **above the gate the boundary always loses, by ≥ 2**, with no exception out to
m = 200 — a clean ×2.5 push past the m≤79–80 frontier.

## 3. Tie check (Task 4) — the only ties are a finite list

Scanning **all** shapes (no gate), a *genuine* boundary tie (`J*` contains an index
`>b`) occurs at exactly:

| c | b-values with a genuine boundary tie | max b | shape of `J*` |
|---|---|---|---|
| 1 | **{1, 2, 4}** | 4 | `{3}` or `{3,5}` (a-odd) / `{0,2}` (a-even, b=1) |
| 2 | **{2}** | 2 | `{0,4}` |
| 3 | **{6}** | 6 | `{3,5,7,9}` (a≡1, b≡2 mod 4) |

These match the proof's excluded families verbatim (§2 `b∉{1,2,4}`; §3 `b≥3`; §4
`c=3` box-interior `b≥6/9`). **The mechanism is transparent:** a tie is a *boundary*
tie only while `b` is small enough that part of the interior box `{j₀,…,j₀+2c}` pokes
out past `b`. E.g. `c=3, a` odd has `J*={3,5,7,9}`; at `b=6` the indices `7=b+1,
9=b+3` are boundary (→ tie), but at `b=10,14,…` (same `a≡1,b≡2 mod4`) the whole box
sits in `[0,b]` and the boundary strictly loses. So there are **no new ties** at large
`b` — the interior theorems need no correction.

---

## 4. The structural reason (Task 3) — CORRECTION + the right statement

**CODE.md's hypothesis "v₂(M_{b+i}) is bounded O(1)" is FALSE.** Measured maxima over
the sweep:

| c | max v₂(M_{b+i}), any i | where | max v₂(top M_{b+c}) | max v₂(M₀) |
|---|---|---|---|---|
| 1 | 7 | (a,b)=(129,128): M_{b+1}=b=128=2⁷ | 7 | 22 |
| 2 | 12 | (272,32) | 6 | 13 |
| 3 | 17 | (239,98) | 8 | 23 |

Even the `a`-free top `M_{b+c}` is unbounded (it contains `b`-factors: `b=2^k ⟹ v₂≥k`).
So one cannot pin the boundary below by bounding `v₂(M_{b+i})` alone.

**The reason the boundary loses is the proof's own engine, not boundedness.** Work with
`Δ(b+i) := val(b+i) − val(0)` (parity `≡ b+i`). The exact reduction (hook Lemma H +
Lemma T + Kummer) writes
> `Δ(b+i) = (positive linear in b) − 2·s₂(·) + 2·carries(·) − 2·v₂(N_i)`,
i.e. *positive-linear-in-b + carries − a deficit*, where the deficit `2v₂(N_i)` is `2·v₂`
of an `a`-linear factor. **Lemma F (factor-in-product)** closes it: that very factor
`P+t₀` is one term of the consecutive-integer product `∏_{t=1}^{Q}(P+t)` whose 2-adic
valuation the carry-sum *builds* (`carries + v₂(Q!) = v₂∏`, Kummer), and `Q≥4`
consecutive integers always carry a spare even factor. So the deficit is dominated *by
construction*. The boundary index `j=b+i` is large, so its `+(b+i)` linear term plus the
binomial carries swamp the lone `a`-linear deficit; the interior min sits low near
`j₀∈{0,3}`. **That asymmetry — large index + carry-built valuation vs. a single linear
deficit — is why the boundary always loses, uniformly in c.**

For the prover: `Δval ≥ 2` in-gate is the live margin; you are not racing a tie. The
`c=3` residual (`v₂(N₂) ≥ v₂(P+1)−1` and the `a`-odd two-factor Lemma F) is the *only*
piece between here and a full hand proof, and the data says it holds with the same ≥2
slack (m≤200). The finite tie list above is the complete set of base cases.

---

## Files
`threerow-boundary/jobA_widesweep.py` (c≤3, m≤200, gated sweep + boundedness probe +
MN cross-check), `threerow-boundary/jobA_c45.py` (c=4,5 via Jacobi–Trudi determinant),
`threerow-boundary/boundary_forms.py` (closed-form fits), proof
`proofs/2026-06-15-boundary-lemma-threerow.md`.
