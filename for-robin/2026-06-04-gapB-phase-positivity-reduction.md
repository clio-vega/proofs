# Gap B (two-row d=4 fibre vanishing): a cleaner, uniform reduction — and the new gap

**Date:** 2026-06-04 (prove session)
**Proof doc:** `~/projects/proofs/2026-06-04-gapB-two-row-d4.tex` (compiles, 5 pp.)

## TL;DR
I did **not** close the endpoint inequality `D_m>0` that PROVE.md targeted. Instead I
found a route that makes that endpoint-only framing unnecessary: the whole two-row law
(**all** shapes, not just the balanced endpoint) reduces to a single **phase positivity**

> `(PP)  Q_k := Im(r_k \bar r_{k-1}) > 0`  for `1≤k≤m`, all `m≥4`.

Everything around `(PP)` is now proved rigorously; `(PP)` itself is reduced to one clean
inequality family with an identified mechanism, but is **not yet fully proved** (one
uniform estimate remains). Verified exactly / 50-digit for `m≤419`.

## The chain (all rigorous except the last bullet)
- `f(u)=u²+(1+i)u+1`, `r_l=[u^l]f^m`, two-row `G_{(a,b)}(i)=r_b−r_{b-1}`, `b≤m`.
- **Phase criterion (Lemma):** `Q_b ≠ 0 ⇒ r_b ≠ r_{b-1} ⇒ G_{(a,b)}≠0`.
  (If `r_b=r_{b-1}` then `C_b=r_b\bar r_{b-1}=|r_{b-1}|²∈ℝ`, so `Q_b=0`. Contrapositive.)
  This is *weaker* than `|r_b|≠|r_{b-1}|` (the old N-monotonicity target) but still suffices
  — and it covers every shape at once, no special endpoint.
- **Casoratian recurrence:** multiplying the r-recurrence by `\bar r_k`,
  `(k+1)C_{k+1}=(1+i)α_k N_k + β_k \bar C_k`, i.e.
  `(k+1)P_{k+1}=α_k N_k+β_k P_k` and **`(k+1)Q_{k+1}=α_k N_k − β_k Q_k`**
  (`α_k=m−k`, `β_k=2m−k+1`). This is exactly the discrete Wronskian of `(r_l, \bar r_l)`.
- **`P_k>0` for all `1≤k≤m` (proved, clean induction):** `P_1=m>0`; if `P_k>0` then
  `C_k≠0` so `r_k≠0`, `N_k>0`, and `(k+1)P_{k+1}=α_kN_k+β_kP_k>0`. (Gives `N_k>0` for free.)
- **Reduction:** since `Q_1=m>0`, `(PP) ⇔ (B_k): β_k Q_k < α_k N_k` for `1≤k≤m−1`
  (because `(B_k) ⇔ Q_{k+1}>0`). Set `x_k = β_k Q_k/(α_k N_k)`; `(B_k) ⇔ x_k<1`.
  Base: `x_1 = 1/(m−1)` (so `<1` iff `m≥3`; `=1` at `m=2` → the `(2,2)` exception). Also
  `Q_2 = m²(m−2)`.
- **Finite cases `m≤3`:** done directly. NB `(3,3)` has `Q_3=0` (phase criterion mute) but
  `G_{(3,3)}=1+2i≠0` — a genuine "harmless zero" of `Q` that is **not** a theorem exception.
  The only `x_k≥1` ever are `(m,k)=(2,1)` and `(3,2)`, recording `(2,2)` [real] and `(3,3)`
  [harmless] respectively. For `m≥4`, `x_k ≤ 0.78` (max `0.7798` at `(6,5)`).

## The remaining gap — and its mechanism (this is the interesting part)
`x_k<1` is governed by the exact 2-D system (state `(μ_k=N_k/N_{k-1}, x_k)`):
`x_{k+1} = [β_{k+1}α_k/(α_{k+1}(k+1))]·(1−x_k)/μ_{k+1}` coupled to the Gram relation for `μ`.
First eqn has the form `x_{k+1}=a_k(1−x_k)`, whose fixed point `a_k/(1+a_k)` is *always* `<1`.
The only difficulty: near the top `a_k>1` (since `α_{k+1}=m−k−1` is small), so we need the
`(1−x_k)` factor small, i.e. a matching **lower** bound on `x_k`.

In the distance-from-top variable `j=m−k`, `y_j:=x_{m−j}`, the system degenerates (as
`m→∞`, `μ→1`) to the limiting recursion **`y_{j-1} = j/(j-1)·(1−y_j)`**, which has an exact
**self-propagating invariant `½ < y_j < (j+1)/(2j)`** (Lemma in the doc) driving `y_j` from
`½` (large j) up to `y_1→¾ < 1`. So `¾` is where the endpoint `x_{m-1}` limits, explaining
the numerics. I verified the limiting recursion matches exact `x_k` near the top to ~1e-2
(shrinking toward `j=1`) at `m=300`, and the band holds for `j=1..8`.

**What's left:** a *uniform-in-m* control of the corrections in the exact system to the
limiting recursion, on the bounded top window where `a_k>1`, plus an easy bulk estimate
`x_k≤½` for large `j` (where `a_k<1` makes `x_{k+1}=a_k(1−x_k)<a_k<1` immediate). Both are
inequalities purely in `(μ_k,x_k)`; the friction is the square-root coupling through `μ`.

## Why I think this supersedes the old PROVE.md framing
The old route (`N_l` strictly increasing) needed the endpoint `D_m>0` = balanced shape
`(m,m)`, with the obstruction being the *phase* of `r_{m-1}\bar r_{m-2}`. The new route makes
that phase the **direct object** (`Q_k`), proves the surrounding scaffold unconditionally,
and reveals the universal `¾` profile. If you want to push it: the cleanest target is now
"`x_k<1` for `m≥4`" via the 2-D system — quite possibly a half-page once the right
two-sided envelope on `(μ,x)` near the top is written down. I ran out of session before
nailing the uniform square-root bound.

— Clio
