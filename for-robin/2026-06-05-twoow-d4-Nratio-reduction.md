# Two-row d=4 fibre law now hangs on ONE norm-ratio inequality (2026-06-05)

**Setup.** `f(u) = u^2 + (1+i)u + 1`, `r_l = [u^l] f^m`, `N_l = |r_l|^2`,
`C_l = r_l * conj(r_{l-1}) = P_l + i Q_l`. For a two-row shape `(a,b) ⊢ 2m`
the d=4 fibre value is `G_{(a,b)}(i) = r_b − r_{b-1}`.

**The law.** `G_{(a,b)}(i) = 0 ⟺ (a,b) = (2,2)`.

## What advanced today

The whole law is equivalent to the phase positivity `Q_l > 0` for
`1 ≤ l ≤ m−1`, i.e. `x_l := β_l Q_l / (α_l N_l) < 1` (with `α_l = m−l`,
`β_l = 2m−l+1`). Last session this was framed as a *downward-in-j* band that
did not self-propagate. I replaced it with a **forward induction in l** driven
by an exact **master identity**:

    E / N_l = (j−1) B^2 (s_l + j/B)^2 + j (x_l − 1) f_2(x_l),
    j = m−l,  B = m+j+1,  s_l = P_l/N_l,  f_2 = m(m+1) − 2j(j−1) + j(j−1)x_l > 0,
    and  E > 0  ⟺  x_{l+1} < 1.

Since `x_l < 1` (the inductive hypothesis), the second term is ≤ 0, so each
step reduces to a one-line inequality `(PP_l)`. Together with the
Cauchy–Schwarz **equality** `P_l^2 + Q_l^2 = N_l N_{l-1}` (a lower bound on
`s_l` becomes an upper bound on `μ_l = N_l/N_{l-1}`), the induction closes in
three regimes: bulk (`a_l ≤ 1`, unconditional), entry strip
(`1 < a_l ≤ 3/2`, needs `x_l > 1/3`), deep band (`a_l > 3/2`, needs the master
identity + the norm-ratio bound + `x_l ≥ 1/2`).

I also nailed the **sharp limiting profile** `ℓ_j = (2j+1)/(4j)` (the exact
solution of `ℓ_{j−1} = (j/(j−1))(1−ℓ_j)`, with `x_{m−j} ↗ ℓ_j` from below),
and **proved the m→∞ law** `E(j) = (4j²−1)/(16j) > 0` for all j.

## The single remaining gap (finite / algebraic)

Everything reduces to two exactly-verified finite inequalities. The soft one
(`x_l > 1/3` on the entry strip, `x_l ≥ 1/2` on the top — both follow from
`x_l` increasing in l) is essentially boundary bookkeeping. **The load-bearing
one is:**

> **Lemma 1 (norm-ratio bound).** For `m ≥ 4`, `1 ≤ j ≤ m−1` with `2m−3j > 0`,
> and `l = m−j`,
>
>     μ_l = N_l / N_{l-1}  ≤  1 + 3(2j+1) / (2m − 3j).

Verified exactly: **1 violation in 119798** cases over `m ≤ 600` — the lone
exception `(m,j) = (5,1)` is off-band. Asymptotically tight:
`μ_l − 1 ~ 3(2j+1)/(2m)`. The full induction certificate (including the closing
rational inequality) passes with 0 failures to `m = 700`.

Caveat I want to flag: `N_l = Norm_{Q(i)/Q}(r_l)` is the norm of an order-2
holonomic sequence, so `{N_l}` is order-4 holonomic (symmetric-square type).
There is therefore **no second-order Riccati** for `μ_l`, and no single
discriminant gives monotonicity — Lemma 1 has to be proved as the genuine
norm-ratio envelope it is.

## The ask

**Do you know a slick proof of a norm-ratio upper bound for
`|r_l|^2` where `r_l = [u^l] f^m`?** Concretely: `N_l/N_{l-1} ≤ 1 + O((2j+1)/m)`
with `j = m−l`, for the specific `f = u^2 + (1+i)u + 1`. The N-recurrence is
`(l+1)^2 N_{l+1} = 2(m−l)^2 N_l + (2m−l+1)^2 N_{l-1} + 2(m−l)(2m−l+1)(P_l−Q_l)`
with `P_l^2 + Q_l^2 = N_l N_{l-1}`. Any clean envelope of this shape closes the
entire two-row d=4 law.

Full writeup (6 pp, compiles): `2026-06-05-gapB-two-row-Nratio.tex` in the
proofs repo. It supersedes `2026-06-04-gapB-two-row-d4.tex` — note that the old
"bounded window in j" remark there is now **false** (the band window grows like
`√m`).
