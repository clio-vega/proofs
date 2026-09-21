# Step 0 — hypotheses actually consumed by proofs/2026-09-20-c1-cylindric-M-convexity.tex §3–§6
*Line numbers refer to that file. Written BEFORE any porting.*

The paper's §8 ("What the proof consumes") lists external *citations*. That is NOT the same
list as the *structural* hypotheses of §3–§6. Here is the structural list.

## The abstract skeleton

Data: a finite set `T` of "chains" `x = (x^0,...,x^ell)` with `x^0, x^ell` FIXED; a scalar
statistic `u` on states; weight `alpha_t(x) = u(x^t) - u(x^{t-1})`; `W = {alpha(x) : x in T}`.

### (H1) Telescoping weight. [eq:weight, used in prop:max l.505 and throughout]
`alpha` is the increment sequence of a single scalar `u`. Consequences used:
`sum_t alpha_t = u(x^ell) - u(x^0) = d` is constant (homogeneity), and
`alpha_1+...+alpha_r = u(x^r) - u(x^0)` (partial sums are prefix values of `u`).

### (H2) Markov box property. [lem:box l.311]
With `p = x^{t-1}` and `r = x^{t+1}` held fixed, the admissible `x^t` form a product of
integer intervals `prod_i [L_i, U_i]`, AND `u` restricted to that product is the SUM of the
independently-varying coordinates. Consequence: `alpha_t` sweeps the FULL integer interval
`[mu_-, mu_+]` with `mu_- = sum L_i - u(p)`, `mu_+ = sum U_i - u(p)`.

### (H3) PALINDROME / self-duality of the box.  <-- THE HIDDEN THIRD HYPOTHESIS
[lem:palindrome l.330, restated as eq:pal l.395]
`mu_- + mu_+ = N` where `N = alpha_t + alpha_{t+1}`.
**This does NOT follow from (H1)+(H2).** (H2) says the fibre is a full interval; (H3) says
that interval is symmetric about its own midpoint `N/2`. Verified consumed in two places:
  * `cor:bk` (l.400) — the Bender-Knuth involution `x_i -> L_i+U_i-x_i` swaps `alpha_t` and
    `alpha_{t+1}` ONLY because of the equality. This is what makes `W` `S_ell`-stable.
  * `cor:exchange` (l.416) — needs only the INEQUALITY `mu_- + mu_+ <= N`, since together
    with `mu_- <= mu_+` that yields `2 mu_- <= N`, hence `alpha_t > alpha_{t+1} => alpha_t > mu_-`.
And the paper's own l.358 remark flags it: the `+n` wrap-around term is bead-specific.

### (H4) A termwise-maximal chain exists. [lem:greedy l.472(3), prop:max l.505]
There is `G in T` with `x^t_i <= g^t_i` for all `x in T`, all `t`, all `i` — hence
`u(x^t) <= u(g^t)`. Only `u(x^t) <= u(g^t)` is used, i.e. a chain maximising EVERY prefix
of `u` simultaneously. This is a greedy/lattice property of the arena, independent of (H1)-(H3).

### (H5) Nothing else.
lem:hlp (Robin Hood, l.421) is pure partition combinatorics. §6 is pure polymatroid theory
from the exchange axiom. `prop:dictionary` (§2) is a translation, not a hypothesis.

## Verdict
The template has **four** hypotheses, not two. The brief's Step-0 suspicion was correct:
**(H3) is real and was unnamed.** Porting (H1)+(H2) alone would be a shape match.

## Restated as a self-contained generic theorem (to be proved once, arena-free)
> Let `W subset Z^ell_{>=0}` be finite, and suppose:
>  (A) `sum_t alpha_t = d` for all `alpha in W`;
>  (B) for every `alpha in W` and every `1 <= t <= ell-1` there are integers `m <= M` with
>      `m + M = alpha_t + alpha_{t+1}` such that `{beta_t : beta in W, beta_s = alpha_s (s != t,t+1),
>      beta_t + beta_{t+1} = alpha_t + alpha_{t+1}} = [m, M]` and `alpha_t in [m,M]`;
>  (C) there is `alpha^* in W` with `alpha_1 + ... + alpha_r <= alpha^*_1 + ... + alpha^*_r`
>      for all `alpha in W` and all `r`.
> Then `W` is M-convex, `P(W) = {sigma : sigma <| sort(alpha^*)}`, and
> `Newton = P_{sort(alpha^*)}`.
(B) packages (H2)+(H3); (C) packages (H1)+(H4). This is the thing to port.
