# Pre-registered per-arm predictions — written 2026-10-17 c2, BEFORE these arms finished

Theory: for a one-variable chi, Psi_chi(M) = prod_{m in M} chi(w_m), the renormalised
structure constants are regular on

  * every LONG-root divisor H_ab = {w_a w_b = 1}   iff  (COB)  chi(x)chi(1/x) = g(x)^{-1},
      g(x) = (1-q)^2/((1-q x^h)(1-q x^{-h})),  h=2 in type C;
  * every SHORT-root divisor S_m^eps = {w_m = eps}, eps^h = 1,  iff  (SR) chi(eps) = 1,
      because the short-root collision multiplier is measured to be exactly 1.

(Both subject to the degenerate case of Thm C(iii): a (lam,mu) with no M'' >= lam,mu is
regular whatever chi does.)

| arm   | chi                | (COB) | (SR)                | predicted long-root poles | predicted short-root poles |
|-------|--------------------|-------|---------------------|---------------------------|----------------------------|
| id    | 1                  | FAIL  | holds (chi=1)       | the baseline 13           | none                       |
| psi   | (1-qv^2)/(1-q)     | HOLDS | holds (psi(+-1)=1)  | NONE                      | none                       |
| x2psi | v^2 (1-qv^2)/(1-q) | HOLDS | holds ((+-1)^2=1)   | NONE                      | none                       |
| xpsi  | v (1-qv^2)/(1-q)   | HOLDS | FAILS at -1 only    | NONE                      | at w_i + 1 only            |
| x     | v                  | FAIL  | FAILS at -1 only    | the baseline 13           | at w_i + 1 only            |
| sq    | ((1-qv^2)/(1-q))^2 | FAIL  | holds               | EXACTLY the baseline 13   | none                       |
| h1    | (1-qv)/(1-q)       | FAIL  | FAILS at -1 only    | the baseline 13           | at w_i + 1 only            |
| unrm  | (1-qv^2)           | FAIL  | FAILS at BOTH +-1   | the baseline 13           | at w_i+1 AND w_i-1         |

Sharpest of these, because it is a SET prediction and not a count: **sq must reproduce
exactly the same 13 triples as id**, since Thm C(iii) makes the pole locus depend on chi
only through the predicate (kappa*g == 1), which fails for both.
And **unrm is the only arm predicted to produce a (w_i - 1) pole.**

Already in at the time of writing: id (13, long-root only), psi (0), xpsi (57, all w_i+1,
zero long-root). The five rows above for x2psi / x / sq / h1 / unrm were unmeasured.
