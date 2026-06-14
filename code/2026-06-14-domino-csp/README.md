# 2026-06-14 domino-csp — does the order law live on the CTT domino-tableau CSP?

Companion to connection `oddcontent-falsified-domino-csp-survives`. **Verdict: DEAD.**

The Colmenarejo–Tenner–Thompson domino CSP (arXiv 2602.23343) reads the **2-quotient**
(`|SDT(λ)| = C(m,|q0|)·f^{q0}·f^{q1}`, 0 failures); the order law `ord_{q=−1}G_λ=⌊|2-core|/2⌋`
reads the **2-core**. Complementary halves of the partition — orthogonal by construction.
On CTT's whole domain (domino-tileable shapes) the order law is identically 0; the discriminating
staircases `δ_k` (their own 2-cores) carry no domino tableaux at all.

- `2026-06-14-jobA-domino-csp.py` — A1 CTT CSP grounded on binary words; A2 shape+maj fit
  (`(N,N)` two-row rectangle, maj convention `top_gt_sum_i`); A3 tileable ⟺ order-law-trivial;
  A4 staircases outside the domain; A5 `|SDT|` factors through the 2-quotient.
- `FINDINGS.md` — full verdict and the redirect (core-reading domino statistic = spin/cospin,
  i.e. the LLT / Littlewood-t=2 face, not CTT's maj).

Pure Python (sympy); self-checking.
