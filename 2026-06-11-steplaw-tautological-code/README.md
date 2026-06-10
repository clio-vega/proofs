# 2026-06-11 — the d=4 step law is tautological on J*-pairs (verification code)

Companion to `proofs/2026-06-11-steplaw-M-half.md`.

- `jobM_explore.py` — main exploration. Defines `D_j = 2^j M_j = <s_λ, p1^{2(m-j)}(p1^2-p2)^j>`,
  finds toggle pairs, checks the (tautological) a=3 step law `Δv2(M)=-4` on genuine J*-pairs.
- `job1_tie_census.py` — M_j / χ_b / val(j) / J* machinery (exact Murnaghan–Nakayama, no Sage).
- `symfunc.py` — pure-Python power-sum symmetric-function inner products (cross-checks Tools A,B,C).
- `core_quotient.py` — abacus d-core / d-quotient.

Key verified facts (exact):
- Tool A `D_j = <s_λ,p1^{2(m-j)}(p1^2-p2)^j> = 2^j M_j` : 159/159 (m≤6)
- Tool B `A_λ(x) = Σ_r C(m,r) R_r (x+2)^r`            : 159/159
- Tool C `M_j = Σ_μ (#vertical-2-strip chains λ→μ) f^μ`: 159/159
- Theorem 1 (step law forced by val-constancy)        : 38/38 a=3 J*-pairs, m≤12
- counterexamples to a non-circular reading           : 2470 (j∈J*, bit3=0, j+8∉J*), m≤12
- coarse engine `Φ/2^e ≡ z^{j0}(1+z)^g`, g≥2 on ties  : 1624/1624 ties, m≤12

Run: `python3 jobM_explore.py`
