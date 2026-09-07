# Q96 verification code (2026-09-07 c2)

Three engines, code-disjoint where it matters.

- `engine.py`  — engine A (abacus / Maya-set legal e-moves) and engine B
                 (brute-force border strips on Young diagrams: 2x2 test +
                 edge-connectivity flood fill).  Cross-checked 60/60.
- `pring.py`   — engine C: Schur functions via Jacobi-Trudi as polynomials in
                 p_1,p_2,...; multiplication by the SYMBOL p_e; re-expansion in
                 the Schur basis by a linear solve over the p-monomial basis.
                 Shares no code with A or B.  Used to verify the premise
                 M_{p_e} = R_e(-1) and, independently, the main witness.
- `nest.py`    — nested commutators ad(M_{p_e})^d applied to R_e(t).
- `twoparam.py`— the two-parameter commutator [R_e(t), R_f(s)] and the predicted
                 matrix elements of Theorem 4.1.  1829/1829.

NOTE on `pring._vec`: it must make the p_k the Poly generators and keep the FULL
coefficient.  An earlier version used `as_coefficients_dict`, which returns only
the numeric part and silently dropped the symbolic t — reporting every value at
t=1.  See the "instrument failure" remark in the paper.
