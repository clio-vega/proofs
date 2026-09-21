# code-q209 — affine Stanley symmetric polynomials, direct from (def-affine-stanley)

Companion code for `../2026-09-21-c1-affine-stanley-exchange.tex`.

- `affstan.py`   — affine symmetric group in window notation; cyclically decreasing elements;
                   cyclically decreasing factorisations. Implements WZZ (def-affine-stanley)
                   directly, by dynamic programming over chains in the affine weak order.
                   No theorem from the paper is used anywhere in it.
- `verify_all.py`— reproduces every number quoted in section 7 of the paper in one run.
- `generic.py`   — exhaustive search over homogeneous subsets of the simplex; the source of
                   the counterexample in Proposition 5.5 ((H4) is independent).
- `STEP0-hypotheses.md` — the hypothesis extraction from the cylindric paper that this
                   session began with; (H3) and (H4) are named there for the first time.

    python3 verify_all.py
